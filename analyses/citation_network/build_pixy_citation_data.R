#!/usr/bin/env Rscript

# citing-paper metadata for korunes & samuk 2021 pixy (doi 10.1111/1755-0998.13326)
# writes subject_counts.csv (subject_area, n_citing_papers) and
# pixy_citations_per_year.csv (year, citations_added, cumulative), feeds figures/Figure_S1_Citations.R
#
# env vars:
#   OPENALEX_API_KEY, OPENALEX_EMAIL  polite pool / rate limits
#   FOCAL_DOI         default 10.1111/1755-0998.13326
#   OUT_DIR           default pixy_citer_topic_network_output
#   MAX_CITERS        default Inf (number to cap)
#   SUBJECT_LEVEL     domain, field, subfield, topic (default subfield)
#   INSTALL_MISSING_PACKAGES   default TRUE

##########
# configuration
##########
# OPENALEX_EMAIL joins the polite pool, passed as mailto. optional
# https://docs.openalex.org/how-to-use-the-api/rate-limits-and-authentication

env_chr <- function(name, default) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) default else x
}

env_flag <- function(name, default = FALSE) {
  x <- toupper(env_chr(name, ifelse(default, "TRUE", "FALSE")))
  x %in% c("TRUE", "T", "1", "YES", "Y")
}

env_num_or_inf <- function(name, default) {
  x <- env_chr(name, as.character(default))
  if (tolower(x) %in% c("inf", "infinite", "all")) return(Inf)
  y <- suppressWarnings(as.numeric(x))
  if (is.na(y)) default else y
}

target_doi <- env_chr("FOCAL_DOI", "10.1111/1755-0998.13326")
out_dir <- env_chr("OUT_DIR", "pixy_citer_topic_network_output")
max_citers <- env_num_or_inf("MAX_CITERS", Inf)
subject_level <- tolower(env_chr("SUBJECT_LEVEL", "subfield"))
install_missing_packages <- env_flag("INSTALL_MISSING_PACKAGES", TRUE)

allowed_subject_levels <- c("domain", "field", "subfield", "topic")
if (!subject_level %in% allowed_subject_levels) {
  stop(
    "SUBJECT_LEVEL must be one of: ", paste(allowed_subject_levels, collapse = ", "),
    call. = FALSE
  )
}

request_sleep_seconds <- ifelse(nzchar(Sys.getenv("OPENALEX_API_KEY")), 0.05, 0.25)

##########
# package setup
##########
required_packages <- c(
  "httr2", "jsonlite", "dplyr", "purrr", "tibble", "tidyr",
  "stringr", "readr"
)

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  if (install_missing_packages) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, repos = "https://cloud.r-project.org")
  } else {
    stop(
      "Missing packages: ", paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
}

invisible(lapply(required_packages, library, character.only = TRUE))

##########
# helpers
##########
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

clean_doi <- function(x) {
  x |>
    stringr::str_trim() |>
    stringr::str_remove("^https?://(dx\\.)?doi\\.org/") |>
    stringr::str_remove("^doi:")
}

short_oa_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_remove(x, "^https?://openalex\\.org/")
  x[is.na(x) | x == "NULL" | x == ""] <- NA_character_
  x
}

as_chr1 <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0) return(default)
  y <- suppressWarnings(as.character(x[[1]]))
  if (is.na(y) || y == "NULL") default else y
}

as_int1 <- function(x, default = NA_integer_) {
  if (is.null(x) || length(x) == 0) return(default)
  y <- suppressWarnings(as.integer(x[[1]]))
  if (is.na(y)) default else y
}

primary_subject_area <- function(work) {
  pt <- work$primary_topic %||% list()
  sa <- switch(
    subject_level,
    domain   = as_chr1(pt$domain$display_name,   NA_character_),
    field    = as_chr1(pt$field$display_name,    NA_character_),
    subfield = as_chr1(pt$subfield$display_name, NA_character_),
    topic    = as_chr1(pt$display_name,          NA_character_)
  )
  if (is.na(sa) || sa == "") "Unclassified" else sa
}

##########
# openalex api
##########
openalex_base <- "https://api.openalex.org"
openalex_api_key <- Sys.getenv("OPENALEX_API_KEY")
openalex_email <- Sys.getenv("OPENALEX_EMAIL")

work_select <- paste(
  c("id", "doi", "publication_year", "primary_topic"),
  collapse = ","
)

oa_get <- function(path, query = list()) {
  Sys.sleep(request_sleep_seconds)

  if (nzchar(openalex_api_key)) query$api_key <- openalex_api_key
  if (nzchar(openalex_email)) query$mailto <- openalex_email

  req <- httr2::request(paste0(openalex_base, path))
  req <- httr2::req_user_agent(
    req,
    paste0(
      "pixy-citation-data-rscript/0.3 ",
      ifelse(nzchar(openalex_email), paste0("(", openalex_email, ")"), "")
    )
  )
  if (length(query) > 0) {
    req <- do.call(httr2::req_url_query, c(list(req), query))
  }
  req <- httr2::req_retry(req, max_tries = 4)
  req <- httr2::req_error(req, is_error = function(resp) FALSE)

  resp <- httr2::req_perform(req)
  status <- httr2::resp_status(resp)
  if (status >= 400) {
    body_txt <- tryCatch(httr2::resp_body_string(resp), error = function(e) "")
    stop(
      "OpenAlex request failed with HTTP ", status, ".\n",
      "If you see 403/409/rate-limit, set OPENALEX_API_KEY before running.\n",
      body_txt,
      call. = FALSE
    )
  }

  httr2::resp_body_json(resp, simplifyVector = FALSE)
}

fetch_work_by_doi <- function(doi) {
  doi <- clean_doi(doi)
  message("Looking up focal paper by DOI: ", doi)
  res <- oa_get(
    "/works",
    list(
      filter = paste0("doi:https://doi.org/", doi),
      per_page = 1,
      select = "id,doi,publication_year,cited_by_count"
    )
  )
  if (length(res$results) == 0) {
    stop("No OpenAlex work found for DOI: ", doi, call. = FALSE)
  }
  res$results[[1]]
}

fetch_citing_works <- function(focal_id, max_records = Inf) {
  focal_id <- short_oa_id(focal_id)
  cursor <- "*"
  out <- list()
  expected_count <- NA_integer_

  repeat {
    remaining <- if (is.infinite(max_records)) 200 else max_records - length(out)
    per_page <- min(200, remaining)
    if (per_page <= 0) break

    message("Fetching works that cite ", focal_id, " ... currently have ", length(out))
    res <- oa_get(
      "/works",
      list(
        filter = paste0("cites:", focal_id, ",from_publication_date:2021-01-01"),
        sort = "publication_date:desc",
        per_page = per_page,
        cursor = cursor,
        select = work_select
      )
    )

    if (!is.null(res$meta$count)) expected_count <- as.integer(res$meta$count)
    if (length(res$results) == 0) break
    out <- c(out, res$results)

    if (!is.na(expected_count)) {
      message("OpenAlex reports ", expected_count, " total citing works; downloaded ", length(out), ".")
    }

    next_cursor <- res$meta$next_cursor %||% NA_character_
    if (is.na(next_cursor) || next_cursor == cursor) break
    cursor <- next_cursor

    if (!is.infinite(max_records) && length(out) >= max_records) break
  }

  out
}

##########
# build outputs
##########
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

focal_work <- fetch_work_by_doi(target_doi)
focal_id <- short_oa_id(focal_work$id)

message("Focal OpenAlex ID: ", focal_id)
message("OpenAlex cited_by_count: ", focal_work$cited_by_count %||% NA_integer_)
message("Subject level: ", subject_level)

citing_works <- fetch_citing_works(focal_id, max_records = max_citers)

citer_tbl <- tibble::tibble(
  subject_area = purrr::map_chr(citing_works, primary_subject_area),
  publication_year = purrr::map_int(citing_works, ~ as_int1(.x$publication_year))
)

# subject_counts.csv
subject_counts <- citer_tbl |>
  dplyr::count(subject_area, sort = TRUE, name = "n_citing_papers")

readr::write_csv(subject_counts, file.path(out_dir, "subject_counts.csv"))

# pixy_citations_per_year.csv
per_year <- citer_tbl |>
  dplyr::filter(!is.na(publication_year)) |>
  dplyr::count(publication_year, name = "citations_added") |>
  dplyr::arrange(publication_year) |>
  dplyr::rename(year = publication_year) |>
  dplyr::mutate(cumulative = cumsum(citations_added))

readr::write_csv(per_year, file.path(out_dir, "pixy_citations_per_year.csv"))

message("Done. Wrote outputs to: ", normalizePath(out_dir))
message("Top subject areas:")
print(subject_counts |> dplyr::slice_head(n = 20), n = 20)
message("Citations per year:")
print(per_year, n = nrow(per_year))
