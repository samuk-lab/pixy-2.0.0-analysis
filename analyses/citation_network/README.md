# citation_network → Figure S1

Pulls citing-paper metadata for the original pixy paper (Korunes & Samuk 2021,
doi 10.1111/1755-0998.13326) from the OpenAlex API, and reduces it to the two tables
that Figure S1 plots.

## Overview

`build_pixy_citation_data.R` looks up the pixy DOI on OpenAlex, pages through every
work that cites it, and tallies the citing papers by year and by OpenAlex subfield.
It writes, into `pixy_citer_topic_network_output/`:

- `pixy_citations_per_year.csv` — citations added per year, and the cumulative total
- `subject_counts.csv` — citing papers per subfield

Both feed `../../figures/Figure_S1_Citations.R`.

## Scripts

| Script | Purpose |
|--------|---------|
| `build_pixy_citation_data.R` | Query OpenAlex for works citing the pixy DOI and write the two CSVs above |

## Run

```sh
Rscript build_pixy_citation_data.R
```

The script needs internet access (it queries `api.openalex.org`). No key is required;
setting `OPENALEX_EMAIL` in the environment joins OpenAlex's polite pool for better
rate limits. Other options — the focal DOI, subject level, and record cap — are
environment variables documented in the script header. Because the citation count
grows over time, a rerun gives slightly different numbers than the shipped CSVs, which
back the published Figure S1.

## Dependencies

R packages: `httr2`, `jsonlite`, `dplyr`, `purrr`, `tibble`, `tidyr`, `stringr`,
`readr` (installed automatically unless `INSTALL_MISSING_PACKAGES=FALSE`).
