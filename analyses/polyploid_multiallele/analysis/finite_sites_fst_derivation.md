# Finite-sites expectation for Hudson's FST under the two-population split model

**Date:** 2026-07-14
**Purpose:** supply the θ-dependent expectation that the *multiallelic* Hudson FST should
be compared against in Figure 4. The existing `theoretical_fst` in
`11_dxy_fst_theta_sweep.R` is the infinite-sites (coalescent) value and is the correct
target for the *biallelic* estimator only.

**Status:** derived, and validated against four independent checks (§5) — check 4 now
against the full n = 2000 rerun rather than a single smoke-test replicate. Supersedes the
approximate `g(D) = D/(1 + 4D/3)` treatment in `HANDOFF_fst_multiallelic_2026-07-14.md`,
which predicted a −19% effect at θ = 0.1. The correct value is −1.6%.

**Update 2026-07-15.** §5.4 has been restated: its original ratio test passes at n = 1
and fails at n = 2000 for a reason that has nothing to do with this derivation, which the
rerun confirms. If you are here to check whether the finite-sites expectation holds, read
§5.4.1 — it does, within noise, at every θ × ploidy. Do not resurrect the ratio test.

---

## 1. Model and notation

The simulated model (vcfsim, `dxy_*_2pop_theta*` arms): a single ancestral population of
size `alpha * Ne` splits `T` generations ago into two populations each of size `Ne`, with
no subsequent migration. Mutation follows Jukes & Cantor (1969): every site mutates to
each of the three alternative nucleotides at equal rate, with total substitution rate `mu`
per site per generation.

| symbol | meaning |
|---|---|
| `Ne` | 1,720,600 (per-population effective size) |
| `T` | 1e6 generations (split time) |
| `alpha` | 2 (ancestral population is `alpha * Ne`; vcfsim doubles it) |
| `tau` | `T / (2 * Ne)` = 0.290596 (split time in coalescent units) |
| `theta` | `4 * Ne * mu` |
| `mu` | total substitution rate per site per generation |

Note `theta` is defined with the *per-population* `Ne`, matching the sweep's
`theta_nom = 4 * Ne * mu`.

---

## 2. The master equation

> **Prior art — read this before claiming novelty (added 2026-07-15).** An earlier version of
> this document presented the Laplace-transform step below as "the whole trick", which
> implied it was new here. It is not. Deriving F_ST from the distribution of pairwise
> coalescence times, and finding that it depends on the mutation rate, is the standard
> framework:
>
> - **Rousset (1996)** relates the probability of identity *in state* to the distribution of
>   coalescence times and computes F_IS and F_ST from it, explaining their behaviour "at high
>   mutation rate loci" (stepwise mutation, island model).
> - **Wilkinson-Herbots (1998)** calculates "the Laplace transform of the distribution of the
>   coalescence time of a pair of genes" and uses it to obtain F_ST, its θ → 0 limit, *and*
>   its derivative with respect to the mutation rate, concluding that F_ST "can depend
>   strongly on the mutation rate" (infinite alleles; island and stepping-stone models). §2
>   and the §5.3 check below are both instances of this.
> - **Charlesworth (1998)** states the qualitative core for DNA sequence data: F_ST "is
>   strongly influenced by the level of within-population diversity".
> - **de Jong et al. (2024)** is the closest to our setting: it advocates Hudson's
>   `F_ST = (dxy − pi_xy)/dxy` computed on *observed sequence dissimilarity*, and notes that
>   F_ST maps onto coalescent units only "in case of recent population splits, when novel
>   mutations are negligible", with "no universal relationship" otherwise. That is the
>   two-estimand result of this document, stated qualitatively.
>
> What is actually ours is the instance, not the method: the closed form for a clean
> two-population split (no migration) with a doubled ancestral population, under Jukes–Cantor
> (a k = 4 allele model) rather than infinite-alleles or stepwise mutation, evaluated for
> Hudson's ratio-of-averages estimator. Its value is that it turns de Jong et al.'s
> qualitative condition into the specific number Figure 4 is judged against (−1.6% at
> θ = 0.1). Frame it that way in the manuscript — "following Rousset (1996) and
> Wilkinson-Herbots (1998)" — not as a new derivation.
>
> Searched 2026-07-15 (PubMed, Consensus, web). Island and stepping-stone models are covered
> thoroughly by the above; a closed form for the *no-migration split model under JC69*
> specifically was not found, but absence of a hit is weak evidence — this is a 50-year-old
> literature and much of it is in books (e.g. Wakeley's *Coalescent Theory*) that these
> searches do not index. Do not assert novelty in print without a librarian-grade check.

Everything below follows from one observation.

Under Jukes–Cantor, the probability that two sequences separated by a **total path length**
`t` differ at a given site is

```
P(differ | t) = (3/4) * (1 - exp(-(4/3) * mu * t))            (Jukes & Cantor 1969)
```

Two lineages whose most recent common ancestor lived `T_c` generations ago are separated by
a path of length `t = 2 * T_c` (one branch each). Substituting, and writing
`beta = (8/3) * mu`:

```
P(differ | T_c) = (3/4) * (1 - exp(-beta * T_c))
```

Taking expectations over the coalescence-time distribution and using the law of total
expectation gives the **master equation**:

```
E[pi] = (3/4) * (1 - L(beta))         where   L(beta) = E[exp(-beta * T_c)]
```

`L` is the Laplace transform of the coalescence-time distribution, evaluated at `beta` —
the same object Wilkinson-Herbots (1998) computes for the island and stepping-stone models.

This is the whole trick (standard, not new — see the prior-art note above). Every
estimand — π within a population, dxy between populations — is the *same* function of
`beta`; they differ only in the distribution of `T_c`. Infinite-sites results are the `beta -> 0` linearisation of this same expression
(`1 - exp(-beta*T_c) ≈ beta*T_c`, giving `E[pi] ≈ (3/4)*beta*E[T_c]`, i.e. proportional to
mean coalescence time). Finite-sites corrections are *not* an extra factor bolted on; they
are what you get by not linearising.

A useful identity, with `a = 1/(2*Ne)` the pairwise coalescence rate within a population:

```
beta / a = (8/3) * mu * 2 * Ne = (16/3) * Ne * mu = (4/3) * theta
```

so `beta = a * (4/3) * theta`. Every `theta` below enters through this one relation.

---

## 3. The two coalescence-time distributions

### 3.1 dxy — lineages in different populations

They cannot coalesce before the split. `T_c = T + T_a`, where `T_a` is exponential with
rate `1/(2 * alpha * Ne) = a/alpha` (coalescence in the ancestral population).

```
L_dxy = exp(-beta*T) * (a/alpha) / ((a/alpha) + beta)
      = exp(-(4/3)*theta*tau) / (1 + alpha*(4/3)*theta)

E[dxy] = (3/4) * (1 - exp(-(4/3)*theta*tau) / (1 + alpha*(4/3)*theta))
```

With `alpha = 2` this is exactly the script's existing `theoretical_dxy`:

```r
theoretical_dxy = (3/4) * (1 - exp(-8 * mu * split_time / 3) / (1 + (8/3) * theta_nom))
```

This answers the standing question about the `8/3`: it is `alpha * (4/3)` with
`alpha = 2`. The doubled ancestral population is why the coefficient is 8/3
rather than 4/3, as the script's comment asserts. The ancestral coalescence rate
is halved, so the ancestral polymorphism term saturates at twice the θ scale. Reproducing
the script's formula from first principles anchors the rest of the derivation:
the same machinery, applied to a different `T_c`, gives `E[pi_w]`.

### 3.2 pi_within — lineages in the same population

Now coalescence *can* happen before the split, at rate `a`. Splitting on whether it does:

```
L_pi_w = a*(1 - exp(-(a+beta)*T))/(a+beta)          [coalesced before T]
       + exp(-(a+beta)*T) * (a/alpha)/((a/alpha)+beta)   [survived to the ancestral pop]
```

In θ and τ, writing `k = 1 + (4/3)*theta`:

```
L_pi_w = (1 - exp(-tau*k))/k + exp(-tau*k) / (1 + alpha*(4/3)*theta)

E[pi_w] = (3/4) * (1 - L_pi_w)
```

### 3.3 The result

```
E[FST]_finite = 1 - E[pi_w] / E[dxy]
```

This is Hudson's FST as a ratio of expected dissimilarities, which is what pixy's
ratio-of-averages estimator targets (Hudson, Slatkin & Maddison 1992; Bhatia et al. 2013).

---

## 4. Numerical values

```
  theta   E[FST] infinite   E[FST] finite   rel. diff
  0.005       0.23696          0.23676        -0.08%
  0.010       0.23696          0.23657        -0.16%
  0.025       0.23696          0.23599        -0.41%
  0.050       0.23696          0.23502        -0.82%
  0.100       0.23696          0.23309        -1.63%
```

The finite-sites effect on FST is small — 1.6% at θ = 0.1, not the ~19% the handoff's
approximation predicted. The intuition: finite-sites saturation compresses `pi_w` and `dxy`
*in nearly the same proportion*, and FST is their ratio, so most of the effect cancels.
Only the residual asymmetry survives, because `dxy` (with its longer expected path) sits
further up the saturation curve than `pi_w` and so is compressed slightly more.

Why the handoff's approximation failed: `g(D) = D/(1 + 4D/3)` is the correct saturation map
for a **single panmictic population**, where `T_c` is exponential (it is exactly §5.1
below). Applying it to `dxy` assumes `T_c` is exponential there too — but between
populations `T_c` is *shifted* by the deterministic split time `T`, and a shifted
exponential saturates far less than a plain one. The tell was already in the handoff: its
own dxy came out 2.7% off the script's value. The derivation above reproduces the script's
dxy exactly.

---

## 5. Validation

All four checks pass (`verify_all.py`; reproduced in §7).

### 5.1 The panmictic limit reproduces Tajima (1996)

As `tau -> inf` (no split; single population) `L_pi_w -> 1/k` and

```
E[pi_w] -> (3/4)*(1 - 1/(1 + (4/3)*theta)) = theta / (1 + 4*theta/3)
```

This is **exactly** the result on p. 1458 of Tajima (1996), obtained there by substituting
`n = 2` into his Equation 6, and which he notes agrees with Equation 15 of Tajima (1983).
Against his Table 1 (Jukes–Cantor without rate variation):

| θ | Tajima (1996) Table 1 | this derivation |
|---|---|---|
| 0.005 | 0.0050 | 0.004967 |
| 0.01 | 0.0099 | 0.009868 |
| 0.02 | 0.0195 | 0.019481 |
| 0.05 | 0.0469 | 0.046875 |
| 0.1 | 0.0882 | 0.088235 |

Tajima assumes a panmictic population at equilibrium, which is precisely the `tau -> inf`
limit of the structured model — so his result is a special case of §3.2, and it matches to
every published digit.

### 5.2 E[dxy] reproduces the existing script exactly

At all five θ, to < 1e-12. See §3.1.

### 5.3 The θ → 0 limit reproduces the coalescent FST

`E[FST]_finite -> 0.236959581` as θ → 1e-9, against the infinite-sites
`(tau + (alpha-1)*(1 - exp(-tau)))/(tau + alpha) = 0.236959605`.

Equivalently, `E[pi_w] -> theta * T_W/(2Ne)` and `E[dxy] -> theta * T_B/(2Ne)` with
`T_W/(2Ne) = 1 + (alpha-1)*exp(-tau)` and `T_B/(2Ne) = tau + alpha`, recovering
`FST = 1 - T_W/T_B` — Slatkin's (1991) expression of FST as a ratio of mean coalescence
times, which is the form the script already uses.

### 5.4 It predicts the observed multiallelic FST

**Restated 2026-07-15, after the n = 2000 rerun landed. The original form of this check was
a ratio test; it does not survive at n = 2000, and it should not. It is superseded by the
two separate comparisons below, which pass. Do not reinstate the ratio test.**

#### What the check originally said, and why it was the wrong test

The first version compared, at θ = 0.100 on the n = 1 smoke test:

```
predicted  E[FST]_finite / E[FST]_infinite = 0.983689
observed   multiallelic  / biallelic       = 0.983429   |diff| = 0.00026
```

and read the 0.026-pp agreement as confirmation. The stated reasoning — that a ratio
cancels genealogy noise, because both estimators see the same tree — is true. The error is
in the *denominator*: comparing `multi/bi` against `E_fin/E_inf` silently assumes
`bi = E_inf` **exactly**. It does not. The biallelic estimator sags below the infinite-sites
value at high θ (§5.4.2). At n = 1 that sag (~0.2%) was buried in replicate noise and the
test passed by luck. At n = 2000 the noise is gone and the same test reports z = 17 to 33
at θ = 0.100 — an apparent catastrophic failure of a derivation that is in fact correct.

The discrepancy is *arithmetically exactly* the sag. Diploid, θ = 0.100:

```
obs / pred                 = 0.985195 / 0.983684 = 1.001536
(multi/E_fin) / (bi/E_inf) = 0.999897 / 0.998362 = 1.001537
```

The ratio test conflates two quantities and cannot distinguish "the finite-sites
expectation is wrong" from "the biallelic estimator is slightly biased". Test them apart.

#### 5.4.1 Multiallelic vs E[FST]_finite — PASSES

Pooled ratio-of-sums, 2000 replicates × 10 windows per cell; bootstrap SE over replicates
(400 draws, replicate = resampling unit). `multi - E_finite`:

| ploidy | θ = 0.005 | θ = 0.050 | θ = 0.100 |
|---|---|---|---|
| diploid | +0.000089 ± 0.000520 | -0.000231 ± 0.000161 | **-0.000030 ± 0.000114** |
| tetraploid | -0.000991 ± 0.000452 | -0.000023 ± 0.000159 | **+0.000037 ± 0.000105** |
| hexaploid | -0.000755 ± 0.000455 | +0.000145 ± 0.000145 | **-0.000166 ± 0.000103** |
| octoploid | -0.000266 ± 0.000432 | -0.000007 ± 0.000147 | **-0.000007 ± 0.000100** |

No trend in θ. At θ = 0.100 — where the finite-sites correction is largest and the SE is
smallest, i.e. exactly where a wrong derivation would show up — every ploidy is within
~1.6 SE of zero. The worst of all 20 cells is 2.2 SE (tetraploid, θ = 0.005), which is what
20 cells of replicate noise look like. The derivation is confirmed out-of-sample at
n = 2000, on the stronger evidence the ratio test was only pretending to give.

#### 5.4.2 Biallelic vs E[FST]_infinite — a small, real, θ-dependent sag

`bi - E_infinite` at θ = 0.100:

| ploidy | deviation | z |
|---|---|---|
| diploid | -0.000394 ± 0.000121 | -3.3 |
| tetraploid | -0.000626 ± 0.000109 | -5.7 |
| hexaploid | -0.000947 ± 0.000110 | -8.6 |
| octoploid | -0.000837 ± 0.000107 | -7.8 |

This is the residual homoplasy the handoff already anticipated (-0.16% at θ = 0.1), now
measured with error bars and significant. Biallelic filtering removes only *visible*
homoplasy; back and parallel mutation keep a site biallelic and survive the filter.

**New, and worth following up:** the sag grows with ploidy (diploid -0.00039 → hexaploid
-0.00095). More sampled chromosomes make more multi-hit sites *visible*, so more are
filtered out and the surviving biallelic set is more strongly conditioned. The biallelic
filter's effectiveness is therefore sample-size dependent, not only mutation-model
dependent. That sharpens the open non-JC69 question: "biallelic FST ≈ coalescent FST" is
contingent on both, and this sweep can already see one of the two.

#### Summary

Each estimator is unbiased for its own estimand, to different precision. The multiallelic
estimator sits on E[FST]_finite within noise everywhere. The biallelic estimator sits on
E[FST]_infinite to within 0.2–0.4%, with a significant sag at high θ. Both statements are
now reported in `figures/figs/inline_stats.md` §5.

Every number in §5.4 is reproduced by **`analysis/validate_fst_rerun.py`** (run it from
`analysis/`). It does not implement the ratio test, on purpose.

---

## 6. What this means for the paper

1. **Carry both expectations — in the analysis.** Biallelic → infinite-sites (flat 0.237);
   multiallelic → finite-sites (§4). Judging multiallelic FST against the flat 0.237 would
   make a correct estimator look biased. `dxy_fst_theta_summary.tsv` and the inline stats
   do exactly this: each estimator is tested against its own target.

   Figure 4 deliberately does not (decided 2026-07-15). Its FST row plots the single
   finite-sites curve, because the π and dxy rows already define "truth" as the
   finite-sites observed quantity — the π reference is E[π] = 0.088 at θ = 0.1, not θ
   itself — and `E[FST]_finite = 1 - E[pi_w]/E[dxy]` is built from those same two
   expectations. One line per panel, consistent across rows.

   Note the symmetry of the trap, and state the consequence in the caption: judging
   multiallelic FST against the flat 0.237 makes a correct estimator look biased, and
   judging **biallelic** FST against the finite-sites curve does the same thing in
   reverse — it reads as ~1.6% high at θ = 0.1 because it is recovering the *coalescent*
   FST, not because it is broken. No single reference line flatters both. The choice is
   which estimand the figure is about; Figure 4's answer is the one its other two rows
   already assume.
2. **The multiallelic FST is only mildly θ-dependent** (−1.6% at θ = 0.1). This is a
   *better* story than the handoff assumed: the two estimands do not diverge dramatically,
   and the headline asymmetry stands — the same site exclusion that costs π and dxy
   double-digit percentages costs FST ~1%, because FST is a ratio and the saturation
   largely cancels.
3. **The cancellation is quantified here, but do not claim nobody had derived it.** An
   earlier draft of this document said exactly that, and it is wrong — see the prior-art
   note in §2. That F_ST depends on the mutation rate, and that the dependence runs through
   within-population diversity, is established (Rousset 1996; Charlesworth 1998;
   Wilkinson-Herbots 1998), and de Jong et al. (2024) states the coalescent-versus-observed
   distinction for this exact estimator. What §3 adds is the closed form for *our* model, so
   the sweep has a reference line. Claim that, not priority.
4. **Scope honestly.** This is a JC69 result. Under ts/tv bias or rate heterogeneity the
   master equation still holds but `L(beta)` changes. Tajima (1996) also gives the
   gamma-rate-variation case, and reports that with small shape parameter the effect on
   E[π] is *substantial* — which is a concrete reason to expect the open non-JC69 question
   to matter, and a reference for it.

---

## 7. R implementation

To drop into `11_dxy_fst_theta_sweep.R`, replacing the single `theoretical_fst` column:

```r
# Infinite-sites (coalescent) FST — Slatkin (1991); the correct target for the
# BIALLELIC estimator, which conditions away multi-hit sites and so recovers the
# ratio of mean coalescence times.
fst_tau   = split_time / (2 * Ne),
fst_alpha = 2,
theoretical_fst_infinite = (fst_tau + (fst_alpha - 1) * (1 - exp(-fst_tau))) /
                           (fst_tau + fst_alpha),

# Finite-sites FST under Jukes & Cantor (1969) — the correct target for the
# MULTIALLELIC estimator, which measures observed dissimilarity and therefore
# saturates. Both terms are (3/4)(1 - L(beta)) where L is the Laplace transform of
# the pairwise coalescence time; see analysis/finite_sites_fst_derivation.md.
# Sanity: as theta -> 0 this converges to theoretical_fst_infinite.
fs_k    = 1 + (4/3) * theta_nom,
L_dxy   = exp(-(4/3) * theta_nom * fst_tau) / (1 + fst_alpha * (4/3) * theta_nom),
L_pi_w  = (1 - exp(-fst_tau * fs_k)) / fs_k +
            exp(-fst_tau * fs_k) / (1 + fst_alpha * (4/3) * theta_nom),
theoretical_fst_finite = 1 - (1 - L_pi_w) / (1 - L_dxy),
```

Estimator → expectation mapping:

| estimator column | expectation |
|---|---|
| `new` (biallelic) | `theoretical_fst_infinite` |
| `new_multi` (multiallelic) | `theoretical_fst_finite` |
| `new_multi_fstbi` (old mislabelled column) | `theoretical_fst_infinite` |

Note `theoretical_dxy` needs no change — §3.1 confirms it is already correct.

---

## 8. References

Bibliographic metadata below was verified against PubMed records; DOIs link to the
originals. Rousset (1996), Wilkinson-Herbots (1998), Charlesworth (1998) and de Jong et al.
(2024) were added 2026-07-15 after a prior-art search — see the note in §2. Their absence
from the first version of this list is what let it overstate the derivation's novelty.

- **Jukes, T.H. & Cantor, C.R. (1969).** Evolution of protein molecules. In: Munro, H.N.
  (ed.) *Mammalian Protein Metabolism*, Vol. III, pp. 21–132. Academic Press, New York.
  — The substitution model; source of `P(differ | t) = (3/4)(1 - exp(-(4/3) mu t))`.
  *(Book chapter, no DOI; citation details are the conventional ones and were not
  machine-verified — worth a check against the physical volume before submission.)*

- **Tajima, F. (1983).** Evolutionary relationship of DNA sequences in finite populations.
  *Genetics* 105:437–460. [10.1093/genetics/105.2.437](https://doi.org/10.1093/genetics/105.2.437)
  — Equation 15: the original `E(pi) = theta/(1 + 4 theta/3)`.

- **Tajima, F. (1996).** The amount of DNA polymorphism maintained in a finite population
  when the neutral mutation rate varies among sites. *Genetics* 143:1457–1465.
  [10.1093/genetics/143.3.1457](https://doi.org/10.1093/genetics/143.3.1457)
  — E(π) under JC69 finite sites (p. 1458, from Eq. 6 with n = 2) and Table 1; also the
  gamma-rate-variation extension relevant to the open non-JC69 question.

- **Slatkin, M. (1991).** Inbreeding coefficients and coalescence times. *Genetical
  Research* 58:167–175. [10.1017/s0016672300029827](https://doi.org/10.1017/s0016672300029827)
  — FST as the ratio of average coalescence times; the basis of the infinite-sites
  expectation already in the script.

- **Rousset, F. (1996).** Equilibrium values of measures of population subdivision for
  stepwise mutation processes. *Genetics* 142:1357–1362.
  [10.1093/genetics/142.4.1357](https://doi.org/10.1093/genetics/142.4.1357)
  — Relates identity *in state* to the distribution of coalescence times and computes F_ST
  from it; explains F_ST at high mutation rate. **Prior art for §2's approach.**

- **Wilkinson-Herbots, H.M. (1998).** Genealogy and subpopulation differentiation under
  various models of population structure. *Journal of Mathematical Biology* 37:535–585.
  [10.1007/s002850050140](https://doi.org/10.1007/s002850050140)
  — Computes the Laplace transform of the pairwise coalescence time and from it F_ST, its
  θ → 0 limit, and its derivative with respect to the mutation rate; finds F_ST "can depend
  strongly on the mutation rate". **This is §2 and §5.3, for island/stepping-stone models
  under infinite alleles.** Not indexed in PubMed; bibliographic details from the publisher
  page, verify page range before submission.

- **Charlesworth, B. (1998).** Measures of divergence between populations and the effect of
  forces that reduce variability. *Molecular Biology and Evolution* 15:538–543.
  [10.1093/oxfordjournals.molbev.a025953](https://doi.org/10.1093/oxfordjournals.molbev.a025953)
  — Compares F_ST definitions for DNA sequence data; F_ST "is strongly influenced by the
  level of within-population diversity". The qualitative core of this document's result.

- **de Jong, M.J. et al. (2024).** Calculating and interpreting FST in the genomics era.
  *bioRxiv* 2024.09.24.614506.
  [10.1101/2024.09.24.614506](https://doi.org/10.1101/2024.09.24.614506)
  — Advocates Hudson's `F_ST = (dxy − pi_xy)/dxy` on observed sequence dissimilarity, and
  notes F_ST converts to coalescent units only for recent splits "when novel mutations are
  negligible", with no universal F_ST-to-split-time relationship otherwise. **The closest
  prior art to this document: same estimator, same split model, our result stated
  qualitatively.** Preprint — check for a published version before submission.

- **Hudson, R.R., Slatkin, M. & Maddison, W.P. (1992).** Estimation of levels of gene flow
  from DNA sequence data. *Genetics* 132:583–589.
  [10.1093/genetics/132.2.583](https://doi.org/10.1093/genetics/132.2.583)
  — The Hudson FST estimator that pixy implements.

- **Bhatia, G., Patterson, N., Sankararaman, S. & Price, A.L. (2013).** Estimating and
  interpreting FST: the impact of rare variants. *Genome Research* 23:1514–1521.
  [10.1101/gr.154831.113](https://doi.org/10.1101/gr.154831.113)
  — Attributes the Hudson estimator to Hudson et al. (1992) and recommends the
  ratio-of-averages combination across sites; "Therefore, we recommend using a ratio of
  averages." This is what pixy computes and what the ΣN/ΣD pooling in the analysis
  reproduces.

- **Nei, M. & Li, W.-H. (1979).** Mathematical model for studying genetic variation in
  terms of restriction endonucleases. *PNAS* 76:5269–5273.
  [10.1073/pnas.76.10.5269](https://doi.org/10.1073/pnas.76.10.5269)
  — Defines nucleotide diversity and the net-divergence correction for ancestral
  polymorphism.

- **Takahata, N. & Nei, M. (1985).** Gene genealogy and variance of interpopulational
  nucleotide differences. *Genetics* 110:325–344.
  [10.1093/genetics/110.2.325](https://doi.org/10.1093/genetics/110.2.325)
  — Coalescent treatment of between-population differences with ancestral polymorphism;
  the classical basis for the `T + T_a` decomposition in §3.1.

- **Cutter, A.D., Jovelin, R. & Dey, A. (2013).** Molecular hyperdiversity and evolution in
  very large populations. *Molecular Ecology* 22:2074–2095.
  [10.1111/mec.12281](https://doi.org/10.1111/mec.12281)
  — Why this regime matters empirically: reports ~10% of segregating sites in
  *C. brenneri* carrying three or four nucleotide variants, and argues for finite-sites
  models over a post-hoc Jukes–Cantor multiple-hits correction. Directly supports the
  paper's case for multiallelic-aware estimators.
