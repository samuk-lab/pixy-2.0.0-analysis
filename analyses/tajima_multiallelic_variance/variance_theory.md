# Finite-sites Tajima's D: mutation counting and the variance term

_Rewritten 2026-07-15. The previous version of this file (git history) argued that the
mutation-count substitution — which it called "Route A" — **was not a fix**, and
recommended a simulation-calibrated variance ("Route B") instead. **That conclusion was
wrong, and it was wrong because of a bug.** It is retracted here; see §6 for the full
retraction and why the error was not detectable from the document itself. The shipped
pixy branch `multiallelic-mutation-count-theta-d` implements the mutation-count
estimator, i.e. exactly what the old §4 argued against._

## 1. Tajima's D and its variance

Tajima's *D* contrasts two estimators of θ = 4Nₑμ — nucleotide diversity π (mean
pairwise differences) and Watterson's θ_W = S / a₁, with a₁ = Σ_{i=1}^{n−1} 1/i — and
standardises the difference (Tajima 1989):

    D = (π − θ_W) / sqrt( V̂(π − θ_W) ),
    V̂(π − θ_W) = e₁·S + e₂·S(S − 1),
    e₁ = c₁ / a₁,  e₂ = c₂ / (a₁² + a₂),  a₂ = Σ 1/i².

## 2. What the estimator actually counts

Both terms of V̂ are method-of-moments plug-ins. Under the infinite-sites coalescent
(Kimura 1969; Watterson 1975) the number of mutations on the genealogy has moments

    E[S] = a₁ θ,     E[S(S − 1)] = (a₁² + a₂) θ²

(the second from the coalescent tree-length moments E[L] = 2a₁, Var(L) = 4a₂; Wakeley
2009). So S/a₁ and S(S−1)/(a₁²+a₂) are unbiased for θ and θ², and substituting them into
Var(π − θ_W) = e₁'θ + e₂'θ² gives V̂. These are statements about the number of
mutations on the genealogy, not the number of variable positions. Under infinite sites
each mutation makes a new site, so "segregating sites" = "mutations" and the distinction
is invisible.

## 3. Recurrent mutation splits sites from mutations

pixy's polyploid validation simulates a finite-sites Jukes–Cantor process (Jukes &
Cantor 1969): a site can be hit repeatedly and carry 3–4 alleles. Then

    S_sites = number of segregating sites,
    M       = Σ_sites (k − 1)   (k = alleles at a site).

M is exactly **Tajima's `s*`, the minimum (parsimony) number of mutations per site**
(Tajima 1996, Eq 14), and it — not S_sites — is the quantity whose expectation tracks the
genealogical mutation count. Tajima 1996 Eq 15 gives its closed-form expectation under
JC69. So the analytic fix is to substitute M for S everywhere: θ_W → M/a₁ and V̂ → e₁M +
e₂M(M−1). Because every biallelic site has k−1 = 1, this is byte-identical to Tajima 1989
on biallelic data and departs only at multiallelic sites.

A multiallelic site is *by definition* a departure from infinite sites. The `k−1` rule is
what you use once you have left that model, not a consequence of it — the correct
citation is Tajima 1996's `s*`, not "the infinite-sites assumption."

## 4. The mutation-count substitution works — the evidence

**θ_W hits a published closed form.** `tajima1996_expectations.py` codes Tajima 1996's
E(s\*) (Eq 15) and reproduces the paper's Table 1. The pipeline's η θ_W matches
E(s\*)/a₁ to ~1e-5 across ploidy (se ≈ 6e-4, i.e. ~0.02 SE), at 100–120k replicates:

| ploidy | η θ_W     | E(s\*)/a₁ | diff   |
|--------|-----------|-----------|--------|
| 2n     | 0.0362678 | 0.0362599 | +8e-6  |
| 4n     | 0.0362088 | 0.0362078 | +1e-6  |
| 6n     | 0.0361522 | 0.0361649 | −1.3e-5|
| 8n     | 0.0361260 | 0.0361302 | −4e-6  |

This is an external anchor: an estimator landing on a closed form from the literature,
not a simulation matching itself. (`analysis/thetaw_tajimaD_summary.tsv`.)

**π is exact.** pixy's multiallelic π equals tskit's unbiased pairwise π to ~1e-16 on
shared msprime JC69 matrices. Both numerator arms are independently validated.

**The variance term is fine — measured, not assumed.** See §5.

## 5. The variance term, measured at true polyploidy (`variance_validate.py`, 2026-07-15)

The `e₂·M(M−1)` term treats two mutations sharing one site as two independent units. That
is a real approximation and it was, until now, the one component never externally checked.
`variance_validate.py` checks it by driving msprime straight to a `GenotypeArray` (no VCF,
so the reader bug of §6 cannot apply), at r = 0 — the only regime where Tajima's variance
is derived and sd(D) = 1 is the correct target — with an **infinite-sites arm as a positive
control**, where η ≡ S and the substitution is a literal no-op.

2000 reps/cell, 10 kb, Ne/μ as in the Figure 2/3 arms, se(sd) ≈ 0.014 (r=0) and
≈ 0.0025 (r=μ). Output: `results/variance_sweep.tsv`.

| ploidy | r | model    | multi. frac | sd(D_old) | sd(D_new) | mean(D_old) | mean(D_new) |
|--------|---|----------|-------------|-----------|-----------|-------------|-------------|
| 2n     | 0 | infinite | 0.000       | 0.8768    | 0.8768    | −0.0979     | −0.0979     |
| 4n     | 0 | infinite | 0.000       | 0.8703    | 0.8703    | −0.1137     | −0.1137     |
| 6n     | 0 | infinite | 0.000       | 0.8771    | 0.8771    | −0.1325     | −0.1325     |
| 8n     | 0 | infinite | 0.000       | 0.8898    | 0.8898    | −0.1010     | −0.1010     |
| 2n     | 0 | jc69     | 0.0352      | 0.8656    | 0.8290    | +0.0205     | −0.1195     |
| 4n     | 0 | jc69     | 0.0459      | 0.8902    | 0.8389    | +0.0855     | −0.0834     |
| 6n     | 0 | jc69     | 0.0508      | 0.9022    | 0.8463    | +0.0633     | −0.1146     |
| 8n     | 0 | jc69     | 0.0546      | 0.8995    | 0.8380    | +0.0419     | −0.1434     |
| 2n     | μ | infinite | 0.000       | 0.1560    | 0.1560    | −0.0007     | −0.0007     |
| 4n     | μ | infinite | 0.000       | 0.1552    | 0.1552    | +0.0041     | +0.0041     |
| 6n     | μ | infinite | 0.000       | 0.1559    | 0.1559    | −0.0048     | −0.0048     |
| 8n     | μ | infinite | 0.000       | 0.1524    | 0.1524    | −0.0004     | −0.0004     |
| 2n     | μ | jc69     | 0.0380      | 0.1552    | 0.1499    | **+0.1159** | **−0.0349** |
| 4n     | μ | jc69     | 0.0481      | 0.1530    | 0.1451    | **+0.1461** | **−0.0295** |
| 6n     | μ | jc69     | 0.0532      | 0.1572    | 0.1484    | **+0.1683** | **−0.0193** |
| 8n     | μ | jc69     | 0.0573      | 0.1572    | 0.1484    | **+0.1805** | **−0.0168** |

The bolded r=μ rows are the Figure 2/3 regime and carry the headline result on their own:
site counting rises **+0.116 → +0.181** with ploidy against an infinite-sites control that
sits at **0.000**, while mutation counting sits at **−0.035 → −0.017**, flat. Those
mutation-count values independently reproduce the pipeline's own `new_multi` D
(−0.033/−0.025/−0.019/−0.016) and the `calibrate_neutral_tajimaD.py` floor
(−0.033/−0.025/−0.020/−0.015) to ~0.002, from a harness that shares no code with either.

Three findings, in order of importance.

**(a) Tajima's D is under-dispersed under INFINITE sites, in BOTH regimes.** The control
gives sd ≈ 0.88 at r = 0 and ≈ 0.155 at r = μ, with no homoplasy and the patch inactive. D
is a ratio of a mean-zero numerator to a random denominator, so neither E[D] = 0 nor
Var(D) = 1 holds exactly (Tajima 1989 uses a beta approximation for exactly this reason;
the control's mean(D) ≈ −0.11 at r = 0 is the same effect — note it washes out to ≈ 0.000
at r = μ). **This retires the old §4(c) claim** that the measured "sd ≈ 0.85 at r = 0" was
evidence of a problem with the substitution: it is baseline Tajima behaviour, present where
the substitution does nothing at all. It also corrects the old numbers themselves — 0.85 /
0.19 came from the truncated run *and* at θ = 0.025 rather than the figures' 0.0378.

**(b) The η variance costs 2–6% of sd, and does not grow with ploidy.** Against the
infinite-sites control at r = 0, sd(D_new) is −5.5% / −3.6% / −3.5% / −5.8% at 2n/4n/6n/8n.
Site counting drifts a comparable amount the *other* way (−1.3% / +2.3% / +2.9% / +1.1%).
Multiallelic fraction climbs 3.5% → 5.5% across ploidy while the gap does not trend, so
the independence approximation in `e₂·M(M−1)` is not the driver. For scale: recombination
shrinks sd(D) from ≈ 0.88 to ≈ 0.155 — a ~5.7× effect that applies to stock pixy and every
other Tajima's D implementation equally. The η approximation is a rounding error next to
the error already in the formula.

**(c) The numerator correction is confirmed against an external ideal.** In the Figure 2/3
regime (r = μ), the no-homoplasy control sits at **0.000**; site-counting D is inflated to
**+0.116 / +0.146 / +0.168 / +0.181** (2n→8n), growing with ploidy, while mutation-counting
D sits at **−0.035 / −0.030 / −0.019 / −0.017**, flat. Unlike the D floor in
`neutral_tajimaD_reference.tsv` (which is produced by `eta_stats`, the same machinery under
test, and so can only check the *implementation*), the infinite-sites arm is a genuine
external reference: it is what D does when there is no homoplasy to correct for.

**Verdict: the variance term is not a barrier to shipping.** The README's "open question"
is closed. What remains true is that *absolute* D values are not calibrated test
statistics — but that was true before this change and is true of every implementation; see
§7.

## 6. Retraction: why the old §4 said the opposite

The previous version of this file concluded the substitution "is not a fix," on three
grounds. All three fail:

**Old §4(a) — "a monotone downshift, not a centering correction."** This was the load-bearing
empirical claim, and it is an **artifact of a bug**. Its evidence was `results/run.log`
(2026-06-26 21:06), produced by `analyze.py`, which called
`allel.read_vcf(path, fields=["calldata/GT"])` **without `numbers=`**. That defaults to
diploid and silently truncates every polyploid genotype to its first 2 alleles. Every
"polyploid" cell in that run was diploid. With no polyploidy there is almost no
multiallelism, so η ≈ S and the substitution *can only* nudge D slightly negative — which
the document then read as "a downshift with nothing to fix." The tell is visible in the
log itself: `D_old` at r=μ was −0.016 / −0.003 / −0.008 / −0.003 across 2n–8n — **flat, no
ploidy trend** — where the correct re-sim (`eta_validate.py`, fixed reader) gives +0.069 /
+0.132 / +0.157 / +0.180. *The signal the substitution exists to correct was absent from
the data used to reject it.* This file was written at 21:39, 33 minutes after that log;
the reader bug was found the next morning.

**Old §4(b) — "Watterson under finite sites estimates a composite parameter"
(Roychoudhury & Wakeley 2010).** True, and unchanged. But it is an argument about the
*magnitude* of θ_W, not against mutation counting: it predicts exactly the ~4% shortfall
below 4Nₑμ that we measure and that Tajima 1996's E(s\*) predicts analytically. A known,
quantified, flat-in-ploidy floor with a closed form is not a defect. It is also not
avoidable — no exact sampling formula exists for general finite-allele models (Bhaskar,
Kamm & Song 2012).

**Old §4(c) — "recombination invalidates the analytic variance regardless."** True, and
unchanged — but it indicts stock pixy and every other Tajima's D equally, so it is not an
argument against this change. The specific number it leaned on (sd ≈ 0.85 at r=0) is now
shown by §5(a) to be baseline behaviour, not a symptom.

**The methodological lesson.** The bug was invisible from inside the analysis: every cell
ran, every number was plausible, and the conclusion was self-consistent. What exposed it
was a *cross-check against a different measurement of the same quantity* — the re-sim's
ploidy trend in `D_old`, which the truncated run could not reproduce. Two guards now
exist: `variance_validate.py` never writes a VCF (the reader cannot enter the loop), and
it carries an infinite-sites arm whose invariant (`D_old == D_new` exactly, since η ≡ S)
fails loudly if the harness or the estimator is wrong.

## 7. What is still NOT claimed

The mutation-count estimator is an **estimator** fix, not a **test** fix. Under finite
sites and recombination, an absolute D value is not a calibrated statistic:

- E[D] under neutrality is not 0 — it is the finite-sites floor (~−0.01 to −0.04 here,
  deepening with missingness), plus the standardisation's own ~−0.11 offset at r = 0.
- sd(D) is not 1 — it is ~0.88 at r = 0 and ~0.155 at r = μ.

So D is sound for **comparison across windows or groups computed the same way**, which is
what pixy is for and what Figures 2–3 do. It is not a p-value. A calibrated finite-sites
neutrality test needs a simulated null matched on θ̂, ploidy, missingness and *per-window
recombination* (the old "Route B"), and the parked work in this folder found no ρ-free
shortcut: block-studentised self-norm is anti-conservative and blows up at r=0; a
genome-wide MAD z-score over-calls cold spots (FPR to 0.65) and goes blind in hot ones.
That remains a separate research problem.

## References

- Achaz, G. (2008). Testing for neutrality in samples with sequencing errors.
  *Genetics* 179: 1409–1424.
- Bhaskar, A., Kamm, J.A. & Song, Y.S. (2012). Approximate sampling formulae for
  general finite-alleles models of mutation. *Adv. Appl. Probab.* 44: 408–428
  (arXiv:1109.2386).
- Jukes, T.H. & Cantor, C.R. (1969). Evolution of protein molecules. In
  *Mammalian Protein Metabolism* (ed. Munro), pp. 21–132. Academic Press.
- Kimura, M. (1969). The number of heterozygous nucleotide sites maintained in a
  finite population due to steady flux of mutations. *Genetics* 61: 893–903.
- Misawa, K. & Tajima, F. (1997). Finite-site estimators of θ and modified
  neutrality statistics. *Genetics* 147: 1959–1964.
- Roychoudhury, A. & Wakeley, J. (2010). Sufficiency of the number of segregating
  sites in the limit under finite-sites mutation. *Theoretical Population
  Biology* 78: 118–122.
- Tajima, F. (1989). Statistical method for testing the neutral mutation
  hypothesis by DNA polymorphism. *Genetics* 123: 585–595.
- Tajima, F. (1996). The amount of DNA polymorphism maintained in a finite
  population when the neutral mutation rate varies among sites. *Genetics*
  143: 1457–1465.
- Watterson, G.A. (1975). On the number of segregating sites in genetical models
  without recombination. *Theoretical Population Biology* 7: 256–276.
- Wakeley, J. (2009). *Coalescent Theory: An Introduction.* Roberts & Company.
- Yang, Z. (1996). Finite-sites models and the distribution of the number of
  segregating sites under multiple-hit substitution. *J. Mol. Evol.* 43.
