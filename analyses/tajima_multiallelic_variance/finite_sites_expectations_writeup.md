# Finite-sites expectations for Watterson's θ and Tajima's *D* under multiallelic, polyploid sampling

## Summary

In multiallelic-aware mode, pixy estimates nucleotide diversity (π) from allelic
pairwise differences but estimates Watterson's θ from the number of segregating
*sites*. Under an infinite-sites mutation model these two currencies are the same
integer, so the choice does not matter. Under a finite-sites model with recurrent
mutation they diverge, and the imbalance biases Tajima's *D* upward. The bias
grows as more gene copies are sampled (ploidy × sample size). Here, we re-derive
both estimators in terms of the number of *mutations* rather than the number of
*sites*, and we give the finite-sites expectations of θ̂ and *D* — the reference
lines the simulated polyploid data are compared against — using the closed-form
Jukes–Cantor results of Tajima (1996). This note records the derivation and the
reasoning behind each step.

## 1. The two θ estimators are equivalent only under infinite sites

Tajima's (1989) test contrasts two method-of-moments estimators of the scaled
mutation parameter θ = 4*N*<sub>e</sub>μ. The first is nucleotide diversity,

$$\hat{\pi} = \text{mean number of pairwise differences},$$

and the second is Watterson's (1975) estimator,

$$\hat{\theta}_W = \frac{S}{a_1}, \qquad a_1 = \sum_{i=1}^{n-1}\frac{1}{i},$$

where *S* is the number of segregating sites in a sample of *n* gene copies. Under
neutrality both are unbiased for θ, so their standardised difference,

$$D = \frac{\hat{\pi} - \hat{\theta}_W}{\sqrt{\widehat{\mathrm{Var}}(\hat{\pi}-\hat{\theta}_W)}},
\qquad
\widehat{\mathrm{Var}}(\hat{\pi}-\hat{\theta}_W) = e_1 S + e_2 S(S-1),$$

has expectation zero. The variance constants are $e_1 = c_1/a_1$ and
$e_2 = c_2/(a_1^2 + a_2)$ with $a_2 = \sum_{i=1}^{n-1} 1/i^2$ (Tajima 1989).

Both terms of the variance are moment plug-ins, and they rest on the coalescent
moments of the number of mutations carried by the genealogy (Kimura 1969; Watterson
1975):

$$E[S] = a_1\theta, \qquad E[S(S-1)] = (a_1^2 + a_2)\theta^2,$$

the second following from the tree-length moments $E[L] = 2a_1$, $\mathrm{Var}(L) = 4a_2$
(Wakeley 2009). These are statements about the number of mutations placed on the
genealogy, not about the number of variable positions in the alignment. Under
infinite sites every mutation lands on a fresh position, so "segregating sites" and
"mutations" are the same count and the distinction never surfaces. This is the
assumption the finite-sites case breaks.

## 2. Recurrent mutation separates sites from mutations

pixy's polyploid validation simulates a finite-sites Jukes–Cantor process (Jukes &
Cantor 1969) in which a single position can be struck more than once and can
therefore segregate for three or four alleles. Once that happens, two counts
diverge:

$$S = \#\{\text{segregating sites}\}, \qquad
\eta = \sum_{\text{sites}} (k_i - 1),$$

where $k_i$ is the number of alleles observed at site *i*. A site with $k$ alleles
records a minimum of $k-1$ mutations under the no-homoplasy reading, so η is the
parsimony-minimum mutation count over the sample. This is exactly Tajima's (1996)
$s^\*$, defined there as "the number of nucleotides minus one" summed over sites
($s^\* = q_2 + 2q_3 + 3q_4$, his Eq. 14).

The consequence for pixy is a mismatch in counting currency. Multiallelic-aware π is
built from allelic pairwise differences, so a three- or four-allele site contributes
its full heterozygosity — that is, π already counts mutations. Watterson's θ as
originally implemented counts sites, tallying one unit per variable position
regardless of allele number. Recurrent mutation therefore raises π without raising
$\hat{\theta}_W$ by the same amount, so $\hat{\pi} > \hat{\theta}_W$ and *D* is biased
positive. A site is more likely to pick up a third or fourth allele as more lineages
are sampled, so the bias grows with ploidy. This is the source of the positive,
ploidy-increasing Tajima's *D* seen across the polyploid arms of the manuscript.

## 3. The correction: count mutations, not sites

The fix is to use η wherever the segregating-site count *S* appeared, so that both
arms of *D* count the same thing. Watterson's estimator becomes a per-site-pooled
mutation count,

$$\hat{\theta}_W^{\eta} = \sum_i \frac{k_i - 1}{a_1(n_i)},$$

evaluated at the number of gene copies $n_i$ actually observed at each site — which
matters under missing data and mixed ploidy, where $n_i$ varies among sites — and
the Tajima variance takes the same substitution,

$$D_\eta = \frac{\hat{\pi} - \hat{\theta}_W^{\eta}}
{\sqrt{\; \sum_n \big[\, e_1(n)\,\eta_n + e_2(n)\,\eta_n(\eta_n - 1) \,\big]}},$$

where $\eta_n$ is the total mutation count among sites with $n$ observed copies. The
construction is self-consistent by design: every biallelic site has $k-1 = 1$, so
$\eta = S$ and both estimators reduce to the classical Tajima (1989) forms on
biallelic data; only multiallelic sites change. The implementation reuses pixy's own
harmonic-sum and variance routines with the per-site-class mutation count in place
of the site count (`eta_tajima.py`).

We note two limitations that bound what this correction does. First, it is a
consistent estimator, not an unbiased finite-sites test statistic. Because
$\eta \ge S$ always, $\hat{\theta}_W^{\eta} \ge \hat{\theta}_W$, and the contrast
$\hat{\pi} - \hat{\theta}_W^{\eta}$ only moves downward relative to the stock
statistic. On simulated JC69 data, which contains homoplasy, $D_\eta$ does not
re-centre exactly at zero. That residual is the finite-sites floor rather than a
defect of the estimator. Measured against a no-homoplasy infinite-sites control, in
which every site is biallelic so that $\eta \equiv S$ and the correction is inactive
by construction, the site-counting statistic is inflated by +0.12 to +0.18 and rises
with ploidy, whereas the mutation-count form sits within 0.04 of the control at every
ploidy and is flat in ploidy. We note that $s^\*$ is a finite-sites, parsimony
construct (Tajima 1996): under strict infinite sites no site carries more than two
alleles, $\eta \equiv S$, and the two estimators coincide exactly. It is where that
model fails that the mutation count and the site count part company.
Second, under finite sites Watterson estimates a model-dependent composite rather
than θ itself. Roychoudhury & Wakeley (2010) show that under finite-sites mutation
*S* (and likewise $s^\*$) remains a moment estimator, but of a
mutation-model-dependent quantity rather than a universal θ. A genuinely unbiased
finite-sites *D* would have to commit to a mutation model and derive $E[\hat{\pi}]$,
$E[s^\*]$, and their covariance under it (Tajima 1996; Misawa & Tajima 1997), and
even then the standardised ratio is not exactly unbiased. No general, model-free
closed form exists, because exact sampling formulae are unknown for general
finite-allele mutation models (Yang 1996; Bhaskar, Kamm & Song 2012). This is why
the corrected estimator, though consistent, leaves a small irreducible finite-sites
offset, which we treat below as the reference line rather than as a bias to remove.

## 4. Finite-sites expectations: the closed forms

The reference lines the polyploid figures are compared against are the expectations
of these estimators under the same finite-sites JC69 model used to simulate the
data. All three follow from Tajima (1996), and we reproduced each one exactly
against his Table 1 (`tajima1996_expectations.py`).

For a sample of *n* gene copies, define the probability that a site is occupied by
only *m* specified nucleotides, written with rising factorials
$x^{(n)} = \Gamma(x+n)/\Gamma(x)$:

$$P_m(\theta, n) = \frac{(m\theta/3)^{(n)}}{(4\theta/3)^{(n)}}.$$

The three quantities of interest are then

$$E[\pi] = \frac{\theta}{1 + \tfrac{4}{3}\theta} \qquad\text{(Tajima 1996, Eq. 8)},$$

$$E[s] = 1 - 4\,P_1(\theta,n) \qquad\text{(Eq. 6)},$$

$$E[s^\*] = 3 - 4\,P_3(\theta,n) \qquad\text{(Eq. 15)},$$

where $E[s]$ is the expected proportion of segregating sites and $E[s^\*]$ the
expected minimum number of mutations per site. The expectations of the two Watterson
estimators follow by dividing by the harmonic number:

$$E\!\left[\hat{\theta}_W\right] = \frac{E[s]}{a_1(n)}
\quad\text{(site count)}, \qquad
E\!\left[\hat{\theta}_W^{\eta}\right] = \frac{E[s^\*]}{a_1(n)}
\quad\text{(mutation count).}$$

Three properties of these expressions carry the interpretation. First, $E[\pi]$
saturates but does not depend on sample size: the $\theta/(1+\tfrac43\theta)$ form is
the JC69 correction for multiple hits on a single pairwise lineage, a function of θ
alone. It gives the closed-form inversion we use to set the simulation-calibrated
null, $\hat{\theta}_\pi = \hat{\pi}/(1 - \tfrac43\hat{\pi})$. Second, both θ
expectations fall below θ, and the mutation-count one falls less. Because $s^\*$
counts the minimum mutations implied by the observed alleles, it still misses
recurrent hits that leave no allelic trace, so $E[s^\*]/a_1$ sits a few percent below
the true $4N_e\mu$. In the polyploid sims this is the ~4% deficit of
$\hat{\theta}_W^{\eta}$ relative to the target — an irreducible homoplasy floor, not a
coding error, consistent with the absence of a general finite-allele closed form
(Yang 1996; Bhaskar, Kamm & Song 2012). The site-count expectation $E[s]/a_1$ lies
further below still and, unlike the mutation-count form, declines with *n* as sites
saturate to a fixed allele count. Third, the expectations must be evaluated at the
effective sample size. Both Watterson forms carry $a_1(n)$ and both per-site
expectations depend on *n*, so under missing data they must be evaluated at the
reduced per-site copy number $n_\text{eff}$ actually observed, not at the nominal
ploidy × sample size. vcfsim missingness drops whole individuals, so $n_\text{eff}$
falls in steps of the ploidy and the θ reference slopes upward with missingness;
evaluating the closed form at $n_\text{eff}$ reproduces that slope analytically, with
no re-simulation.

## 5. The Tajima's *D* reference line

Tajima's *D* has no convenient finite-sites closed form, for the reasons in §3: it
is a ratio of a saturating numerator and a mutation-model-dependent denominator whose
joint distribution has no general analytic sampling formula. Three further
complications rule out the Tajima (1989) analytic variance as the reference even in
principle. First, that variance is derived for a single non-recombining locus,
whereas the polyploid arms simulate recombination at rate $r = \mu$; averaging over
many independent genealogies within a window shrinks the standard deviation of *D*
several-fold (measured sd ≈ 0.88 at $r=0$ vs ≈ 0.16 at $r=\mu$). Second, the
departure from unit variance is not specific to finite sites. Under an
infinite-sites control, in which recurrent mutation is absent and the correction is
inactive, we measure the same dispersion (sd ≈ 0.88 at $r=0$ and ≈ 0.16 at
$r=\mu$) together with a mean of ≈ −0.11 at $r=0$: *D* is a ratio whose denominator
is itself a random variable, and Tajima (1989) introduced a beta approximation for
precisely this reason. The departure is therefore a property of the statistic, not a
finite-sites artifact, and it is present whether or not the correction is applied.
Third, the finite-sites numerator and denominator are each biased in
mutation-model-specific ways.

We therefore obtain the *D* reference by simulation calibration rather than by
formula. A neutral coalescent null is simulated in msprime under a JC69
discrete-genome model, so that recurrent mutation and multiallelic sites arise as
they do in the data, matched to the window length, recombination rate, ploidy,
sample size, and missingness of each cell, at θ̂ estimated from the data through the
JC69 π-inversion above. The reference is the mean of the corrected η-statistic over
that null. The property that makes this work is that both the bias and the variance
of the contrast under recurrent mutation are absorbed into the simulated null, so
neither arm has to be individually unbiased and the analytic-variance and
composite-parameter problems of §3 do not arise. The cost is a commitment to a
mutation model and to θ̂, which §3 argues is intrinsic to any finite-sites *D*.

The result is that the correct expected *D* under finite sites is not zero but a
small negative floor (−0.033 at 2n, shallowing to −0.016 at 8n). This floor
reflects the homoplasy asymmetry between the saturating π numerator and the
minimum-mutation denominator. The corrected pipeline estimates sit on this floor. This agreement is a
consistency check on the implementation rather than an independent calibration,
because the floor is itself the mean of the corrected statistic over the null; the
estimator's external validation is its exact agreement with the Tajima (1996) closed
form reported in §6. What the floor does establish is that the large, ploidy-growing
positive *D* of the stock multiallelic statistic (≈ +0.07 to +0.18) was the
site-vs-mutation counting artifact, and once it is removed the residual is the
expected finite-sites floor rather than a residual bias. We note two subtleties in
setting the reference. The coalescent timescale in the calibration must be set so
that the simulated θ equals $4N_e\mu$ — msprime's top-level `ploidy` argument scales
the timescale, and omitting it silently halves θ. And, as for θ, the floor must be
evaluated per (ploidy × missingness) cell, because it moves with the effective
sample size.

## 6. Validation

The closed forms of §4 reproduce every entry of Tajima's (1996) Table 1 for $E[\pi]$,
$E[s]/a_1$, and $E[s^\*]/a_1$ to four decimal places across θ ∈ {0.01, …, 0.1} and
*n* ∈ {20, 200}. On the polyploid sims the corrected estimator's simulated mean
matches the analytic $E[s^\*]/a_1(n)$ to ≤ 4 × 10⁻⁵ across ploidy, so the estimator's
theoretical target is an exact closed form and not merely a simulation. The
mutation-count Tajima's *D* is flat in ploidy and sits on the simulated finite-sites
floor, whereas the site-count statistic rises monotonically with ploidy — the
signature of the counting mismatch this correction removes.

## References

Achaz, G. (2008). Testing for neutrality in samples with sequencing errors.
*Genetics* 179: 1409–1424.

Bhaskar, A., Kamm, J. A. & Song, Y. S. (2012). Approximate sampling formulae for
general finite-alleles models of mutation. *Advances in Applied Probability* 44:
408–428. (arXiv:1109.2386.)

Jukes, T. H. & Cantor, C. R. (1969). Evolution of protein molecules. In H. N. Munro
(ed.), *Mammalian Protein Metabolism*, pp. 21–132. Academic Press, New York.

Kimura, M. (1969). The number of heterozygous nucleotide sites maintained in a
finite population due to steady flux of mutations. *Genetics* 61: 893–903.

Misawa, K. & Tajima, F. (1997). Finite-site estimators of the number of nucleotide
substitutions and modified neutrality statistics. *Genetics* 147: 1959–1964.

Roychoudhury, A. & Wakeley, J. (2010). Sufficiency of the number of segregating
sites in the limit under finite-sites mutation. *Theoretical Population Biology* 78:
118–122.

Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by
DNA polymorphism. *Genetics* 123: 585–595.

Tajima, F. (1996). The amount of DNA polymorphism maintained in a finite population
when the neutral mutation rate varies among sites. *Genetics* 143: 1457–1465.
doi:10.1093/genetics/143.3.1457.

Wakeley, J. (2009). *Coalescent Theory: An Introduction.* Roberts & Company,
Greenwood Village, Colorado.

Watterson, G. A. (1975). On the number of segregating sites in genetical models
without recombination. *Theoretical Population Biology* 7: 256–276.

Yang, Z. (1996). Finite-sites models and the distribution of the number of
segregating sites under multiple-hit substitution. *Journal of Molecular Evolution*
43. [volume/pages?]
