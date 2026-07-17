"""
closed-form finite-sites (Jukes-Cantor, constant rate) expectations from
Tajima (1996), Genetics 143:1457-1465. doi:10.1093/genetics/143.3.1457

analytic targets for pixy's watterson's theta estimators:

  * E(s)/a1(n)   -- site-counting watterson estimator (Eq 6). pixy's multiallelic
                    theta_W, the e_thetaW_jc69 line in 07_thetaw_tajimaD.R.
  * E(s*)/a1(n)  -- mutation-counting (eta = sum k-1) estimator (Eq 15). the
                    branch multiallelic-mutation-count-theta-d.
  * E(pi)        -- JC69 saturation of pairwise differences (Eq 8).

Tajima's s* is our eta: "the minimum number of mutations is the number of
nucleotides minus one", s* = q2 + 2 q3 + 3 q4 (Eq 14), E(s*) = 3 - 4 p_ijk
(Eq 15), p_ijk = prob a particular nucleotide is absent from the sample.

subset-occupancy prob (site uses only m specified nucleotides) via rising
factorials x^(n) = Gamma(x+n)/Gamma(x):

    P_m(theta, n) = (m*theta/3)^(n) / (4*theta/3)^(n)

    E(s)  = 1 - 4 * P_1     (monomorphic for one nucleotide)
    E(s*) = 3 - 4 * P_3     (one nucleotide absent)

theta = 4*N*v (N effective size, v mutation rate per site per generation);
n = sampled gene copies. verified against Tajima 1996 Table 1 (run as script). pure stdlib.
"""

from __future__ import annotations

from math import exp, lgamma


def a1(n: int) -> float:
    """Harmonic sum a1(n) = sum_{j=1}^{n-1} 1/j."""
    return sum(1.0 / j for j in range(1, n))


def _p_subset(m: int, theta: float, n: int) -> float:
    """Prob a site is occupied only by ``m`` specified nucleotides: (m*theta/3)^(n)/(4*theta/3)^(n)."""
    return exp(
        lgamma(4 * theta / 3) + lgamma(m * theta / 3 + n)
        - lgamma(4 * theta / 3 + n) - lgamma(m * theta / 3)
    )


def e_pi_jc69(theta: float) -> float:
    """E(pi) under JC69 (Tajima 1996 Eq 8): theta / (1 + 4*theta/3)."""
    return theta / (1.0 + 4.0 * theta / 3.0)


def e_s_jc69(theta: float, n: int) -> float:
    """E(s), expected proportion of segregating sites (Tajima 1996 Eq 6): 1 - 4*P_1."""
    return 1.0 - 4.0 * _p_subset(1, theta, n)


def e_sstar_jc69(theta: float, n: int) -> float:
    """E(s*), expected minimum mutations per site (Tajima 1996 Eq 15): 3 - 4*P_3."""
    return 3.0 - 4.0 * _p_subset(3, theta, n)


def e_thetaW_s_jc69(theta: float, n: int) -> float:
    """Expectation of the SITE-count Watterson estimator: E(s)/a1(n)."""
    return e_s_jc69(theta, n) / a1(n)


def e_thetaW_sstar_jc69(theta: float, n: int) -> float:
    """Expectation of the MUTATION-count (eta) Watterson estimator: E(s*)/a1(n)."""
    return e_sstar_jc69(theta, n) / a1(n)


if __name__ == "__main__":
    # reproduce Tajima 1996 Table 1 (JC69, no rate variation) as a self-test
    table1 = {
        (0.01, 20): (0.0099, 0.0098, 0.0099),
        (0.02, 20): (0.0195, 0.0192, 0.0196),
        (0.05, 20): (0.0469, 0.0451, 0.0474),
        (0.05, 200): (0.0469, 0.0429, 0.0469),
        (0.1, 20): (0.0882, 0.0817, 0.0900),
        (0.1, 200): (0.0882, 0.0744, 0.0883),
    }
    print(f"{'theta':>6}{'n':>5} | {'E(pi)':>7}{'[paper]':>9} "
          f"{'E(s)/a1':>8}{'[paper]':>9} {'E(s*)/a1':>9}{'[paper]':>9}")
    ok = True
    for (th, n), (pi_p, s_p, ss_p) in table1.items():
        pi, s, ss = e_pi_jc69(th), e_thetaW_s_jc69(th, n), e_thetaW_sstar_jc69(th, n)
        ok &= (round(pi, 4) == pi_p and round(s, 4) == s_p and round(ss, 4) == ss_p)
        print(f"{th:>6}{n:>5} | {pi:>7.4f}{pi_p:>9} {s:>8.4f}{s_p:>9} {ss:>9.4f}{ss_p:>9}")
    print("\nTable 1 reproduced exactly." if ok else "\nMISMATCH vs Table 1!")
