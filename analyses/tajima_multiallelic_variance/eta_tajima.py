"""
mutation-count multiallelic watterson's theta and tajima's d.

pixy's multiallelic estimators count segregating sites (one per site, any k), but
multiallelic pi counts allelic differences; the mismatch biases d positive with
ploidy. fix counts mutations instead of sites:

    eta = sum_i (k_i - 1)        k_i = observed alleles at variant site i

used in both the watterson numerator and the tajima-d variance:

    theta_W_eta = sum_i (k_i - 1) / a1(n_i)
    D_eta       = (raw_pi - theta_W_eta) / sqrt( sum_n e1(n) eta_n + e2(n) eta_n (eta_n - 1) )

reuses pixy's variance machinery (calc_tajima_d_stdev, _harmonic_sum) with the
per-class site count replaced by the per-class mutation count, so it reduces to
biallelic pixy when every k_i = 2 (eta_n == s_n). raw_pi is pixy's multiallelic pi.
"""

from __future__ import annotations

from collections import Counter
from typing import Dict, Optional

import allel
import numpy as np
from allel import GenotypeArray

from pixy.calc import _harmonic_sum, calc_tajima_d_stdev


def eta_stats(gt: GenotypeArray) -> Optional[Dict[str, float]]:
    """mutation-count theta_W and tajima's d for one genotype array."""
    ac = gt.count_alleles()
    num_sites = int(np.count_nonzero(ac.sum(axis=1)))
    if num_sites == 0:
        return None

    raw_pi = float(np.sum(allel.mean_pairwise_difference(ac=ac, fill=0)))

    # variant sites only; per site n_i = observed chromosomes, k_i = observed alleles
    variant = ac[ac[:, 1:].sum(axis=1) != 0]
    eta_counts: Counter = Counter()  # n -> total mutations (sum of k-1)
    for row in variant:
        n_i = int(row.sum())
        k_i = int(np.count_nonzero(row))
        if n_i >= 2 and k_i >= 2:
            eta_counts[n_i] += (k_i - 1)

    watterson_eta = 0.0
    for n_i, eta_n in eta_counts.items():
        a1 = float(_harmonic_sum(n_i))
        if a1 > 0:
            watterson_eta += eta_n / a1

    # same denominator as pixy, mutation count in place of site count. since #160 this
    # is one evaluation at the pooled mean n, not a per-class sum (see calc_tajima_d_stdev).
    # vcfsim masks a fixed count of individuals per site, so n is constant and the two
    # forms coincide exactly
    per_site_n = ac.sum(axis=1)
    d_stdev_eta = calc_tajima_d_stdev(
        total_allele_count=int(per_site_n[per_site_n > 0].sum()),
        num_sites=num_sites,
        num_mutations=int(sum(eta_counts.values())),
    )
    if d_stdev_eta > 0 and np.isfinite(raw_pi) and np.isfinite(watterson_eta):
        d_eta = (raw_pi - watterson_eta) / d_stdev_eta
    else:
        d_eta = float("nan")

    return {
        "num_sites": float(num_sites),
        "raw_pi": raw_pi,
        "watterson_eta": watterson_eta,
        "thetaw_eta": watterson_eta / num_sites,
        "tajimaD_eta": d_eta,
    }
