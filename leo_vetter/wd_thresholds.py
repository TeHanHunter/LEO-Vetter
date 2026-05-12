"""LEO-Vetter thresholds and per-test overrides tuned for white-dwarf hosts.

Drop-in companion to leo_vetter.thresholds, kept fully additive: the FGKM
defaults in thresholds.py are untouched, and this module is opt-in via
``check_thresholds_wd(metrics, case)``.

White-dwarf hosts violate two assumptions baked into the standard LEO-Vetter
test suite:

1. **R_planet << R_star.** For WDs, R_WD ~ 0.005-0.018 R_sun ~ 0.4-1.4 R_jup,
   so any companion larger than Mars-sized has R_p >= R_*. The transit
   chord scales with R_p + R_*, not R_* alone. WD 1856b has
   R_p / R_WD ~ 7-8 and is itself a published planet at b ~ 0.78. This breaks:
     - the V-shape test, which uses (b + R_p/R_*) > 1.5 (FGKM-grazing only
       satisfies this when *abnormally* grazing); for WDs the typical
       value is ~5-15 and the test fires for every plausible candidate
     - the stellar-density expected-qtran formula in parameters.get_q,
       which derives qtran for a point-mass companion crossing R_* alone
2. **TESS 200 s cadence smears short transits.** WD transits are 5-20 min,
   so ingress + egress get blurred into a near-V shape regardless of true
   underlying geometry, and BLS-fitted depth is biased low for any transit
   with fewer than ~5 cadences across the dip. Any threshold that assumes
   a clean trapezoidal U-shape or trusts the absolute depth value should
   be loosened.

Per-target host parameters (M_WD, R_WD from Gaia parallax + Teff) can be
plugged into the `star` dict that LEO-Vetter already accepts; this file
provides the threshold layer on top.

References:
  Kunimoto et al. 2025, AJ 170, 280 — LEO-Vetter
  Vanderburg et al. 2020, Nature 585, 363 — WD 1856 b
"""
from __future__ import annotations

import numpy as np

from leo_vetter.thresholds import (
    weak,
    invalid_transits,
    bad_shape,
    non_unique,
    chases,
    dmm,
    single_event,
    bad_fit,
    sinusoidal,
    asymmetric,
    chi,
    data_gapped,
    odd_even,
    large,
    secondary,
    offset,
)


_R_JUP_PER_R_SUN = 0.10054  # 1 R_jup / 1 R_sun


wd_thresholds = {
    # ---- inherited from FGKM defaults, kept identical ----
    "MES": 6.2,
    "SHP": 0.6,
    "MS1": 0.2,
    "MS2": 0.8,
    "MS3": 0.8,
    "chases": 0.78,
    "DMM": 1.5,
    "max_SES_to_MES": 0.8,
    "AIC1": -60,
    "AIC2": -30,
    "SWEET": 15,
    "ASYM": 8,
    "CHI": 7.8,
    "frac_gap": 0.5,
    "MS4": 0,
    "MS5": -1,
    "MS6": -1,
    "offset": 15,

    # ---- WD-specific changes ----

    # Long-period WD planet candidates only get ~2-3 transits in one 27 d
    # TESS sector (P=10 d → 2 transits). Drop to 2 so we don't reject those.
    # Multi-sector merge restores statistical confidence later.
    "N_transit": 2,

    # V-shape test effectively disabled. For WD planets R_p / R_* is
    # typically 5-15 and (b + R_p/R_*) > 1.5 fires unconditionally; this
    # is a property of WD geometry, not a vetting signal. The 200 s
    # cadence smear also makes the observed shape near-V even for clean
    # box transits, so we cannot trust U-vs-V discrimination at this
    # cadence. Re-enable for higher-cadence follow-up data only.
    "V_shape": 1.0e6,

    # Planet/BD boundary at ~2 R_jup ≈ 22.4 R_earth. Default LEO threshold
    # of 22 R_earth coincidentally matches our boundary; explicit here.
    "size": 22.4,

    # ---- new WD-specific parameters ----

    # Maximum plausible companion radius in R_jup for chord-crossing
    # envelope (the 'fork-the-physics' part of unphysical_duration_wd).
    # Above ~2 R_jup we are in the inflated-BD or stellar regime and the
    # candidate belongs to a different classification class anyway.
    "R_comp_max_R_jup": 2.0,

    # The unphysical_duration FA check compares observed q to expected q
    # for the largest plausible companion at b=0. Loosened from FGKM
    # default 0.6 to 0.4 to allow grazing WD transits (b > 0 shortens
    # the observed q relative to the chord-crossing maximum).
    "q_ratio_min": 0.4,
}


def vshaped_wd(metrics, thresholds):
    """Disabled for WD hosts (see module docstring). Always returns False."""
    message = "FP: V-shape (disabled for WD hosts)"
    if isinstance(metrics, dict):
        return False, message
    return np.zeros(len(metrics), dtype=bool), message


def _qtran_max_wd_chord(metrics, thresholds):
    """Expected maximum observable qtran (= dur/period) for a WD host at b=0
    with the largest plausible companion (R_comp_max = ``thresholds['R_comp_max_R_jup']`` R_jup).

    Derived from T_dur(b=0) = 2 (R_* + R_p_max) / v_orb, where v_orb = 2 π a / P:

        qtran_max = T_dur / P = (1 + R_p_max / R_*) / (π × a/R_*)
                              = (1 + R_p_max / R_WD) / (π × aRs)

    Uses metrics['aRs'] = a/R_* already computed by parameters.derived_parameters
    and metrics['Rs'] (in R_sun) — needs to be present in metrics; integration
    adapter must inject it (see twirl/vetting/leo_adapter.py, to be written).
    """
    Rs_rsun = metrics.get("Rs")
    if Rs_rsun is None:
        # Fall back to canonical R_WD if the adapter didn't inject Rs
        Rs_rsun = 0.013
    r_pmax_over_rs = thresholds["R_comp_max_R_jup"] * _R_JUP_PER_R_SUN / Rs_rsun
    aRs = metrics["aRs"]
    return (1.0 + r_pmax_over_rs) / (np.pi * aRs)


def unphysical_duration_wd(metrics, thresholds):
    """WD-adapted version of leo_vetter.thresholds.unphysical_duration.

    Replaces the MS-density expected-qtran (parameters.get_q, which assumes
    a point-like companion crossing R_* alone) with a chord-sum expectation
    that accounts for R_p ~ R_*. Other clauses (low aRs, trap_qtran > 0.5)
    are kept as-is.

    Triggers FA only when the *observed* qtran exceeds what is physically
    possible for any companion in the planet-regime (R_p <= R_comp_max);
    grazing transits with short observed qtran are allowed.
    """
    message = "FA: unphysical transit orbit (WD chord)"
    low_aRs = (metrics["transit_aRs"] < 1.5) | (metrics["aRs"] < 2)
    qtran_max = _qtran_max_wd_chord(metrics, thresholds)
    long_duration = (
        (qtran_max / metrics["trap_qtran"] < thresholds["q_ratio_min"])
        | np.isnan(metrics["sig_sec"])
        | (metrics["trap_qtran"] > 0.5)
    )
    return (low_aRs | long_duration), message


def check_thresholds_wd(metrics, case, verbose=False, thresholds=None):
    """WD-host version of check_thresholds.

    Same test inventory as the FGKM default, with WD-adapted versions for
    ``unphysical_duration`` and ``vshaped`` and the WD thresholds dict.
    """
    if thresholds is None:
        thresholds = wd_thresholds

    if case == "FA":
        tests = [
            weak,
            invalid_transits,
            bad_shape,
            non_unique,
            chases,
            dmm,
            single_event,
            bad_fit,
            sinusoidal,
            unphysical_duration_wd,
            asymmetric,
            chi,
            data_gapped,
        ]
    elif case == "FP":
        tests = [
            odd_even,
            vshaped_wd,
            large,
            secondary,
        ]
        if isinstance(metrics, dict) and "offset_qual" in metrics:
            tests.append(offset)
    else:
        raise ValueError("Case must be FA or FP")

    if isinstance(metrics, dict):
        mask = False
    else:
        mask = np.zeros(len(metrics), dtype=bool)
    for test in tests:
        flag, message = test(metrics, thresholds)
        mask |= flag
        if isinstance(metrics, dict) and verbose and flag:
            print(message)
    if isinstance(metrics, dict) and verbose and not mask:
        print(f"Passed all {case} tests")
    return mask
