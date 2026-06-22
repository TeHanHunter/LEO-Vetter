# WD Host Tuning Notes

This branch is a TWIRL reference fork for applying LEO-Vetter to white-dwarf
hosts. It keeps the stock FGKM-dwarf behavior unchanged and adds one opt-in
module:

- `leo_vetter/wd_thresholds.py`

Use `check_thresholds_wd(metrics, case)` instead of
`leo_vetter.thresholds.check_thresholds(metrics, case)` when vetting WD-host
candidates.

## Minimal Usage

```python
from leo_vetter.wd_thresholds import check_thresholds_wd

fa = check_thresholds_wd(tlc.metrics, "FA")
fp = check_thresholds_wd(tlc.metrics, "FP")
is_pc = not fa and not fp
```

For TWIRL S56 testing, the caller also uses the WD metric path and injects a
WD radius into the metrics table before threshold checks:

```python
tlc.compute_flux_metrics(star, cap_b=False)
tlc.metrics["Rs"] = 0.013
```

Replace the fallback `0.013 R_sun` with per-target WD radii when they are
available.

## What Changed

Only `leo_vetter/wd_thresholds.py` was added. The stock
`leo_vetter/thresholds.py` defaults and standard `check_thresholds` entry point
are untouched.

WD threshold preset:

- `N_transit = 2`: allows single-sector long-period WD candidates with only two
  observed transits.
- `V_shape = 1e6`: effectively disables the standard V-shape rejection for WD
  hosts.
- `size = 22.4 R_earth`: makes the large-object boundary explicit at about
  `2 R_jup`.
- `R_comp_max_R_jup = 2.0`: maximum companion radius used by the WD
  chord-duration envelope.
- `q_ratio_min = 0.4`: loosens the observed-vs-maximum duration ratio for
  grazing WD geometries.

WD override functions:

- `vshaped_wd`: returns no rejection. WD planets naturally have
  `R_p / R_star >> 1`, and 200 s TESS cadence smears ingress/egress enough that
  the stock V-shape rule rejects plausible WD candidates.
- `unphysical_duration_wd`: replaces the main-sequence point-companion duration
  relation with a WD chord-sum limit using `R_star + R_comp_max`.
- `invalid_transits_wd`: keeps the pruned transit-count requirement but drops
  the post-pruning MES requirement. The normal `weak` test still checks the
  pre-pruning MES floor.
- `non_unique_wd`: keeps the MS1 and MS2 uniqueness checks but drops the MS3
  positive-feature clause, which is fragile for current WD/TESS detrended light
  curves.

Inherited tests are otherwise unchanged: `weak`, `bad_shape`, `chases`, `dmm`,
`single_event`, `bad_fit`, `sinusoidal`, `asymmetric`, `chi`, `data_gapped`,
`odd_even`, `large`, `secondary`, and optional `offset`.

## Rationale

White-dwarf transit vetting violates two assumptions built into the default
FGKM settings:

1. For WDs, physically interesting companions can be much larger than the host
   radius, so `R_p << R_star` is not valid.
2. TESS 200 s cadence samples 5-20 minute WD transits with only a few cadences,
   so per-event pruning and U-shape/V-shape discrimination are too aggressive.

The WD tuning is therefore intentionally permissive. It is meant as a first-pass
triage layer before human inspection and follow-up validation, not as a final
automated planet decision rule.

## Current Check

This tuning was used in the TWIRL S56 WD-host vetting experiment. It passes the
WD 1856 benchmark as `PC` under the TWIRL S56 test setup and is intended as a
reference implementation for Franklin's independent WD tuning.
