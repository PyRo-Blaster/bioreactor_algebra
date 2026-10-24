# bioreactor_algebra

Design notes for setting CHO fed-batch culture parameters on new 7 L and 15 L bioreactors, using multivariable optimization: scaled input → response maps, SVD of the sensitivity matrix J (stiff and sloppy modes), and canonical analysis of the quadratic surface B.

## Contents

| File | What it is |
| --- | --- |
| [designs/01_original_design_svd.md](designs/01_original_design_svd.md) | v1 design: DSD with vessel as a factor (~21 runs), forward selection, bootstrap CIs, SVD-based PAR/NOR |
| [designs/02_alternative_design_bayesian.md](designs/02_alternative_design_bayesian.md) | Alternative: 13-run DSD + 3-run vessel bridge, horseshoe shrinkage, CVaR setpoint, posterior PAR/NOR |
| [designs/03_review_of_alternative.md](designs/03_review_of_alternative.md) | Mathematical review of the alternative, including defects shared by both designs |
| [scripts/dsd_check.py](scripts/dsd_check.py) | Checks for the review: DSD rank and aliasing, estimability at the corners, pass-gate arithmetic |

Live editable v1 with figures: [Claude Doc](https://claude.ai/artifact/MEE9WVGSyKkfqxD31KCbKT).

## Next steps

- Merge the two designs into a corrected v2.
- Python pipeline: JSON schema → scaling → fit → SVD / canonical analysis → decisions, tested on synthetic data with known modes.
