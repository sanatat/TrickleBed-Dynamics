# v1.8 Robust optimal experimental design

**All observations are synthetic.** Neither empirical reactor calibration nor industrial validation is claimed.

Baseline flows: (1.0, 1.35); local flows: (1.6, 1.9); robust flows: (1.6, 1.9).

| Design | Worst-case log det FIM | Mean holdout RMSE |
|---|---:|---:|
| baseline | 5.3617 | 9.7733927e-05 |
| local | 9.6692 | 9.7467676e-05 |
| robust | 9.6692 | 9.7467676e-05 |

Optimization: evaluate the local Fisher matrix on five explicitly selected parameter pairs and maximize the **minimum** log determinant; no test data participates in design selection.
Testing: two parameter settings and one Gaussian noise realization each, with a disjoint holdout flow of 1.75×; refit the two parameters for each simulated training pair.

This demonstrates a computational methodology, not physical predictiveness. Optimized Fisher information is not equivalent to guaranteed lower RMSE.

Important limitations: the 5-point scenario set is not a probability distribution or confidence region; true parameters and noise process are unmeasured; no operational feasibility constraint; one noise trial per scenario gives highly uncertain error summaries; wetting and chemistry are not newly calibrated.
