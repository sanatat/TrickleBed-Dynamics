# v1.6: multi-flow identifiability audit

All observations are **synthetic**, generated with the same tracer model family used for estimation. This is not real-data validation.

Reference synthetic values: exchange=0.0025 (model exchange units), axial dispersion=1.4e-05 m²/s. Gaussian noise SD=8e-05.

Single-flow fit: [0.002527764655549253, 1.0000001337744467e-07]; Jacobian condition: 2.23e+03, local parameter correlation: 0.986.

Two-flow fit: [0.0023215395552692395, 1.4135506959747976e-06]; Jacobian condition: 129, local parameter correlation: 0.925.

Holdout RMSE (flow +60%): single 9.022539e-05, multi 8.024873e-05, weak baseline 0.0008893682.

## Interpretation

A lower local Jacobian condition indicates more balanced information for the two fitted log-parameters, but does not prove structural identifiability. The held-out flow uses the same generative model, so this tests numerical generalization under idealized assumptions, **not** physical predictiveness. Do not claim calibrated Ni/Al2O3 kinetics or published-data validation.

The forward RTD model and the reactive solver share a dispersion operator, but this script does not infer reaction kinetics or validate thermal behavior.
