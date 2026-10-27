# v2.2 — Temporal holdout and model comparison

**Source**: Sales et al. (2005), Figure 5, DOI: 10.1590/S0104-66322005000300012.

**Critical provenance**: Data are 13 *approximate manually read points* inherited from v2.0, not author-supplied raw observations or pixel-calibrated coordinates. The assumed reading errors are not metrological uncertainty estimates.

## What was tested

- Predict outlet from approximate measured inlet by convolution with normalized residence-time kernel.
- Fit only outlet points through 750 s (9 points), and evaluate on four later points **from the same experimental curve**.
- Compare flexible gamma RTD and two-parameter delayed exponential RTD. This is a temporal split, not independent validation on another reactor or flow rate.

## Results

| Model | Training RMSE | Temporal holdout RMSE |
|---|---:|---:|
| Gamma | 0.004756 | 0.010685 |
| Exponential | 0.010205 | 0.010546 |
| No residence time | 0.284278 | 0.021225 |

**Interpretation limitations**: the inlet is known through the held-out interval, and each model sees early values of the same decay curve. With only four test points and uncertain digitization, these scores do not establish generalization or reliable physical parameter values. The models are empirical RTD kernels, not full CFD or the reactive catalyst-pore model.

## Next steps needed for experimental validation

1. Acquire high-resolution source figure and digitize all distinct measurement markers with pixel-to-axis calibration, archiving the selected coordinates and image.
2. Obtain raw data from authors if possible.
3. Measure/obtain separate inlet and outlet traces at a second flow rate; fit one condition and predict the other with a physically supported flow-scaling law.
4. Report RTD normalization, tracer recovery, instrument lag, and identifiability diagnostics.
