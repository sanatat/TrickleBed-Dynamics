# Sales 2005 Figure 5 transfer-function study (v2.1)

**Source:** Sales et al. (2005), DOI 10.1590/S0104-66322005000300012, Fig. 5 (trickle bed, QL=5 L/h, QG=300 L/h).

**Data status:** CSV inherited from v2.0 consists of 13 manually approximated points from a printed graph. NOT original data, NOT calibrated WebPlotDigitizer points. Figure 6 uses a different flow (QL=10 L/h) and is not a held-out validation of these parameters.

## Method

Starting at normalized steady concentration 1, the measured inlet step-down is treated as an exogenous input. Its negative numerical derivative is convolved with a delayed gamma residence-time kernel to predict the outlet response. Gamma CDF bin differences enforce kernel non-negativity and unit integral. Fit uses scipy.optimize.least_squares.

## Fit to approximate points

- Gamma mean residence time: 163.241 s
- Gamma shape: 6.497
- Additional delay: 0.000 s
- RMSE on approximate points: 0.00713
- Zero-residence (outlet=inlet) RMSE: 0.23683

**Not physical validation:** These estimates are highly sensitive to hand-reading, baseline, image quality and model misspecification. The 35 perturbations are a stress test, NOT confidence intervals. This reduced-order residence-time kernel is distinct from the published multi-parameter diffusion/adsorption model and the earlier reactive four-compartment model.

## Next requirements

Digitize original figure using pixel-calibrated x/y axes; report point selection and extraction errors. Fit directly to those points, compare competing residence-time kernels, and test across another condition with physically justified flow scaling. Full reactor validation also requires thermal, kinetic and phase data.
