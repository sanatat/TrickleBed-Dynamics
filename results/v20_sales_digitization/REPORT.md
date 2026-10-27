# Sales (2005) published-figure extraction — v2.0 feasibility audit

**Important provenance:** Approximate values manually read from the low-resolution journal figure shown on PDF page 8 (printed p.450). **Not raw data, not calibrated high-resolution WebPlotDigitizer extraction and not independent experimental validation.** Numerical precision of CSV time values is not evidence of measurement precision.

Source: https://doi.org/10.1590/S0104-66322005000300012 (Figure 5: liquid 5 L/h; Figure 6: liquid 10 L/h; gas 300 L/h). Page 8 of the PDF displays these figures and Table 3. Figure 5 shows experimental inlet and outlet; Figure 6 shows experimental response versus authors' model.

## Approximate digitization and feasibility result

- Figure 5: 13 manually chosen approximate points per curve; Figure 6: 19 approximate points.
- Coordinates were read from a displayed page image without calibration against high-resolution image pixel positions. Conservative subjective reading uncertainty ~±45 s and ±0.06 normalized fraction; errors could be larger in overlapping lines and the tails.
- Figure 6 was fitted to the empirical shifted-gamma survival function `Q(k, max(t-delay,0)/scale)`, with `scipy.special.gammaincc`.
- Fit delay 136.4 s, shape 1.66, scale 89.1 s; in-sample RMSE to **approximate hand-read points** 0.0040 (unitless). These are not identified reactor parameters, nor a fit of the v1.8/v1.9 multiphase PDE to original measured RTD.
- Figure 5 and Figure 6 correspond to **different liquid flow rates**; treat them as distinct experiments. Figure 5 inlet is not an instantaneous pulse, so comparing normalized decreases as impulse E(t) without deconvolution would be incorrect.
- No numeric raw time-series files are supplied in the published article; Table 3 values are author-inferred model parameters, not raw observations.

## Next scientific action

1. Obtain original PDF image at high resolution; calibrate axes with WebPlotDigitizer, digitize individual experimental *markers*, especially Figure 6; export coordinates and perform manual spot checks.
2. Request the authors' original time/conductivity readings, including inlet and outlet.
3. Implement input-output transfer response for the model using measured inlet shape and account for tracer adsorption and liquid-solid exchange. Use Figure 5 for one condition, Figure 6 for a second, without pooling flow-rate-specific parameters.
4. Quantify sensitivity to pixel extraction, adsorption, instrument delay, and baseline; reserve another experiment for holdout validation.

## Reproduce

`python sales2005_figure_audit_v20.py` (reads approx CSV, refits and writes diagnostic plot).

**Scientific status:** exploratory published-figure digitization feasibility, not full published-data quantitative validation.
