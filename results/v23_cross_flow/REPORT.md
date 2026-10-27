# v2.3 — Cross-flow RTD assessment (no refitting)

Source: Sales et al. (2005), Table 3 and Figs 5–6. DOI: 10.1590/S0104-66322005000300012.

## Published values vs assumptions
- **Published fitted experimental-derived Table 3**: at gas velocity 0.0213 m/s, liquid velocity 0.001/0.002 m/s corresponds to holdup 0.17/0.21.
- Bed height 0.80 m (published apparatus). Nominal hydraulic times `height × holdup / superficial liquid velocity`: 136.0 s (5 L/h) and 84.0 s (10 L/h), ratio 0.6176. These are idealized estimates, not measured tracer means.
- v2.2 Figure 5 Gamma mean 163.21 s; shape 6.51, estimated from 13 **manually approximated** graphical readings. No Figure 6 fitting is performed.
- Figure 6 input assumed to be an ideal step-down at t=0. This assumption cannot be confirmed from digital inlet measurements for that trial. Figure 6 points were manually estimated in v2.0, not pixel-calibrated or raw.

## Held-out flow comparison (approximate, illustrative only)
- Scaled model approximate RMSE = 0.41502 normalized fraction.
- Unscaled baseline approximate RMSE = 0.26441.
- The errors include image-reading and inlet-forcing mismatch, and cannot establish the accuracy of the physical reactor model.

## Next evidence requirements
1. Obtain the author's original time-concentration data or a high-resolution source plot.
2. Digitize each experiment with saved pixel coordinates and axis calibration, and assess digitization repeatability.
3. Get measured inlet signal for Figure 6; do not substitute a step unless documented.
4. Validate cross-flow without refitting and report experimental measurement uncertainty.
