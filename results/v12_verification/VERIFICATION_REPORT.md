# v1.2 targeted numerical verification

**Model:** Uncalibrated one-dimensional four-compartment finite-volume prototype. NOT experimentally validated.

## Improvements

- Three independent BDF relative/absolute tolerance settings, explicitly splitting integration at the imposed feed step.
- All-compartment carbon- and hydrogen-equivalent differential balances checked against external boundary flux.
- Explicit energy RHS identity and diagnostics for positivity and pore-saturation bounds.

## Results

- Maximum equivalent flux budget residual: **5.551e-17** (concentration×velocity units; per unit frontal area).
- Max relative outlet variation vs tightest tolerance: **6.265e-06**.
- Lowest species inventory encountered at final state among runs: **0.0001219** mol per m³ bed.

The finite-volume inventory flux identity checks implementation consistency; it is not a fit to physical experiment. Energy RHS comparison only checks the written equation against its derivative, not the accuracy of thermal constitutive laws. The physical model still lacks heat capacity variation, latent heats, pressure closure, genuine wetting correlations, two-isomer kinetics, measured calibration, and 2012 Figure 2 quantitative match.
