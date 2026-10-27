# v1.7 — Local optimal design for synthetic tracer experiments

This is a reproducible **in-silico demonstration**, not experimental validation.

Candidate flow multipliers: (0.65, 0.8, 1.0, 1.2, 1.35, 1.6, 1.9)
Pre-specified baseline: (1.0, 1.35)
Optimal pair: (1.6, 1.9)
Baseline log det FIM: 11.0724
Optimal log det FIM: 15.0443
Increase: 3.9719
Fitted baseline log10 parameter conditioning: 13.192
Fitted optimal log10 parameter conditioning: 15.748
Holdout RMSE, baseline: 9.4176516e-05
Holdout RMSE, optimal: 7.7574813e-05

## Limitations
- The design objective uses sensitivities around assumed true parameters; real-world nominal parameters are unknown.
- Fisher information is a local approximation that depends on sampling schedule and a known constant noise assumption.
- Testing uses one synthetic noise realization and the same model family used to generate the data.
- No external physical RTD data, catalyst kinetic calibration, phase equilibrium, or industrial thermal validation.
- A larger information determinant does not guarantee lower hold-out error in a given noisy experiment.
- The current flow multipliers are relative to a demonstration base superficial velocity; feasibility and safe operation are not evaluated.
