# v2.8 — Physical discrepancy diagnosis, Sales (2005) Figure 5

## What is actually observed

The uploaded Sales et al. (2005) PDF contains experiments. Figure 5 (printed p.450) gives inlet and outlet negative-step tracer traces at 5 L/h liquid. v2.6 extracted **59 graph-resampled points** for each curve from PDF vector paths. These are **not original apparatus time series** and their neighboring points are not independent measurements.

## Root-cause study, not retrofit confirmation

We held literature-linked values (hL=0.17, f_e=0.38, Pe=4.84, D_eff=7.65e-10 m²/s from paper prose, k_LS chosen within the reported range) constant unless the ablation specifies otherwise. The published article does not supply the pore-accessibility fraction, the adsorption relaxation constant, and the packed-bed porosity used by our numerical model. The model is motivated by **Eqs. 7–9**, not an exact implementation of **Eqs. 13–15**.

Fixed-choice ablation results (do not treat these as equally calibrated models):

| Scenario | Graphical RMSE | Model input-output 50% lag (s) |
|---|---:|---:|
| axial_only_reference | 0.01487 | 144.4 |
| full_heterogeneous_reference | 0.04938 | 208.2 |
| no_adsorption | 0.03897 | 196.5 |
| low_pore_access_0p25 | 0.01342 | 160.5 |
| high_intraparticle_diffusion_table2 | 0.05368 | 217.7 |
| low_surface_exchange_0p1 | 0.04868 | 205.0 |
| fast_adsorption_relaxation | 0.04946 | 208.5 |

Observed graph-derived inlet/output 50% lag: **151.0 s**.

## Intraparticle storage accounting (specific to model assumptions)

Chosen mobile capacity hL = **0.1700 m³/m³ bed**. With assumed bed porosity 0.4, f_e 0.38, particle porosity 0.33, adsorption density 1100.0 kg/m³ and K_ads 8.2e-05 m³/kg, **assumed** fully-accessible inner storage capacity = **0.0958 m³ equivalent/m³ bed**, i.e., **56.4%** of mobile storage. In a hypothetical instantaneous-equilibrium limit these capacities would yield an advective lag of approximately 213 s rather than 136 s for mobile liquid alone. These analytical times are *not fitted experimental RTD means*. Their difference pinpoints a plausible excess-storage mechanism in our assumed closures, not a proven experimental root cause.

## Diagnostic single-parameter fitting (NOT validation)

Varying an **assumed** additional pore accessibility factor (between 0 and 1) changes the amount of pore storage and transfer area. The best member of our coarse tested grid is **0.10**, RMSE **0.01233**. This is an *in-sample* diagnostic parameter, not physically identified wetting or porosity. Only 7 trial values were evaluated, not a general optimum or an uncertainty interval.

## Temporal holdout diagnostic (still same experimental trace)

Coarse-grid accessibility chosen using only the transition **300–720 s**: **0.10**, training RMSE **0.01552**; RMSE on subsequent **750–1050 s** points **0.00533**. This is a *within-curve extrapolation* diagnostic, not an independent run nor instrument-level external validation. The Figure-5 input trace remains known over the whole interval. Correlation between graph-resampled points reduces statistical independence.

## Accessibility–dispersion trade-off

A coarse 5 × 3 two-parameter grid on the **same graph-derived Figure 5** found a lowest sampled RMSE of **0.00949** at assumed accessibility **0.10** and dispersion multiplier **0.5**. This is not physical parameter identification. The grid is designed to expose how model mismatch may be redistributed among uncertain processes; it is not a global optimum or experimental uncertainty distribution.

## Numerics

Discrete global tracer conservation identity residual **2.78e-17 m/s**, input identity residual **0 m/s**. Conserving the implemented equations does NOT validate closures.

## Limitation and next discriminating experiment

Figure 5 alone does not separately identify true pore accessibility, adsorption capacity, surface exchange and axial dispersion. Published Table 2 and page-450 prose conflict on effective pore diffusivity; resolve with original authors before claiming a unique value. A separate non-adsorbing tracer (KCl is described by the paper) and the adsorption-capable sodium gluconate at the **same flow**, with recorded inlet/outlet, would separate adsorption from hydrodynamic memory better than further tuning the same curve. Figure 6 has no independently measured inlet series in the PDF and cannot validate a Figure-5-derived transfer model without assumptions.

**This is an evidence-based model-discrepancy study, not experimental confirmation of the full reactive Alopaeus CFD model.**
