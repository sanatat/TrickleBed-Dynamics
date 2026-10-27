# v1.1 scientific credibility audit

**Status: numerical demonstrator, not validated against experiment or the 2012 CFD model.**

## What was measured in code

- In a closed system the discrete source terms conserve alkene-equivalents A+P and hydrogen-equivalents H+P. Maximum differential residual across all states tested: **5.898e-17 mol m⁻³ s⁻¹**. This is a code-internal check only.
- In a short transient run with prescribed flow perturbation, changing the *assumed* kinetic scale, pore-transfer coefficient, starting wetness and effective diffusion changes the predicted peak temperature. See `sensitivity.csv`. This is an identifiability warning, not fitted parameter estimates.
- Mesh diagnostic: 12 versus 18 axial cells give an outlet-alkene difference of **0.05024 mol m⁻³** (0.02%). A three-grid comparison by itself is not a verification of grid independence.

## Research-relevant unresolved gaps

1. Full source kinetic constants from Lylykangas et al., including two distinct octene isomers and their isomerization. The existing single-lumped-species kinetic coefficient is not published data.
2. Independently validated wetting hysteresis / capillarity closure and pressure-drop coupling.
3. Vapor-liquid flash and latent heat; temperature-dependent thermophysical properties.
4. Complete species inventory including inert dodecane, accurate gas/liquid volume balances, and pressure equations.
5. Reproducible CFD/profile reference for quantitative comparison, not just visually read paper figures.
6. Stiffness/parameter-identifiability studies with experimental data before making thermal-safety claims.

## Interpretation

A large increase in model test count does not constitute scientific validation. The 2012 work uses transient Eulerian CFD; this project remains a 1D compartment surrogate. Figure 2 of the publication has **not** been quantitatively reproduced. The model uses arbitrary k and heat parameters and must not be used for process design, safety decisions, or industrial predictions.
