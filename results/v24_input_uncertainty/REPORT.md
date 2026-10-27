# v2.4 — Published-figure input uncertainty audit

Source: Sales, Maranhão, Pereira & Abreu (2005), *Experimental dynamic evaluation of three-phase reactors*, https://doi.org/10.1590/S0104-66322005000300012.

## Findings grounded in publication
- Figure 5 displays experimental **entrance**, experimental **exit** and theoretical **exit** for liquid flow **5 L/h**, gas flow **300 L/h**.
- Figure 6 displays experimental and modeled **outlet** for liquid flow **10 L/h**, gas flow **300 L/h**; it does not provide a separate sampled inlet curve for this run.
- The plotted Figure 6 response is decaying, consistent with a step-down-type experimental protocol. The exact forcing time course is not recoverable from the published graph.
- Published Table 3 reports fitted RTD hydrodynamic parameters, not a downloadable time-series dataset.

## What is computed here
- Retains the Gamma parameters fitted to **manual approximations of Figure 5** in v2.2. No refit on Figure 6.
- Transfers the mean residence scale using the published holdup/velocity relationship.
- Compares hypothesized input start times and gradual ramps to the Figure 6 **manual** readings.
- Creates an envelope spanning assumed timings; **not** statistical, pixel, or experimental uncertainty.

## Key caveats
- No accurate pixel-calibrated digitization has been performed in this release. The source PDF exists on the public journal site but was not available as a container file.
- There are no independently measured Figure 6 inlet signal samples available in our data files.
- Any improved RMSE obtained by choosing an input assumption from Figure 6 is an in-sample sensitivity diagnostic, **not** held-out validation.
- This is a separate RTD surrogate, **not** a validation of the reactive multiphase CFD equations.

## Next action to enable genuine cross-flow validation
Obtain author-supplied entrance and outlet time series (both flow rates) or a high-resolution journal PDF that permits reproducible pixel calibration. Archive the raw pixel-to-axis map, digitization error, and sampling protocol; fit only the 5 L/h run and freeze all parameters before predicting 10 L/h.
