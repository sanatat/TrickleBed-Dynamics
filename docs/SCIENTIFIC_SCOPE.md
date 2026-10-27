# Scientific scope, reliability and next steps

## Motivation

The project investigates tractable numerical approximations for coupled transport and tracer dynamics in gas–liquid–solid trickle-bed reactors, alongside practical identifiability and experimental-design questions.

## What software tests establish

Passing unit and regression tests verifies reported internal invariants and expected behavior in a finite set of scenarios. It **does not** independently establish physical correctness, fit quality outside the available source figures, or predictive performance on industrial plants.

## Sources of discrepancy

The original literature describes multi-phase hydrodynamics, heat and mass transfer, intra-pellet diffusion, wetting and reaction kinetics. This repository uses multiple reduced-order implementations, some coefficients from the literature and many explicitly uncertain or synthetic closure parameters. Numerical stability is distinct from physical calibration.

In the v2.7 and v2.8 studies, the heterogeneous tracer model's match to Figure 5 depended strongly on assumed wetted pore accessibility and adsorption. That *does not demonstrate* the measured pore fraction is 10% or 20%. Model structure and parameters are not uniquely identifiable from the single Figure 5 curve.

## Stronger validation would require

1. Independent sodium-gluconate and KCl inlet/outlet time series plus adsorption equilibrium and instrument response.
2. Joint parameter inference with realistic observational error and model discrepancy.
3. Pre-registered held-out flow conditions (especially Figure 6 with a measured input waveform) and uncertainty intervals.
4. Separate evaluation of hydrodynamic RTD and heat/reaction prediction before presenting model as predictive for industrial control.
5. Reconciliation of adsorption/diffusion units and published transfer-function details.

The repository contains explorations, not proof of novelty or industrial readiness.
