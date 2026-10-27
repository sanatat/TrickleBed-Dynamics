# v2.7 — Heterogeneous RTD inspired by Sales et al. (2005)

## Scientific basis

User-supplied Sales et al. (2005), *Experimental Dynamic Evaluation of Three-Phase Reactors*, **printed pages 445–450**, esp. Eqs. **(7)–(9)** for external liquid transport, spherical intrapellet tracer diffusion, adsorption; Eq. (11)/(12) Danckwerts boundaries and internal sphere flux; Tables **2–3** and Figs. **5–6**. We do **not** claim a literal implementation of analytical transfer function Eqs. (13)–(15).

## Computation and comparison

Using Figure 5 vector-extracted graphical traces from v2.6 (QL = **5 L/h**, QG = **300 L/h**), simulate an initially equilibrated negative tracer step with measured inlet forcing, conservative upwind finite volume advection, axial diffusion, wetted external-sphere exchange, radial spherical diffusion and reversible linear adsorption. Analytical adsorption relation is replaced by a first-order relaxation dynamic qdot = k_ads*(K_ads*C-q), whose kinetic coefficient is **assumed**, not claimed from the paper.

| Variant | Number tuned parameters | Graphical RMSE |
|---|---:|---:|
| Axial-only fixed parameters | 0 | 0.014869 |
| Heterogeneous fixed literature/assumed | 0 | 0.049381 |
| Heterogeneous + one fitted axial dispersion multiplier | 1 | 0.049309 |
| Gamma RTD (v2.6) | 3 | 0.007039 |

One-parameter fit: D_ax multiplier = **1.2954**; reached optimization bound: **False**. Lower graph RMSE is not alone proof of improved physics. These are in-sample graphical fits and have unequal numbers of calibrated parameters.

Conservative stationary check max |A*x_equilibrium+b| = **4.44e-16 s^-1**. The **global inventory identity** w@A + outlet_flux = 0 yields maximum residual **2.78e-17 m/s** and inlet residual **0 m/s** for any state. This certifies discrete tracer inventory bookkeeping, not the physical realism of the closures.

## Explicit provenance of parameters

From Table 3: hL=0.17, fe=0.38, Pe=4.84; uL=0.001 m/s, for the low gas-flow row. Table 2 gives k_LS **(6.11–11.81)e-5 m/s**, choose **8e-5 m/s** *within the range*, not explicitly for one Figure 5 run. Page 450 text states Deff = **7.65e-10 m²/s**, whereas **Table 2 displays a different range in 10^-7 units**; these should not be reconciled without original author clarification. We therefore perform a two-order-of-magnitude sensitivity check, not claim either is universally correct. Text gives Kads = 8.20e-5 m³/kg. The adsorption kinetics rate 0.03/s and packed-bed porosity 0.40 are **assumptions**. Other reported dimensions: 0.8m length, particle diameter 4–6mm; choose midpoint 5mm.

The paper's Pe definition Pe=Z*u/(D_ax*epsilon) is applied as an effective bulk flux coefficient using **assumed epsilon**, not as universal dispersion calibration. This discretization and boundary implementation differ from the authors' inverse Laplace method. Film area and interior liquid volume are taken proportional to fe, which is a simplifying closure rather than the article's exact microstructure treatment.

## Why Figure 6 is not claimed as independent validation

Published Figure 6 (QL = 10 L/h) has extracted **experimental outlet markers** but **no paired numerical inlet signal**. Choosing an unmeasured inlet waveform to maximize output agreement would be circular. Thus our v2.7 report is a rigorous **model comparison for Figure 5**, not a fully validated second-flow predictive test. Likewise it is **not** validation of the reactive or exothermic 2012 Alopaeus model.

## Sensitivity

See `parameter_sensitivity.csv`, which varies assumed bed porosity and the disputed Deff scale. The resulting changes are diagnostic *model uncertainty*, not confidence intervals.

## Files and run

- `heterogeneous_rtd_v27.py` — sparse FV physical model, figures/metrics generator.
- `figure5_models.csv` — vector-derived measured figure curves and predictions.
- `figure5_heterogeneous_comparison.png` — graphic.
- `metrics.json` — reproducible results, explicit assumptions and limitations.
- `parameter_sensitivity.csv` — physical assumption stress test.

```bash
python heterogeneous_rtd_v27.py
python -m pytest -q
```

No copied source PDF included in archive; raw source PDF is held separately by the user.
