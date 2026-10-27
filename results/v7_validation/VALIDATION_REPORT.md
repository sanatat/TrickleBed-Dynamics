# v0.7: Axial validation gap assessment (2012 paper)

**Status: NOT quantitatively validated.** The output is model-vs-literature assessment, not a CFD reproduction.

## Paper-reported targets (2012 paper, printed p.402)

- In Test 1, reaction occurs over ~1 m at χ=1 and ~0.25 m at χ=10.
- Temperature rises from 400 K to ~500 K at χ=1; with χ=10 the similar maximum is reached around z=0.3 m.
- Table 2 geometry: 1.5 m height, 3 mm pellet diameter, bed porosity 0.37, pellet porosity 0.5, 20 bar and 400 K.

## What this run actually calculates

- Synthetic four-compartment (mobile gas, mobile liquid, pore gas and pore liquid) reduced-order model.
- 75 axial cells, constant fitted transport proxies, one lumped olefin reaction and a prescribed wetting trajectory.
- Not a temperature-dependent vapour-liquid equilibrium calculation or a resolved hydrodynamic model.
- Uses fixed illustrative inlet concentrations and kinetic parameters, NOT the published mass-fraction input or the fitted Lylykangas constants.
- The 2012 figure uses molar fractions; current plots use mol/m³; direct curve-shape error metrics would be misleading.

## Simulation outputs

| χ | Outlet proxy conversion | Maximum T (K) | ΔT max (K) |
|---|---:|---:|---:|
| 1 | 14.46% | 401.75 | 1.75 |
| 10 | 38.04% | 404.49 | 4.49 |

## Finding

The heat rise is substantially smaller than the ~100 K reported in the paper, and the code does not model all physics needed to measure corresponding axial reaction lengths. **These are material discrepancies, not a validated match.** Matching peak temperature by tuning only a synthetic heat coefficient would be curve-fitting rather than independent validation.

## Missing inputs and next scientific work

1. Numerical kinetic constants for the two isooctene isomers, adsorption constants and temperature dependence in Lylykangas et al. (2003).
2. Hydrodynamic/wetting relation of Lappalainen et al. (2008/2009) including liquid holdup and pressure drop.
3. Phase equilibrium and enthalpy properties for the four-component gas-liquid mixture, evaporation and condensation.
4. Actual 2012 inlet composition conversions and consistent mass, energy and species equations.
5. Reference data digitized from Figure 2 with a documented procedure, followed by grid/time-step sensitivity and numerical error metrics.

**Research claim appropriate for applications:** Implemented and tested a paper-inspired prototype and quantified major gaps with its reference. Do not claim reproduction of Figure 2.
