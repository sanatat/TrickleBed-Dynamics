# v0.9: Table 2 feed coupled to 1D solver

**Status: demonstrator; neither quantitatively validated nor CFD reproduction.**

The 2012 Test 1 inlet *mass flows and mass fractions* now define the inlet
reactive species concentrations used in the 1D transport equations, unlike v0.8,
where those values were only audited alongside separate illustrative model inputs.

## Mapping assumptions

- Liquid volumetric flow = published liquid mass flow / *assumed* liquid density
  (600, 700, 800 kg/m3 sensitivity analysis).
- Gas volumetric flow = total inlet gas molar flow × RT/P using ideal gas equation;
  this is an assumption, NOT a vapor-liquid flash calculation.
- Cross-sectional area = assumed 1 m2; source describes a 1-m-wide 2D test geometry,
  not a verified 1 m2 physical cross-section.
- Species inventories: A=C8H16 proxy, H2, P=C8H18 proxy; inert C12 omitted as dynamic species.
- Holdup, transfer coefficients, kinetics, heat capacity, reaction enthalpy and
  history-dependent wetting are still synthetic and not independently calibrated.
- The dynamics are 1D plug-flow compartmental balances, not Eulerian CFD.

## Data provenance

Published reference: Gorshkova et al. (2012), Table 2 / Test 1 in `tirik.pdf`.
The source paper states near 400→500 K in Figure 2 conditions, but **it is invalid
 to compute a reproduction error against this with the current physical closures**.

Lylykangas et al. (2003) full parameter table and Lappalainen et al. (2008)
complete wetting closure have NOT been reliably extracted. Parameter values
inherited from v0.8 remain clearly illustrative. The source publication 2003
uses two isooctene isomers and isomerization; current code only lumps one alkene.

## Interpretation

Outputs characterize numerical response and uncertainty under *specified assumed
parameters*, not actual industrial reactor conversion, thermal safety, or
experimental validation. Any proximity to a published temperature is incidental.
