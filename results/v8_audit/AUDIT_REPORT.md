# v0.8 source-grounded input audit

**Scientific status: NOT VALIDATED; NOT A REPRODUCTION.** This release deliberately does not tune synthetic kinetic constants to resemble the paper. It makes the original inlet mass flows and compositions machine-readable and audits units/stoichiometric bounds.

## Verified inputs

Gorshkova et al. 2012, **Table 2 / Test 1**: 1.5 m height, 1 m width (test geometry), 3 mm pellet size, bed void fraction 0.37, particle pore fraction 0.5, catalyst mass 1876 kg, 20 bar, 400 K, liquid flow 1.042 kg/s, gas flow 0.02071 kg/s. The original table lists gas octene mass fraction **0.526**, H2 **0.407**, not pure H2 gas. Unlisted fractions are assigned to inert dodecane, as stated in the paper. Molar masses are calculated from chemical formulas, with values in `reference_audit.py`.

### Calculated stoichiometric feasibility (feed-only, not predicted conversion)

- Total octene inlet: **2.8828 mol/s**.
- Total H2 inlet: **4.3762 mol/s**.
- Limiting inlet species: **octene**.
- Absolute upper bound at these feed flows: **100.00%** of total C8 alkene feed, assuming all reactant H2 can be used. This is *not* kinetic conversion or an experimental observation.

### Thermal comparison

Article Test 1 says ~400 → ~500 K (printed p.402); v0.7 synthetic temperature rises are much smaller. See `reference_vs_synthetic.csv` and `thermal_gap.png`. There is no error score against digitized Figure 2 data because the figure has not been numerically digitized and the model plots concentrations rather than the published molar fractions.

### Outstanding requirements for a physically faithful model

1. Verified numeric reaction and adsorption coefficients, temperature dependence and units from Lylykangas et al. (2003). We located a full-text web search hit, but could **not reliably retrieve and inspect its tabulated parameter values**. Values remain unknown, not guessed.
2. Full Lappalainen et al. (2008) empirical wetting closure, including dimensionless groups and range of validity; publisher abstract alone is insufficient.
3. Physically consistent phase equilibrium, latent heat, H2 solubility, reaction heat, heat capacities, pressure and flow-based holdup closures.
4. Separate two-isomer reaction chemistry if reproducing experimental rate parameters in detail.
5. Faithful conversion of inlet flows into model boundary conditions, real figure digitization, mesh/time-step convergence and literature/experimental error bars.

**Applicability:** auditable research-training progress and input sanity checks. Existing solver remains illustrative and should NOT be described as a published-model reproduction.
