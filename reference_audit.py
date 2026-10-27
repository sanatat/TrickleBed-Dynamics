"""Source-grounded reference inlet / energy-gap audit, not CFD reproduction.

Only Table 2 of Gorshkova et al (2012) is used for reactor inputs.
Molar masses are standard chemical formula values, not fitted parameters.
"""
from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = {
  "source":"Gorshkova et al. (2012), Table 2, test 1 (printed p.402)",
  "bed_height_m":1.5,"bed_width_m":1.0,"pellet_diameter_m":0.003,
  "pellet_porosity":0.50,"bed_porosity":0.37,
  "catalyst_mass_kg":1876.0,"pore_diameter_m":1e-5,
  "pressure_bar":20.0,"inlet_T_K":400.0,
  "liquid_mass_flow_kg_s":1.042,"gas_mass_flow_kg_s":0.02071,
  "liquid_mass_fractions":{"octene":0.3,"octane":1e-5,"H2":0.000377},
  "gas_mass_fractions":{"octene":0.526,"octane":0.0,"H2":0.407},
}
# Formula-based standard molar mass approximation; 2012 paper uses isomers.
MW_KG_MOL={"octene":0.112215,"octane":0.114231,"H2":0.00201588,"dodecane":0.170339}

def complete_mass_fractions(partial):
    d=dict(partial)
    remainder=1.-sum(d.values())
    if remainder < -1e-9: raise ValueError('Mass fractions exceed 1')
    d['dodecane']=remainder
    return d

def component_molar_flows(total_mass_flow, partial):
    fractions=complete_mass_fractions(partial)
    return {k:total_mass_flow*x/MW_KG_MOL[k] for k,x in fractions.items()}

def audit():
    l=component_molar_flows(DATA['liquid_mass_flow_kg_s'],DATA['liquid_mass_fractions'])
    g=component_molar_flows(DATA['gas_mass_flow_kg_s'],DATA['gas_mass_fractions'])
    combined={k:l[k]+g[k] for k in MW_KG_MOL}
    # Stoichiometric upper bound for C8H16+H2 -> C8H18 on basis of inlet molar flow.
    max_reaction_mol_s=min(combined['octene'],combined['H2'])
    return {"reference":DATA,"molecular_weights_kg_mol":MW_KG_MOL,
        "liquid_feed_mol_s":l,"gas_feed_mol_s":g,"total_feed_mol_s":combined,
        "stoichiometric_upper_bound_reaction_mol_s":max_reaction_mol_s,
        "limiting_inlet_species":'H2' if combined['H2']<combined['octene'] else 'octene',
        "theoretical_max_octene_conversion_at_inlet_stoichiometry":max_reaction_mol_s/combined['octene'],
        "assumptions":["Unspecified mass fraction assigned to inert dodecane, as named in paper.","Constant formula-based molecular weights.","Upper bound ignores kinetics, phase partitioning, evaporation and recycle."]}

def write_report(out=Path('results/v8_audit')):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=audit()
    (out/'reference_feed_audit.json').write_text(json.dumps(d,indent=2))
    rows=json.loads(Path('results/v7_validation/model_metrics.json').read_text())
    with (out/'reference_vs_synthetic.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['chi','model_peak_T_K','model_dT_K','paper_approx_peak_T_K','reference_uncertainty_note']);writer.writeheader()
        for row in rows:
            writer.writerow({'chi':row['chi'],'model_peak_T_K':row['model_peak_T_K'],
                'model_dT_K':row['model_max_temp_rise_K'],
                'paper_approx_peak_T_K':500,
                'reference_uncertainty_note':'Approximate narrative statement p402; not digitized'})
    x=np.arange(len(rows)); fig,ax=plt.subplots(figsize=(7,4.5))
    ax.bar(x-.18,[r['model_max_temp_rise_K'] for r in rows],width=.36,label='v0.7 model (illustrative)')
    ax.bar(x+.18,[100]*len(rows),width=.36,label='2012 article, approximate')
    ax.set_xticks(x,[f"chi={r['chi']:g}" for r in rows]);ax.set_ylabel('Maximum temperature rise (K)')
    ax.set_title('Reference gap: not a quantitative CFD validation');ax.legend();fig.tight_layout()
    fig.savefig(out/'thermal_gap.png',dpi=175);plt.close(fig)
    report=f"""# v0.8 source-grounded input audit

**Scientific status: NOT VALIDATED; NOT A REPRODUCTION.** This release deliberately does not tune synthetic kinetic constants to resemble the paper. It makes the original inlet mass flows and compositions machine-readable and audits units/stoichiometric bounds.

## Verified inputs

Gorshkova et al. 2012, **Table 2 / Test 1**: 1.5 m height, 1 m width (test geometry), 3 mm pellet size, bed void fraction 0.37, particle pore fraction 0.5, catalyst mass 1876 kg, 20 bar, 400 K, liquid flow 1.042 kg/s, gas flow 0.02071 kg/s. The original table lists gas octene mass fraction **0.526**, H2 **0.407**, not pure H2 gas. Unlisted fractions are assigned to inert dodecane, as stated in the paper. Molar masses are calculated from chemical formulas, with values in `reference_audit.py`.

### Calculated stoichiometric feasibility (feed-only, not predicted conversion)

- Total octene inlet: **{d['total_feed_mol_s']['octene']:.5g} mol/s**.
- Total H2 inlet: **{d['total_feed_mol_s']['H2']:.5g} mol/s**.
- Limiting inlet species: **{d['limiting_inlet_species']}**.
- Absolute upper bound at these feed flows: **{100*d['theoretical_max_octene_conversion_at_inlet_stoichiometry']:.2f}%** of total C8 alkene feed, assuming all reactant H2 can be used. This is *not* kinetic conversion or an experimental observation.

### Thermal comparison

Article Test 1 says ~400 → ~500 K (printed p.402); v0.7 synthetic temperature rises are much smaller. See `reference_vs_synthetic.csv` and `thermal_gap.png`. There is no error score against digitized Figure 2 data because the figure has not been numerically digitized and the model plots concentrations rather than the published molar fractions.

### Outstanding requirements for a physically faithful model

1. Verified numeric reaction and adsorption coefficients, temperature dependence and units from Lylykangas et al. (2003). We located a full-text web search hit, but could **not reliably retrieve and inspect its tabulated parameter values**. Values remain unknown, not guessed.
2. Full Lappalainen et al. (2008) empirical wetting closure, including dimensionless groups and range of validity; publisher abstract alone is insufficient.
3. Physically consistent phase equilibrium, latent heat, H2 solubility, reaction heat, heat capacities, pressure and flow-based holdup closures.
4. Separate two-isomer reaction chemistry if reproducing experimental rate parameters in detail.
5. Faithful conversion of inlet flows into model boundary conditions, real figure digitization, mesh/time-step convergence and literature/experimental error bars.

**Applicability:** auditable research-training progress and input sanity checks. Existing solver remains illustrative and should NOT be described as a published-model reproduction.
"""
    (out/'AUDIT_REPORT.md').write_text(report)
    print(report.split('### Thermal comparison')[0])

if __name__=='__main__':write_report()
