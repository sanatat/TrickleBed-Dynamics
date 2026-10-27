"""Reference-feed-bound model; thermodynamics, holdups and kinetics NOT calibrated.
2012 Gorshkova Table 2 Test 1 data provided via reference_audit.py.
Density assumption and ideal-gas equation produce approximate superficial velocities.
Carbon-number C8 components are tracked as lumped species; inert dodecane is omitted
from differential states (although included in volumetric-flow calculation).
"""
from dataclasses import replace
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reference_audit import audit
from v6_model import Params, simulate, unpack, volumes

R=8.314462618

def map_reference_feed(liquid_density_kg_m3=700., area_m2=1., n=24, chi=1.):
    if liquid_density_kg_m3 <= 0 or area_m2 <= 0: raise ValueError('density and area must be positive')
    a=audit(); d=a['reference']; l=a['liquid_feed_mol_s'];g=a['gas_feed_mol_s']
    q_l=d['liquid_mass_flow_kg_s']/liquid_density_kg_m3
    # Assumes all feed gas species are ideally gaseous at 20 bar and 400 K,
    # not a validated vapor-liquid equilibrium flash.
    q_g=sum(g.values())*R*d['inlet_T_K']/(d['pressure_bar']*1e5)
    li=[l[k]/q_l for k in ('octene','H2','octane')]
    ga=[g[k]/q_g for k in ('octene','H2','octane')]
    p=Params(n=n,chi=chi, li_A=li[0],li_H=li[1],ga_A=ga[0],ga_H=ga[1],
        u_L=q_l/area_m2,u_G=q_g/area_m2,T_in=d['inlet_T_K'],
        t_step=1e9,t_end=1100.,wet_init=.88)
    return p,{'liquid_density_assumption_kg_m3':liquid_density_kg_m3,
        'area_assumption_m2':area_m2, 'liquid_volumetric_flow_m3_s':q_l,
        'gas_ideal_volumetric_flow_m3_s':q_g,
        'liquid_superficial_velocity_m_s':p.u_L,'gas_superficial_velocity_m_s':p.u_G,
        'liquid_inlet_reactive_concentrations_mol_m3':dict(zip(('A','H2','P'),li)),
        'gas_inlet_reactive_concentrations_mol_m3':dict(zip(('A','H2','P'),ga)),
        'known_limits':['All gas feed assumed ideal at 20 bar 400 K; no VLE flash',
        'Liquid density assumed, not reported in Table 2',
        'Area 1 m2 treated as cross section of 2012 2D test geometry; actual depth unspecified',
        'Bulk holdups and kinetic/transport/thermal parameters synthetic',
        'Dodecane contributes gas molar volume and liquid mass but not differential energy/species equations']}

def run(out=Path('results/v9_reference_coupled')):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rows=[]; sample=[]
    for rho in [600.,700.,800.]:
      for chi in [1.,10.]:
        p,inp=map_reference_feed(liquid_density_kg_m3=rho,chi=chi)
        ts,Y=simulate(p)
        N,T,s=unpack(Y[:,-1],p); c=N/volumes(s,p)[:,None,:]
        row={'liquid_density_assumed_kg_m3':rho,'chi':chi,
            'final_time_s':float(ts[-1]),'inlet_A_liq_mol_m3':float(p.li_A),
            'outlet_A_liq_mol_m3':float(c[0,0,-1]),
            'liquid_concentration_drop_fraction':float(1-c[0,0,-1]/p.li_A),
            'peak_temperature_K':float(T.max()),'delta_T_K':float(T.max()-p.T_in),
            'min_inventory':float(N.min()),'min_pore_saturation':float(s.min()),
            'max_pore_saturation':float(s.max())}
        rows.append(row)
        if rho==700.:
          z=(np.arange(p.n)+.5)*p.length/p.n
          sample.append((chi,z,T,c[0,0],c[2,0]))
    (out/'mapped_reference_feed.json').write_text(json.dumps(map_reference_feed()[1],indent=2))
    (out/'scenario_results.json').write_text(json.dumps(rows,indent=2))
    with (out/'scenario_results.csv').open('w',newline='') as f:
      w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fig,axs=plt.subplots(2,1,figsize=(8,7),sharex=True)
    for chi,z,T,liq,wet in sample:
      axs[0].plot(z,T,label=f'chi={chi:g}')
      axs[1].plot(z,liq,label=f'Bulk A, chi={chi:g}')
      axs[1].plot(z,wet,'--',label=f'Pore A, chi={chi:g}')
    axs[0].set_ylabel('Temperature (K)');axs[1].set_ylabel('A concentration (mol/m3)')
    axs[1].set_xlabel('Axial distance (m)')
    for ax in axs:ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('v0.9 Table-2-feed coupled synthetic model — NOT 2012 validation')
    fig.tight_layout();fig.savefig(out/'axial_reference_feed.png',dpi=170);plt.close(fig)
    report="""# v0.9: Table 2 feed coupled to 1D solver

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
"""
    (out/'REFERENCE_COUPLING_REPORT.md').write_text(report)
    print(json.dumps(rows,indent=2))
    return rows
if __name__=='__main__':run()
