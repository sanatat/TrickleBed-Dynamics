"""Qualitative axial-profile audit against text-described 2012 Figure 2.
Does NOT imply reproduction of the original CFD calculations.
"""
from dataclasses import replace
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v6_model import Params,simulate,unpack,volumes

OUT=Path('results/v7_validation'); OUT.mkdir(parents=True,exist_ok=True)
profiles={}; summary=[]
for chi in [1.,10.]:
    # Matching published Table 2 geometry, not its unprovided fitted parameters.
    p=Params(n=75, chi=chi, wet_init=.88, t_step=1e6,t_end=1800.)
    tt,YY=simulate(p)
    N,T,s=unpack(YY[:,-1],p)
    C=N/volumes(s,p)[:,None,:]
    z=(np.arange(p.n)+.5)*p.length/p.n
    # 'reacted region' proxy: axial portion accounting for 90% of total consumption
    # inferred from integrated product formation, not original paper's region definition.
    product=C[2,2] # pore-liquid isooctane proxy
    l_oct=C[0,0]; wet_oct=C[2,0]
    profiles[chi]=(z, l_oct,wet_oct,product,T,s)
    x=1.-float(l_oct[-1]/p.li_A)
    summary.append({'chi':chi,'time_s':float(tt[-1]),'model_outlet_liquid_proxy_conversion':x,
        'model_peak_T_K':float(T.max()),'model_max_temp_rise_K':float(T.max()-p.T_in),
        'model_outlet_liquid_A_mol_m3':float(l_oct[-1]),'model_wetting_mean':float(s.mean())})
    with (OUT/f'axial_chi{int(chi)}.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['z_m','liquid_A_mol_m3','wet_A_mol_m3','wet_P_mol_m3','T_K','pore_saturation'])
        w.writerows(zip(z,l_oct,wet_oct,product,T,s))

fig,axs=plt.subplots(2,1,figsize=(8,8),sharex=True)
for chi,(z,a,aw,prod,T,s) in profiles.items():
    axs[0].plot(z,a,label=f'Liquid isooctene proxy, chi={chi:g}')
    axs[0].plot(z,aw,ls='--',label=f'Pore liquid isooctene proxy, chi={chi:g}')
    axs[1].plot(z,T,label=f'chi={chi:g}')
axs[0].set_ylabel('Isooctene proxy (mol/m3)');axs[0].legend(fontsize=8)
axs[1].set_ylabel('Bed temperature (K)');axs[1].set_xlabel('Axial position from inlet (m)');axs[1].legend()
fig.suptitle('v0.7 synthetic axial comparison — NOT a reproduction of 2012 Figure 2')
fig.tight_layout();fig.savefig(OUT/'axial_profiles.png',dpi=160);plt.close(fig)

# Published quantities transcribed from text and Table 2, not digitized from figure.
ref={'reacted_axial_extent_chi1_approx_m':1.0,'reacted_axial_extent_chi10_approx_m':.25,
     'T_in_K':400,'T_near_exit_chi1_approx_K':500,'T_plateau_location_chi10_approx_m':.3,
     'notes':'Article 2012 page 6 (printed page 402); approximate narrative values, NOT digitized curve data.'}
with (OUT/'paper_text_targets.json').open('w') as f:json.dump(ref,f,indent=2)
with (OUT/'model_metrics.json').open('w') as f:json.dump(summary,f,indent=2)
with (OUT/'VALIDATION_REPORT.md').open('w') as f:
    f.write('''# v0.7: Axial validation gap assessment (2012 paper)\n\n''')
    f.write('''**Status: NOT quantitatively validated.** The output is model-vs-literature assessment, not a CFD reproduction.\n\n''')
    f.write('''## Paper-reported targets (2012 paper, printed p.402)\n\n- In Test 1, reaction occurs over ~1 m at χ=1 and ~0.25 m at χ=10.\n- Temperature rises from 400 K to ~500 K at χ=1; with χ=10 the similar maximum is reached around z=0.3 m.\n- Table 2 geometry: 1.5 m height, 3 mm pellet diameter, bed porosity 0.37, pellet porosity 0.5, 20 bar and 400 K.\n\n''')
    f.write('''## What this run actually calculates\n\n- Synthetic four-compartment (mobile gas, mobile liquid, pore gas and pore liquid) reduced-order model.\n- 75 axial cells, constant fitted transport proxies, one lumped olefin reaction and a prescribed wetting trajectory.\n- Not a temperature-dependent vapour-liquid equilibrium calculation or a resolved hydrodynamic model.\n- Uses fixed illustrative inlet concentrations and kinetic parameters, NOT the published mass-fraction input or the fitted Lylykangas constants.\n- The 2012 figure uses molar fractions; current plots use mol/m³; direct curve-shape error metrics would be misleading.\n\n''')
    f.write('''## Simulation outputs\n\n| χ | Outlet proxy conversion | Maximum T (K) | ΔT max (K) |\n|---|---:|---:|---:|\n''')
    for s in summary:
        f.write(f"| {s['chi']:g} | {100*s['model_outlet_liquid_proxy_conversion']:.2f}% | {s['model_peak_T_K']:.2f} | {s['model_max_temp_rise_K']:.2f} |\n")
    f.write('''\n## Finding\n\nThe heat rise is substantially smaller than the ~100 K reported in the paper, and the code does not model all physics needed to measure corresponding axial reaction lengths. **These are material discrepancies, not a validated match.** Matching peak temperature by tuning only a synthetic heat coefficient would be curve-fitting rather than independent validation.\n\n## Missing inputs and next scientific work\n\n1. Numerical kinetic constants for the two isooctene isomers, adsorption constants and temperature dependence in Lylykangas et al. (2003).\n2. Hydrodynamic/wetting relation of Lappalainen et al. (2008/2009) including liquid holdup and pressure drop.\n3. Phase equilibrium and enthalpy properties for the four-component gas-liquid mixture, evaporation and condensation.\n4. Actual 2012 inlet composition conversions and consistent mass, energy and species equations.\n5. Reference data digitized from Figure 2 with a documented procedure, followed by grid/time-step sensitivity and numerical error metrics.\n\n**Research claim appropriate for applications:** Implemented and tested a paper-inspired prototype and quantified major gaps with its reference. Do not claim reproduction of Figure 2.\n''')
print(json.dumps(summary,indent=2))
