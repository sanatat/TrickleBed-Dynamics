"""Literature-supported activation-energy sensitivity, not calibrated kinetics.

Lylykangas et al. (2003) report Ea(TMP-1)=34 and Ea(TMP-2)=49 kJ/mol.
The existing v0.9 lumped A state cannot represent two isomers separately.
Therefore runs are alternate Ea values applied to a *lumped proxy*, not
physical two-isomer kinetics. Prefactors and adsorption constants stay synthetic.
"""
from pathlib import Path
from dataclasses import replace
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reference_coupled_v9 import map_reference_feed
from v6_model import simulate,unpack,volumes
R=8.314462618
EA={"TMP-1":34000.0,"TMP-2":49000.0} # J/mol, verified Lylykangas abstract

def arrhenius_ratio(Ea,T,Tref=400.):
    if np.any(np.asarray(T)<=0) or Tref<=0 or Ea<0:raise ValueError('invalid temperature/Ea')
    return np.exp(-Ea/R*(1/np.asarray(T)-1/Tref))

def evaluate(out=Path('results/v10_kinetics_provenance')):
    out=Path(out);out.mkdir(exist_ok=True,parents=True)
    data=[]
    for isomer,Ea in EA.items():
        for chi in (1.,10.):
            p,meta=map_reference_feed(liquid_density_kg_m3=700.,chi=chi,n=24)
            p=replace(p,Ea=Ea)
            t,Y=simulate(p)
            N,T,s=unpack(Y[:,-1],p);C=N/volumes(s,p)[:,None,:]
            data.append({'activation_energy_proxy':isomer,'Ea_J_mol':Ea,'chi':chi,
              'peak_temperature_K':float(T.max()),'outlet_bulk_A_mol_m3':float(C[0,0,-1]),
              'outlet_bulk_A_drop_percent':float(100*(1-C[0,0,-1]/p.li_A)),
              'kinetics_status':'SYNTHETIC k_ref and adsorption; literature Ea only'})
    with (out/'results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    temps=np.linspace(360,480,150)
    fig,ax=plt.subplots(figsize=(8,5))
    for iso,Ea in EA.items():ax.plot(temps,arrhenius_ratio(Ea,temps),label=f'{iso}, Ea={Ea/1000:.0f} kJ/mol')
    ax.axvline(400,color='gray',ls=':',label='2012 inlet temperature 400 K')
    ax.set_xlabel('Temperature (K)');ax.set_ylabel('Relative k(T) / k(400 K)')
    ax.set_title('Verified activation energies; unknown absolute rate constants')
    ax.grid(alpha=.25);ax.legend();fig.tight_layout();fig.savefig(out/'arrhenius_relative.png',dpi=170);plt.close(fig)
    (out/'source_provenance.json').write_text(json.dumps({
       'primary_kinetics':'Lylykangas, Rautanen, Krause (2003), AIChE J. 49:1508-1515, DOI 10.1002/aic.690490614',
       'verified_from':'Publisher abstract', 'Ea_TMP1_J_mol':34000,'Ea_TMP2_J_mol':49000,
       'unknown':['k_ref for both isomers','adsorption constants and their units','isomer-resolved feed fractions in 2012 CFD','complete kinetics table'],
       'warning':'The v0.9 model has only one octene species. These are ALTERNATIVE sensitivity cases; not two-isomer simulations.'},indent=2))
    (out/'KINETICS_REPORT.md').write_text("""# v1.0 Literature-supported temperature dependence

**NOT a reproduction of Gorshkova et al. (2012).**

- Lylykangas et al. (2003), AIChE Journal 49:1508-1515: publisher abstract explicitly reports activation energies of **34 kJ/mol for TMP-1** and **49 kJ/mol for TMP-2**, on Ni/Al2O3.
- Computed Arrhenius *relative-rate* curves are physically meaningful conditional on Arrhenius temperature dependence, but they do **not** determine the absolute reaction rates.
- Legacy v0.9 model represents isooctenes as one lumped concentration, so the two Ea scenarios are **proxies** applied separately to the same lumped system, NOT a chemical simulation of both isomers.
- The adsorption coefficients, pre-exponential factors, mass transfer, heat transfer and gas-liquid equilibrium closures remain **uncalibrated illustrations**.
- The model may be especially misleading about absolute thermal rise and reaction conversion. This release does not claim a quantitative fit to Fig. 2 (2012).

## Next scientifically necessary work

1. Obtain the **complete original Lylykangas 2003 paper with parameter table** and confirm units, rate law, and both-isomer constants.
2. Implement separate TMP-1, TMP-2, isooctane states and any measured isomerization pathways.
3. Recover the full wetting closure of Lappalainen (2008) and physical thermal/phase-equilibrium properties.
4. Compare digitized 2012 Fig. 2 values only after closures are verified.
""")
    print(json.dumps(data,indent=2))
    return data
if __name__=='__main__':evaluate()
