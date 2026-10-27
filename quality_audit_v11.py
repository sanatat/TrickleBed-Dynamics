"""Numerical integrity and parameter-uncertainty audit for the uncalibrated v6 model.
No run in this file is experimental validation or reproduction of published CFD.
"""
from dataclasses import replace
from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v6_model import Params,initial,rhs,simulate,unpack,volumes


def residual_inventory(y,p,reaction=False,transfer=True):
    # For no-flow closed bed: A+P and H+P must be conserved by sources.
    dy=rhs(17.,y,p,flow=False,reaction=reaction,transfer=transfer)
    dn,_,_=unpack(dy,p)
    return np.array([np.max(np.abs(dn[:,0,:].sum(axis=0)+dn[:,2,:].sum(axis=0))),
                     np.max(np.abs(dn[:,1,:].sum(axis=0)+dn[:,2,:].sum(axis=0)))])


def metrics(p):
    t,Y=simulate(p)
    N,T,S=unpack(Y[:,-1],p)
    C=N/volumes(S,p)[:,None,:]
    return {'peak_temperature_K':float(T.max()), 'outlet_alkene_mol_m3':float(C[0,0,-1]),
      'outlet_concentration_drop_fraction':float(1-C[0,0,-1]/p.li_A),
      'min_species_inventory_mol_m3bed':float(N.min()),'saturation_min':float(S.min()),
      'saturation_max':float(S.max())}


def run(output=Path('results/v11_credibility')):
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    p=Params(n=12,t_end=210.,t_step=100.,wet_init=.2)
    y=initial(p)
    # perturb compartments at deterministic positive values; avoid artificial cancellation from equilibrium initial conditions
    rng=np.random.default_rng(17); yt=y.copy().reshape(14,p.n)
    yt[:12]*=rng.uniform(.65,1.35,size=(12,p.n))
    residuals={'closed_no_reaction':residual_inventory(yt.ravel(),p,False).tolist(),
        'closed_with_reaction':residual_inventory(yt.ravel(),p,True).tolist()}
    scenarios=[]
    for name,changes in [('baseline',{}),('kinetics_x0.5',{'k_ref':p.k_ref*.5}),
                         ('kinetics_x2',{'k_ref':p.k_ref*2}),('pore_transfer_x0.5',{'k_LW':p.k_LW*.5}),
                         ('pore_transfer_x2',{'k_LW':p.k_LW*2}),
                         ('initially_prewet',{'wet_init':.88}),
                         ('effective_diffusion_x10',{'chi':10.})]:
        row={'scenario':name, **metrics(replace(p,**changes))};scenarios.append(row)
    grid=[]
    for n in (8,12,18):grid.append({'n_cells':n,**metrics(replace(p,n=n))})
    for filename,rows in [('sensitivity.csv',scenarios),('grid_diagnostics.csv',grid)]:
        with (out/filename).open('w',newline='') as f:
            w=csv.DictWriter(f,rows[0].keys());w.writeheader();w.writerows(rows)
    base=scenarios[0]; grid_difference=abs(grid[-1]['outlet_alkene_mol_m3']-grid[-2]['outlet_alkene_mol_m3'])
    audit={'model_status':'uncalibrated 1D synthetic model, NOT 2012 CFD reproduction',
      'max_closed_inventory_source_residual':max(max(row) for row in residuals.values()),
      'residuals_mol_m3_s':residuals,'grid_12_to_18_absolute_outlet_difference_mol_m3':grid_difference,
      'grid_12_to_18_relative_outlet_difference':grid_difference/max(abs(grid[-1]['outlet_alkene_mol_m3']),1e-12),
      'baseline':base,'provenance':'Gorshkova et al. 2012 supplies geometry and phenomenological framework; kinetic/transfer/wetting constants uncalibrated'}
    (out/'audit.json').write_text(json.dumps(audit,indent=2))
    fig,ax=plt.subplots(figsize=(9,4.7));labels=[x['scenario'] for x in scenarios];yy=[x['peak_temperature_K']-p.T_in for x in scenarios]
    ax.barh(labels,yy);ax.invert_yaxis();ax.set_xlabel('Peak temperature rise (K)')
    ax.set_title('Parameter sensitivity — hypothetical coefficients, not experimental evidence')
    fig.tight_layout();fig.savefig(out/'sensitivity_peak_temperature.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4.5))
    ax.plot([v['n_cells'] for v in grid],[v['outlet_alkene_mol_m3'] for v in grid],marker='o')
    ax.set_xlabel('Axial cells');ax.set_ylabel('Outlet alkene concentration (mol/m³)')
    ax.set_title('Mesh diagnostic (not a convergence proof)');fig.tight_layout()
    fig.savefig(out/'mesh_diagnostic.png',dpi=170);plt.close(fig)
    report=f"""# v1.1 scientific credibility audit

**Status: numerical demonstrator, not validated against experiment or the 2012 CFD model.**

## What was measured in code

- In a closed system the discrete source terms conserve alkene-equivalents A+P and hydrogen-equivalents H+P. Maximum differential residual across all states tested: **{audit['max_closed_inventory_source_residual']:.3e} mol m⁻³ s⁻¹**. This is a code-internal check only.
- In a short transient run with prescribed flow perturbation, changing the *assumed* kinetic scale, pore-transfer coefficient, starting wetness and effective diffusion changes the predicted peak temperature. See `sensitivity.csv`. This is an identifiability warning, not fitted parameter estimates.
- Mesh diagnostic: 12 versus 18 axial cells give an outlet-alkene difference of **{grid_difference:.4g} mol m⁻³** ({100*audit['grid_12_to_18_relative_outlet_difference']:.2f}%). A three-grid comparison by itself is not a verification of grid independence.

## Research-relevant unresolved gaps

1. Full source kinetic constants from Lylykangas et al., including two distinct octene isomers and their isomerization. The existing single-lumped-species kinetic coefficient is not published data.
2. Independently validated wetting hysteresis / capillarity closure and pressure-drop coupling.
3. Vapor-liquid flash and latent heat; temperature-dependent thermophysical properties.
4. Complete species inventory including inert dodecane, accurate gas/liquid volume balances, and pressure equations.
5. Reproducible CFD/profile reference for quantitative comparison, not just visually read paper figures.
6. Stiffness/parameter-identifiability studies with experimental data before making thermal-safety claims.

## Interpretation

A large increase in model test count does not constitute scientific validation. The 2012 work uses transient Eulerian CFD; this project remains a 1D compartment surrogate. Figure 2 of the publication has **not** been quantitatively reproduced. The model uses arbitrary k and heat parameters and must not be used for process design, safety decisions, or industrial predictions.
"""
    (out/'CREDIBILITY_REPORT.md').write_text(report)
    return audit
if __name__=='__main__':print(json.dumps(run(),indent=2))
