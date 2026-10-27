"""Evidence-weighted diagnosis of v2.7 heterogeneous RTD over-retention.

Published experiment: Sales et al. (2005) Figure 5, vector-digitized in v2.6.
This study does not claim to recover original, unobserved experimental parameters.
"""
from __future__ import annotations
from dataclasses import replace, asdict
import csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from heterogeneous_rtd_v27 import (Params, build_linear_system, inventory_weights,
                                   check_global_inventory, tracer_response, load_vector_curve,
                                   metrics, effective_dispersion_per_bed)

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'/'v28_root_cause'

def crossing_time(t,c,level=.5):
    """Time when decreasing negative-step trace first crosses level."""
    t=np.asarray(t,float);c=np.asarray(c,float)
    hit=np.flatnonzero(c<=level)
    if len(hit)==0:return float('nan')
    i=int(hit[0]);
    if i==0:return float(t[0])
    if c[i-1]==c[i]:return float(t[i])
    return float(np.interp(level,[c[i],c[i-1]],[t[i],t[i-1]]))

def normal_step_moment(t,y):
    """Mean delay if ideal step from 1 to 0 starts at t=0; diagnostic only.
    Here real inlet is NOT ideal step, so compare two traces' integrals instead.
    """
    return float(np.trapezoid(np.clip(y,0,1),t))

def evaluate():
    OUT.mkdir(parents=True,exist_ok=True)
    t,ci=load_vector_curve('fig5_experimental_entrance')
    to,obs=load_vector_curve('fig5_experimental_exit')
    if not np.allclose(t,to):raise ValueError('Inlet/outlet grids differ')
    base=Params()
    configs={
        'axial_only_reference':(base,False),
        'full_heterogeneous_reference':(base,True),
        'no_adsorption':(replace(base,K_ads_m3_kg=0.),True),
        'low_pore_access_0p25':(replace(base,accessible_pore_fraction=.25),True),
        'high_intraparticle_diffusion_table2':(replace(base,effective_pore_D_m2_s=1.8e-7),True),
        'low_surface_exchange_0p1':(replace(base,interface_scale=.1),True),
        'fast_adsorption_relaxation':(replace(base,adsorption_relaxation_s_inv=.3),True),
    }
    rows=[]; predictions={}
    observed_delay=crossing_time(t,obs)-crossing_time(t,ci)
    for name,(p,pore) in configs.items():
        yy=tracer_response(t,t,ci,p,enable_pore=pore)
        predictions[name]=yy
        rows.append({'variant':name,'rmse':metrics(yy,obs),
                     'observed_t50_delta_s':observed_delay,
                     'model_t50_delta_s':crossing_time(t,yy)-crossing_time(t,ci),
                     'rmse_transition_300_900':metrics(yy[(t>=300)&(t<=900)],obs[(t>=300)&(t<=900)]),
                     'param_source':'fixed publisher + explicitly declared assumptions; no fit'})
    # Profile accessibility as an effective latent volume fraction; a constrained *single*
    # parameter fit on Fig 5. It is a diagnostic correction, not physical identification.
    fit_grid=[.05,.10,.20,.35,.5,.75,1.0]
    prof=[]
    grid_predictions={}
    train=(t>=300)&(t<=720)
    holdout=(t>720)&(t<=1050)
    for x in fit_grid:
        p=replace(base,accessible_pore_fraction=x)
        yy=tracer_response(t,t,ci,p)
        grid_predictions[x]=yy
        prof.append({'accessible_pore_fraction':x,'rmse':metrics(yy,obs),
                     'train_300_720_rmse':metrics(yy[train],obs[train]),
                     'heldout_750_1050_rmse':metrics(yy[holdout],obs[holdout]),
                     'mean_delay_s':crossing_time(t,yy)-crossing_time(t,ci)})
    best=min(prof,key=lambda r:r['rmse'])
    best_train=min(prof,key=lambda r:r['train_300_720_rmse'])
    predictions['best_accessibility_GRID_FIT']=grid_predictions[best['accessible_pore_fraction']]
    predictions['best_accessibility_TRANSITION_TRAIN']=grid_predictions[best_train['accessible_pore_fraction']]
    # Joint-parameter grid is a diagnostic for trade-offs, not a calibrated 2D fit.
    # Evaluated on the same Figure 5 (in-sample), not independent validation.
    joint=[]
    for access in (.05,.1,.2,.4,.8):
        for dscale in (.5,1.,2.):
            if dscale==1. and access in grid_predictions:
                y=grid_predictions[access]
            else:
                y=tracer_response(t,t,ci,replace(base,accessible_pore_fraction=access,
                                                 axial_dispersion_scale=dscale))
            joint.append({'accessible_pore_fraction':access,'dispersion_multiplier':dscale,
                          'rmse':metrics(y,obs)})
    joint_best=min(joint,key=lambda q:q['rmse'])
    with (OUT/'parameter_tradeoff.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(joint[0].keys()));wr.writeheader();wr.writerows(joint)
    # Reaction of the model to 3 geometrical assumptions, no extra fit.
    pore_capacity= (1-base.bed_porosity_assumed)*base.wetting_efficiency*(base.internal_porosity+base.particle_bulk_density*base.K_ads_m3_kg)
    mobile_capacity=base.liquid_holdup
    ratio=pore_capacity/mobile_capacity
    # Analytical equilibrium residence-time estimate is illustrative, not an RTD fit.
    advective_lag=base.length_m*base.liquid_holdup/base.superficial_u_m_s
    equilibrium_lag=base.length_m*(base.liquid_holdup+pore_capacity)/base.superficial_u_m_s
    a,b=build_linear_system(base)
    stationary_res, inlet_flux_err=check_global_inventory(base)
    # report why a fitted fraction should not be over interpreted
    with (OUT/'ablation.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=rows[0].keys());wr.writeheader();wr.writerows(rows)
    with (OUT/'pore_accessibility_profile.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(prof[0].keys()));wr.writeheader();wr.writerows(prof)
    with (OUT/'figure5_diagnostic_predictions.csv').open('w',newline='') as f:
        wr=csv.writer(f);wr.writerow(['time_s','inlet_graphical','outlet_graphical',*predictions])
        wr.writerows(zip(t,ci,obs,*predictions.values()))
    m={
        'source':'Sales 2005 Fig.5 vector-extracted graphical inlet/output (v2.6)',
        'N_graphical_resamples':int(len(t)), 'paper_fig5_flow_L_h':5,
        'observed_outlet_minus_inlet_t50_s':observed_delay,
        'bed_mobile_capacity_per_m3':mobile_capacity,
        'modeled_wetted_intraparticle_plus_adsorption_capacity_per_m3':pore_capacity,
        'assumed_pore_to_mobile_capacity_ratio':ratio,
        'approx_advection_time_mobile_s':advective_lag,
        'approx_time_at_full_pore_equilibrium_s':equilibrium_lag,
        'fixed_parameter_ablation':rows,
        'pore_accessibility_profile':prof,
        'best_profile_point_fitted_to_fig5':best,
        'joint_parameter_tradeoff_best_in_sample':joint_best,
        'joint_parameter_tradeoff_all':joint,
        'best_accessibility_fitted_transition_training_only':best_train,
        'transition_train_definition':'graphical 300-720 s',
        'transition_holdout_definition':'graphical 750-1050 s; same experiment, NOT external validation',
        'inventory_identity_max_residual_m_s':stationary_res,
        'inlet_identity_error_m_s':inlet_flux_err,
        'conclusion':'No unique root cause. Extra pore storage/adsorption in the chosen closures shifts response; plausible accessible volume and adsorption kinetics are unmeasured.',
        'critical_qualification':'Best accessible pore fraction is calibrated on same Fig.5 and cannot be interpreted as observed porosity/wetting or as independent validation.',
    }
    (OUT/'metrics.json').write_text(json.dumps(m,indent=2))
    fig,ax=plt.subplots(figsize=(9,5))
    ax.plot(t,ci,':',c='gray',label='Inlet published trace')
    ax.plot(t,obs,'-',lw=2.0,c='black',alpha=.8,label='Outlet graphical trace (resampled)')
    for name,sty in [('axial_only_reference','--'),('full_heterogeneous_reference','-'),('no_adsorption','-.'),('best_accessibility_GRID_FIT','-')]:
        ax.plot(t,predictions[name],sty,lw=1.8,label=name.replace('_',' '))
    ax.set(xlabel='Time (s)',ylabel='Normalized tracer signal',xlim=(0,1800),ylim=(-.03,1.05),title='Sales (2005) Fig. 5: diagnosis of excess model delay')
    ax.grid(alpha=.2);ax.legend(fontsize=7,loc='upper right');fig.tight_layout();fig.savefig(OUT/'ablation_comparison.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.5))
    ax.plot([i['accessible_pore_fraction'] for i in prof],[i['rmse'] for i in prof],'-o')
    ax.axhline(rows[0]['rmse'],linestyle='--',c='gray',label='No-pore reference RMSE')
    ax.set(xlabel='Assumed accessible pore fraction (NOT measured)',ylabel='RMSE (normalized signal)',title='One-parameter diagnostic profile on same Fig. 5')
    ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(OUT/'accessibility_profile.png',dpi=180);plt.close(fig)
    xx=sorted(set(i['accessible_pore_fraction'] for i in joint))
    yy=sorted(set(i['dispersion_multiplier'] for i in joint))
    zz=np.array([[next(i['rmse'] for i in joint if i['accessible_pore_fraction']==x and i['dispersion_multiplier']==d) for x in xx] for d in yy])
    fig,ax=plt.subplots(figsize=(8,4.4))
    im=ax.imshow(zz,origin='lower',aspect='auto',cmap='viridis',extent=[-.5,len(xx)-.5,-.5,len(yy)-.5])
    ax.set(xticks=np.arange(len(xx)),yticks=np.arange(len(yy)),xticklabels=[str(i) for i in xx],yticklabels=[str(i) for i in yy],xlabel='Assumed accessible pore fraction',ylabel='Dax multiplier',title='In-sample RMSE: accessibility–dispersion trade-off')
    fig.colorbar(im,ax=ax,label='RMSE')
    fig.tight_layout();fig.savefig(OUT/'parameter_tradeoff.png',dpi=180);plt.close(fig)
    report=f'''# v2.8 — Physical discrepancy diagnosis, Sales (2005) Figure 5

## What is actually observed

The uploaded Sales et al. (2005) PDF contains experiments. Figure 5 (printed p.450) gives inlet and outlet negative-step tracer traces at 5 L/h liquid. v2.6 extracted **{len(t)} graph-resampled points** for each curve from PDF vector paths. These are **not original apparatus time series** and their neighboring points are not independent measurements.

## Root-cause study, not retrofit confirmation

We held literature-linked values (hL=0.17, f_e=0.38, Pe=4.84, D_eff=7.65e-10 m²/s from paper prose, k_LS chosen within the reported range) constant unless the ablation specifies otherwise. The published article does not supply the pore-accessibility fraction, the adsorption relaxation constant, and the packed-bed porosity used by our numerical model. The model is motivated by **Eqs. 7–9**, not an exact implementation of **Eqs. 13–15**.

Fixed-choice ablation results (do not treat these as equally calibrated models):

| Scenario | Graphical RMSE | Model input-output 50% lag (s) |
|---|---:|---:|
'''
    for r in rows:
        report+=f"| {r['variant']} | {r['rmse']:.5f} | {r['model_t50_delta_s']:.1f} |\n"
    report+=f'''
Observed graph-derived inlet/output 50% lag: **{observed_delay:.1f} s**.

## Intraparticle storage accounting (specific to model assumptions)

Chosen mobile capacity hL = **{mobile_capacity:.4f} m³/m³ bed**. With assumed bed porosity {base.bed_porosity_assumed}, f_e {base.wetting_efficiency}, particle porosity {base.internal_porosity}, adsorption density {base.particle_bulk_density} kg/m³ and K_ads {base.K_ads_m3_kg} m³/kg, **assumed** fully-accessible inner storage capacity = **{pore_capacity:.4f} m³ equivalent/m³ bed**, i.e., **{ratio:.1%}** of mobile storage. In a hypothetical instantaneous-equilibrium limit these capacities would yield an advective lag of approximately {equilibrium_lag:.0f} s rather than {advective_lag:.0f} s for mobile liquid alone. These analytical times are *not fitted experimental RTD means*. Their difference pinpoints a plausible excess-storage mechanism in our assumed closures, not a proven experimental root cause.

## Diagnostic single-parameter fitting (NOT validation)

Varying an **assumed** additional pore accessibility factor (between 0 and 1) changes the amount of pore storage and transfer area. The best member of our coarse tested grid is **{best['accessible_pore_fraction']:.2f}**, RMSE **{best['rmse']:.5f}**. This is an *in-sample* diagnostic parameter, not physically identified wetting or porosity. Only 7 trial values were evaluated, not a general optimum or an uncertainty interval.

## Temporal holdout diagnostic (still same experimental trace)

Coarse-grid accessibility chosen using only the transition **300–720 s**: **{best_train['accessible_pore_fraction']:.2f}**, training RMSE **{best_train['train_300_720_rmse']:.5f}**; RMSE on subsequent **750–1050 s** points **{best_train['heldout_750_1050_rmse']:.5f}**. This is a *within-curve extrapolation* diagnostic, not an independent run nor instrument-level external validation. The Figure-5 input trace remains known over the whole interval. Correlation between graph-resampled points reduces statistical independence.

## Accessibility–dispersion trade-off

A coarse 5 × 3 two-parameter grid on the **same graph-derived Figure 5** found a lowest sampled RMSE of **{joint_best['rmse']:.5f}** at assumed accessibility **{joint_best['accessible_pore_fraction']:.2f}** and dispersion multiplier **{joint_best['dispersion_multiplier']:.1f}**. This is not physical parameter identification. The grid is designed to expose how model mismatch may be redistributed among uncertain processes; it is not a global optimum or experimental uncertainty distribution.

## Numerics

Discrete global tracer conservation identity residual **{stationary_res:.3g} m/s**, input identity residual **{inlet_flux_err:.3g} m/s**. Conserving the implemented equations does NOT validate closures.

## Limitation and next discriminating experiment

Figure 5 alone does not separately identify true pore accessibility, adsorption capacity, surface exchange and axial dispersion. Published Table 2 and page-450 prose conflict on effective pore diffusivity; resolve with original authors before claiming a unique value. A separate non-adsorbing tracer (KCl is described by the paper) and the adsorption-capable sodium gluconate at the **same flow**, with recorded inlet/outlet, would separate adsorption from hydrodynamic memory better than further tuning the same curve. Figure 6 has no independently measured inlet series in the PDF and cannot validate a Figure-5-derived transfer model without assumptions.

**This is an evidence-based model-discrepancy study, not experimental confirmation of the full reactive Alopaeus CFD model.**
'''
    (OUT/'REPORT.md').write_text(report)
    print(json.dumps({'ablation':{r['variant']:round(r['rmse'],6) for r in rows},'best_grid':best,'capacity_ratio':ratio,'inlet_outlet_lag_s':observed_delay},indent=2))
    return m

if __name__=='__main__':evaluate()
