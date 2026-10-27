"""v2.9: single real graphical tracer vs synthetic two-tracer identifiability audit.

Published paper: Sales et al. (2005), DOI:10.1590/S0104-66322005000300012
Figure 5 has a published sodium-gluconate inlet and outlet; KCl RTD points are
NOT available in this source. Never mix the simulated KCl with real data fitting.

The finite-volume model is based on Eqs. 7-9, NOT the article's exact transfer
function (13-15). Graph vector paths (v2.6) are correlated graphical samples.
"""
from __future__ import annotations
from dataclasses import replace
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from heterogeneous_rtd_v27 import Params, tracer_response, load_vector_curve

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results' / 'v29_dual_tracer'
LOG_A_BOUNDS = (np.log([0.025, 1e-6]), np.log([1.0, 5e-4]))


def parameterized(base: Params, logpar, tracer: str='NaGl') -> Params:
    access, ads = np.exp(logpar)
    if tracer == 'NaGl':
        return replace(base, accessible_pore_fraction=float(access), K_ads_m3_kg=float(ads))
    if tracer == 'KCl_assumed_nonadsorbing':
        # KCl adsorption set to ZERO *by hypothetical experimental design*,
        # and diffusivity is chosen for illustration; neither was measured here.
        return replace(base, accessible_pore_fraction=float(access), K_ads_m3_kg=0.,
                       effective_pore_D_m2_s=1.3e-9)
    raise ValueError('Unknown tracer')


def predicted(t, ci, p, logpar, tracer='NaGl'):
    return tracer_response(t, t, ci, parameterized(p, logpar, tracer))


def sensitivity(t, ci, base, logpar, tracer, h=.065):
    cols=[]
    for k in range(2):
        step=np.zeros(2); step[k]=h
        plus=predicted(t,ci,base,logpar+step,tracer)
        minus=predicted(t,ci,base,logpar-step,tracer)
        cols.append((plus-minus)/(2*h))
    return np.column_stack(cols)


def fisher(J, assumed_sigma=.015):
    return (J.T @ J)/(assumed_sigma**2)


def fiducial_diagnostics(F):
    eigen=np.linalg.eigvalsh(F)
    cond=float(eigen[-1]/eigen[0]) if eigen[0]>1e-14 else float('inf')
    try:
        inv=np.linalg.inv(F)
        corr=float(inv[0,1]/np.sqrt(inv[0,0]*inv[1,1]))
        std=np.sqrt(np.diag(inv)).tolist()
    except np.linalg.LinAlgError:
        corr=float('nan'); std=[float('nan')]*2
    return {'smallest_eigenvalue':float(eigen[0]),'largest_eigenvalue':float(eigen[-1]),
            'condition_number':cond,'log_parameter_correlation':corr,
            'optimistic_crlb_logparam_sd':std,
            'logdet':float(np.sum(np.log(eigen))) if np.all(eigen>0) else float('-inf')}


def fit_experimental(t, cin, cout, base):
    """Exploratory in-sample fit to ONE real sodium-gluconate curve only."""
    bounds=LOG_A_BOUNDS
    subset=(t>=150)&(t<=1260)
    def residual(x):
        return (predicted(t,cin,base,x,'NaGl')-cout)[subset]
    fit=least_squares(residual,np.log([.16,8.2e-5]),bounds=bounds, max_nfev=75,
                      diff_step=2e-2,xtol=3e-4,ftol=3e-4,gtol=3e-4)
    y=predicted(t,cin,base,fit.x,'NaGl')
    return {'result':fit, 'prediction':y,
            'rmse_full':float(np.sqrt(np.mean((y-cout)**2))),
            'rmse_transit':float(np.sqrt(np.mean((y[subset]-cout[subset])**2)))}


def simulated_recovery(t, cin, base, seed=2026, sigma=.01):
    """Synthetic-only controlled experiment. No published KCl observations used."""
    rng=np.random.default_rng(seed)
    truth=np.log([.2,8.2e-5])
    time_mask=(t>=120)&(t<=1260)&(np.arange(t.size)%2==0)
    yn=predicted(t,cin,base,truth,'NaGl')
    yk=predicted(t,cin,base,truth,'KCl_assumed_nonadsorbing')
    noisy_na=yn+ rng.normal(0,sigma,size=t.size)
    noisy_kcl=yk+rng.normal(0,sigma,size=t.size)
    def fit(two):
        def residual(x):
            res=[(predicted(t,cin,base,x,'NaGl')-noisy_na)[time_mask]/sigma]
            if two: res.append((predicted(t,cin,base,x,'KCl_assumed_nonadsorbing')-noisy_kcl)[time_mask]/sigma)
            return np.concatenate(res)
        return least_squares(residual, np.log([.3, 5e-5]), bounds=LOG_A_BOUNDS,
                             max_nfev=60, diff_step=2e-2, ftol=5e-5, xtol=5e-5, gtol=5e-5)
    fa,fb=fit(False),fit(True)
    return {'true':truth, 'na_true':yn, 'kcl_true':yk, 'na_noisy':noisy_na,
            'kcl_noisy':noisy_kcl, 'na_fit':fa,'dual_fit':fb,'train_mask':time_mask,
            'sigma':sigma,'seed':seed}


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    t,ci=load_vector_curve('fig5_experimental_entrance')
    tout,obs=load_vector_curve('fig5_experimental_exit')
    if not np.array_equal(t,tout): raise RuntimeError('Misaligned vectors')
    p=Params()
    fit=fit_experimental(t,ci,obs,p)
    experiment=simulated_recovery(t,ci,p)
    # Sensitivity/Jacobian assessment around a DECLARED hypothetical truth.
    truth=experiment['true']
    # More sparsely sampled than graphical 59-point curves; still hypothetical
    # independent Gaussian noise is an assumption and CRLB is optimistic.
    sel=(t>=120)&(t<=1260)&(np.arange(t.size)%2==0)
    Jn=sensitivity(t,ci,p,truth,'NaGl')[sel]
    Jk=sensitivity(t,ci,p,truth,'KCl_assumed_nonadsorbing')[sel]
    Fn=fisher(Jn); Fk=fisher(Jk)
    designs={'NaGl_only':fiducial_diagnostics(Fn),
             'NaGl_repeated_twice':fiducial_diagnostics(2*Fn),
             'NaGl_plus_assumed_KCl':fiducial_diagnostics(Fn+Fk)}
    # Sensitivity to the unknown baseline accessibility: these are LOCAL,
    # hypothetical Fisher calculations, not confidence intervals for Figure 5.
    robustness=[]
    for trial_access in (.10,.20,.40):
        for trial_ads in (4.1e-5,8.2e-5,1.64e-4):
            trial=np.log([trial_access,trial_ads])
            Jna=sensitivity(t,ci,p,trial,'NaGl')[sel]
            Jaux=sensitivity(t,ci,p,trial,'KCl_assumed_nonadsorbing')[sel]
            dia=fiducial_diagnostics(fisher(Jna))
            dib=fiducial_diagnostics(fisher(Jna)+fisher(Jaux))
            robustness.append({'assumed_pore_accessibility':trial_access,
                               'assumed_NaGl_Kads_m3_kg':trial_ads,
                               'single_min_eigen':dia['smallest_eigenvalue'],
                               'dual_min_eigen':dib['smallest_eigenvalue'],
                               'single_condition':dia['condition_number'],
                               'dual_condition':dib['condition_number'],
                               'min_eigen_gain':dib['smallest_eigenvalue']/dia['smallest_eigenvalue']})
    # Physical storage identifiability at first moment: Eq 16 in source article
    # uses (1-e) f_e (epsilon_in + rho_particle*K_ads). KCl adsorption here ZERO.
    porosity=p.bed_porosity_assumed
    volume=np.pi*(.04/2)**2*.80  # reactor as published: diameter 0.04 m height 0.8 m
    q_l=5/1000/3600
    hL=p.liquid_holdup
    fe=p.wetting_efficiency
    mu_nagl=volume/q_l*(hL + (1-porosity)*fe*(p.internal_porosity+p.particle_bulk_density*p.K_ads_m3_kg))
    mu_kcl=volume/q_l*(hL + (1-porosity)*fe*p.internal_porosity)
    delta=mu_nagl-mu_kcl
    # Important: units/check for equation storage differ in paper's nomenclature;
    # this is a *conditional demonstration*, not reported experimental first moments.
    data={
      'experiment_source':'Sales et al. 2005 Figure 5; Q_L 5 L/h and Q_G 300 L/h; user-provided download.pdf',
      'experiment_data_type':'Published graphical time series, not independent raw instrument readings',
      'kcl_data':'NOT AVAILABLE IN PROVIDED ARTICLE; ALL KCl trajectories SYNTHETIC',
      'rtd_model':'1D axial dispersion + spherical intraparticle diffusion + reversible adsorbate state, v2.8 family',
      'experiment_NaGl_only':{'fit_accessible_pore_fraction':float(np.exp(fit['result'].x[0])),
                              'fit_Kads_m3_kg':float(np.exp(fit['result'].x[1])),
                              'fit_rmse_full':fit['rmse_full'],
                              'fit_rmse_transition':fit['rmse_transit'],
                              'on_parameter_bounds': bool(np.any(np.isclose(fit['result'].x,LOG_A_BOUNDS[0],atol=.035)) or np.any(np.isclose(fit['result'].x,LOG_A_BOUNDS[1],atol=.035))),
                              'optimization_success':bool(fit['result'].success),
                              'identification_warning':'Fit is in-sample and joint parameters may be nonunique. No external experimental validation.'},
      'controlled_synthetic_design':{
       'true_pore_accessibility':float(np.exp(truth[0])),
       'true_NaGl_Kads_m3_kg':float(np.exp(truth[1])),
       'KCl_Kads_assumed_m3_kg':0.,
       'KCl_pore_diffusivity_assumed_m2_s':1.3e-9,
       'synthetic_noise_sigma_normalized_concentration':experiment['sigma'],
       'assumed_fisher_sigma':.015,
       'synthetic_design_seed':experiment['seed'],
       'Fisher_sample_count_per_tracer':int(np.sum(sel)),
       'fisher_diagnostics':designs,
       'hypothetical_robustness_grid':robustness,
       'single_NaGl_synthetic_fit':{
         'pore_access':float(np.exp(experiment['na_fit'].x[0])),
         'Kads':float(np.exp(experiment['na_fit'].x[1]))},
       'dual_synthetic_fit':{
         'pore_access':float(np.exp(experiment['dual_fit'].x[0])),
         'Kads':float(np.exp(experiment['dual_fit'].x[1]))},
       'single_fit_success':bool(experiment['na_fit'].success),
       'dual_fit_success':bool(experiment['dual_fit'].success),
       'single_parameter_log_error':float(np.linalg.norm(experiment['na_fit'].x-truth)),
       'dual_parameter_log_error':float(np.linalg.norm(experiment['dual_fit'].x-truth))},
      'first_moment_conditional_example':{
        'assumed_bed_void_fraction':porosity,
        'published_liquid_holdup':hL,'published_wetting_efficiency':fe,
        'calculated_NaGl_mean_time_s':mu_nagl,
        'calculated_assumed_nonadsorbing_KCl_mean_time_s':mu_kcl,
        'calculated_mean_delay_difference_s':delta,
        'warning':'Computed from Eq.16-style capacity under assumed bed void fraction and KCl adsorption=0; NOT experimental mean times.'},
      'limitations':['KCl data series unavailable in published source or verified accessible data.','Synthetic experiments use the same forward model for data generation and fitting (inverse crime).',
         'Estimated CRLB assumes independent Gaussian noise at selected timepoints; digitized traces are correlated.',
         'NaGl source fit uses only one original flow and is not externally validated.',
         'Adsorption relaxation rate, bed porosity, and accessible pore fraction lack independent experimental constraints.',
         'No inference of real KCl adsorption kinetics or transport parameters is possible.',
         'Not a reactor thermal model; no claim of full Alopaeus CFD reproduction.']}
    (OUT/'metrics.json').write_text(json.dumps(data,indent=2,allow_nan=False))
    with (OUT/'fig5_real_vs_fitted.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s','NaGl_inlet_graphical','NaGl_outlet_graphical','NaGl_heterogeneous_fit']);w.writerows(zip(t,ci,obs,fit['prediction']))
    with (OUT/'hypothetical_two_tracers.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s','NaGl_synthetic','NaGl_synthetic_noisy','KCl_hypothetical','KCl_synthetic_noisy']);w.writerows(zip(t,experiment['na_true'],experiment['na_noisy'],experiment['kcl_true'],experiment['kcl_noisy']))
    with (OUT/'fisher_summary.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['design','condition_number','min_eigenvalue','logdet','parameter_correlation','sigma_log_pore','sigma_log_Kads']);
        for name,r in designs.items():w.writerow([name,r['condition_number'],r['smallest_eigenvalue'],r['logdet'],r['log_parameter_correlation'],*r['optimistic_crlb_logparam_sd']])
    with (OUT/'local_fisher_robustness.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(robustness[0]))
        wr.writeheader();wr.writerows(robustness)
    fig,ax=plt.subplots(figsize=(9,5))
    ax.plot(t,ci,':',color='gray',label='Sales Fig.5 NaGl inlet (graphical)')
    ax.scatter(t,obs,s=13,color='black',label='Sales Fig.5 NaGl outlet (graphical)')
    ax.plot(t,fit['prediction'],lw=1.8,label='Heterogeneous v2.9 fit — same experimental curve')
    ax.set(xlim=(0,1500),ylim=(-.02,1.05),xlabel='Time (s)',ylabel='Normalized tracer signal',title='Real published Figure 5: exploratory in-sample NaGl fit')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'NaGl_experimental_fit.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    ax.plot(t,experiment['na_true'],label='Synthetic NaGl',lw=2)
    ax.plot(t,experiment['kcl_true'],label='Hypothetical nonadsorbing KCl',lw=2)
    ax.scatter(t[sel],experiment['kcl_noisy'][sel],marker='x',s=18,label='KCl synthetic noisy points')
    ax.set(xlim=(0,1500),xlabel='Time (s)',ylabel='Normalized tracer signal',title='CONTROLLED SYNTHETIC EXAMPLE — NOT EXPERIMENTAL KCl DATA')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'synthetic_dual_tracer.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.7))
    names=list(designs)
    x=np.arange(len(names))
    ax.bar(x,[designs[n]['smallest_eigenvalue'] for n in names])
    ax.set(xticks=x,xticklabels=['NaGl','Repeat NaGl','NaGl + KCl*'],ylabel='Smallest eigenvalue of Fisher matrix',title='Hypothetical local identifiability (* KCl simulated)')
    ax.grid(axis='y',alpha=.2);fig.tight_layout();fig.savefig(OUT/'fisher_information.png',dpi=170);plt.close(fig)
    rpt=f'''# TrickleBed-Dynamics v2.9 — Dual-tracer identifiability: evidence-first study

## Source and data availability

Sales et al., *Experimental dynamic evaluation of three-phase reactors*, Brazilian Journal of Chemical Engineering 22(3), 443–452 (2005); DOI 10.1590/S0104-66322005000300012.
The paper mentions an auxiliary **KCl** tracer in the discussion of Eq. (16) on printed page 449. However, its PDF reports **no numerical time-series for KCl**. Figure 5 on printed page 450 presents sodium gluconate experimental inlet and outlet curves at QL=5 L/h; Figure 6 shows sodium gluconate experimental outlet at QL=10 L/h, not KCl. Table 3 reports model-optimized hydrodynamic parameters, not raw KCl response points. We searched for separately released KCl curves but did not verify any publicly accessible original data file. **No experimental KCl data are contained in this package.**

The Figure 5 time series are PDF-vector-derived graphical resampling from v2.6, NOT independent apparatus measurements. No assumption of 59 independent measurements is justified.

## Real sodium gluconate experiment: exploratory fit only

With a two-parameter in-sample fit of effective pore accessibility and sodium gluconate adsorption coefficient (other constants taken from v2.8), the fitted pore accessibility is **{np.exp(fit['result'].x[0]):.4g}**, Kads is **{np.exp(fit['result'].x[1]):.4g} m³/kg** and full-figure graphical RMSE is **{fit['rmse_full']:.5f}** (fit transition RMSE **{fit['rmse_transit']:.5f}**). This is *not independent validation*. In particular the graph-derived data and unverified model closures do not allow interpreting a best-fit accessibility factor as physically measured pore availability. The model uses the finite-volume approximation motivated by the article's Eqs. 7–9, not a validated reproduction of its Laplace transfer functions (13–15).

## Controlled hypothetical dual-tracer experiment (synthetic only)

We assume *for illustration* that KCl is **non-adsorbing** (Kads=0) and has intraparticle D=1.3e-9 m²/s, while NaGl adsorption is given by a positive coefficient. These are NOT published KCl measurements. Both signals are generated by the v2.8 numerical model with added independent Gaussian synthetic noise, then fitted to the same model. This *inverse crime* makes the results only a demonstration of methodology, not a realistic forecast of industrial parameter accuracy.

We ask whether one sodium gluconate trace could separate two log parameters better if an additional tracer with different adsorption were measured, relative to (i) sodium gluconate alone and (ii) taking a second independent replicate of sodium gluconate. A local Fisher-information calculation uses {np.sum(sel)} selected nominal observation times per tracer and an assumed independent normalized-signal sigma of 0.015.

| Synthetic design | Smallest Fisher eigenvalue | Fisher condition number | Correlation of parameters |
|---|---:|---:|---:|
'''
    for name,r in designs.items():
        rpt+=f"| {name} | {r['smallest_eigenvalue']:.3g} | {r['condition_number']:.3g} | {r['log_parameter_correlation']:.3f} |\n"
    rpt+=f'''
The Fisher calculation is **conditional** on a point estimate, the assumption of known KCl transport properties and independent Gaussian noise. The optimistic Cramér–Rao bounds are not uncertainty estimates for the actual published graph.

### Synthetic parameter recovery (one seeded experiment)

- Synthetic true access fraction: **0.20**; NaGl Kads **8.2e-5 m³/kg**.
- Single-tracer synthetic fit: access **{np.exp(experiment['na_fit'].x[0]):.4g}**, Kads **{np.exp(experiment['na_fit'].x[1]):.4g}**.
- Dual-tracer synthetic fit: access **{np.exp(experiment['dual_fit'].x[0]):.4g}**, Kads **{np.exp(experiment['dual_fit'].x[1]):.4g}**.
- Log-parameter error norm: single **{np.linalg.norm(experiment['na_fit'].x-truth):.4g}**, dual **{np.linalg.norm(experiment['dual_fit'].x-truth):.4g}**.

These illustrative errors are one seeded realization; no statistical superiority claim is made based on a single run.

### Hypothetical local sensitivity over nine parameter assumptions

The model also recomputes Fisher diagnostics over a 3 × 3 grid of assumed accessibility and NaGl adsorption constants, without changing or adding experimental data. The minimum-eigenvalue gain ranges from **{min(r['min_eigen_gain'] for r in robustness):.2g}×** to **{max(r['min_eigen_gain'] for r in robustness):.2g}×** over this selected grid. The sampled range does **not** establish global robustness; all values are conditional on a hypothetical KCl response, known auxiliary transport parameters, and the same forward-model family. See `local_fisher_robustness.csv`.

### First moment linkage to Sales Eq. (16)

Assuming bed void fraction {porosity:.2f} (NOT measured in the paper), NaGl Kads=8.2e-5 m³/kg and hypothetically non-adsorbing KCl, Eq.(16)-style storage estimates at QL=5 L/h yield **{mu_nagl:.1f} s** for NaGl and **{mu_kcl:.1f} s** for KCl: difference **{delta:.1f} s**. This is a numerical illustration derived from assumed tracer properties, **not experimentally measured residence times**. A real analysis requires KCl adsorption equilibrium and injection/measurement response as obtained in a separate experiment.

## What would make this publishable / industrially useful?

1. Original separately measured KCl inlet/outlet time series and noise/calibration metadata.
2. Independent KCl and NaGl adsorption equilibrium measurements (or a validated physical basis for assuming KCl nonadsorbing).
3. Replicates and a different-flow holdout, preferably measured on the same bed.
4. Joint model discrepancy/measurement noise estimation, rather than interpreting graph-resampling as independent samples.
5. Comparison to full published transfer function (Sales Eqs. 13–15) and the earlier Alopaeus model without claiming exact CFD reproduction.

## Reproduce

`python dual_tracer_v29.py` and `python -m pytest -q` (run from project root).
'''
    (OUT/'REPORT.md').write_text(rpt)
    print('Experimental NaGl RMSE',fit['rmse_full'])
    print('Fisher diagnostics:')
    for k,v in designs.items():print(k,v)
    print('Synthetic recovered',np.exp(experiment['na_fit'].x),np.exp(experiment['dual_fit'].x))
    print('Report',OUT/'REPORT.md')
    return data


if __name__=='__main__':
    main()
