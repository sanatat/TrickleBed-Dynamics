"""Synthetic multi-flow tracer parameter estimation and identifiability.
Not experimental validation. Forward model is tracer_coupled_v14, whose
axial dispersion closure is also used in reactive_rtd_v15.
"""
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.stats import chi2
from v6_model import Params
from tracer_coupled_v14 import TracerConfig,run_tracer

TRUE=(0.0025,1.4e-5)
SEED=20261026

def setup(n=10):
    base=Params(n=n,t_step=10000,wet_init=.73,wet_adv=.73,wet_rec=.73,cap_tau=140)
    return base

def observed(p, t, theta, noise=0., seed=0):
    c=TracerConfig(k_exchange=theta[0],dispersion=theta[1],pulse_time=25.,pulse_sigma=5.)
    y=run_tracer(p,c,t)['outlet_conc']
    if noise: y=y+np.random.default_rng(seed).normal(0,noise,len(t))
    return y

def estimate(cases, guess=(.002,1e-5)):
    # use log10 coordinates for scale-invariant conditioning
    def residual(x):
        th=10**x
        return np.concatenate([(observed(p,t,th)-y)/sigma for p,t,y,sigma in cases])
    fit=least_squares(residual,np.log10(guess),bounds=(np.log10([1e-5,1e-7]),np.log10([.03,3e-4])),max_nfev=45,ftol=1e-7,xtol=1e-7,gtol=1e-7)
    jac=np.asarray(fit.jac)
    u,sv,vt=np.linalg.svd(jac,full_matrices=False)
    corr=np.linalg.pinv(jac.T@jac)
    correlation=corr[0,1]/np.sqrt(corr[0,0]*corr[1,1]) if np.all(np.diag(corr)>0) else np.nan
    return fit,dict(estimate=(10**fit.x).tolist(),singular_values=sv.tolist(),jacobian_condition=float(sv[0]/sv[-1]),parameter_correlation=float(correlation),cost=float(fit.cost),success=bool(fit.success))

def run(out='results/v16_identifiability'):
    folder=Path(out);folder.mkdir(parents=True,exist_ok=True)
    t=np.linspace(0,240,61)
    p1=setup();p2=replace(p1,u_L=p1.u_L*1.35)
    # Different random streams; no reusing held-out observations in fit.
    noise=8e-5
    y1=observed(p1,t,TRUE,noise,seed=SEED)
    y2=observed(p2,t,TRUE,noise,seed=SEED+1)
    fit1,diag1=estimate([(p1,t,y1,noise)])
    fit2,diag2=estimate([(p1,t,y1,noise),(p2,t,y2,noise)],guess=tuple(diag1['estimate']))
    # Hold-out flow 60% larger: generated independently and never used in calibration.
    p3=replace(p1,u_L=p1.u_L*1.60)
    test=observed(p3,t,TRUE,noise,seed=SEED+2)
    pred1=observed(p3,t,diag1['estimate']);pred2=observed(p3,t,diag2['estimate'])
    baseline=observed(p3,t,(.0002,1e-7))
    def rmse(y): return float(np.sqrt(np.mean((y-test)**2)))
    metrics={'data_type':'fully_synthetic','true_parameters':{'k_exchange':TRUE[0],'D_ax_m2_s':TRUE[1]},'noise_sd':noise,'seed':SEED,'single_flow':diag1,'two_flow':diag2,'holdout':{'flow_multiplier':1.6,'rmse_single':rmse(pred1),'rmse_two':rmse(pred2),'rmse_weak_baseline':rmse(baseline)},'note':'Local SVD and correlation only; confidence intervals from inverse Hessian would require stronger assumptions. These are synthetic identifiability diagnostics, not experimental validation.'}
    (folder/'metrics.json').write_text(json.dumps(metrics,indent=2))
    np.savetxt(folder/'holdout.csv',np.column_stack((t,test,pred1,pred2,baseline)),delimiter=',',header='time_s,synthetic_holdout,single_flow_fit,two_flow_fit,weak_baseline',comments='')
    fig,ax=plt.subplots(figsize=(8,4.5));ax.plot(t,test,'.',ms=3,label='Synthetic holdout');ax.plot(t,pred1,label='Single-flow calibration');ax.plot(t,pred2,label='Two-flow calibration');ax.plot(t,baseline,':',label='Weak baseline');ax.set(xlabel='Time (s)',ylabel='Outlet tracer concentration (arb. units)',title='Untouched flow holdout (+60%)');ax.legend();fig.tight_layout();fig.savefig(folder/'holdout.png',dpi=160);plt.close(fig)
    summary=f"""# v1.6: multi-flow identifiability audit

All observations are **synthetic**, generated with the same tracer model family used for estimation. This is not real-data validation.

Reference synthetic values: exchange={TRUE[0]} (model exchange units), axial dispersion={TRUE[1]:.3g} m²/s. Gaussian noise SD={noise}.

Single-flow fit: {diag1['estimate']}; Jacobian condition: {diag1['jacobian_condition']:.3g}, local parameter correlation: {diag1['parameter_correlation']:.3g}.

Two-flow fit: {diag2['estimate']}; Jacobian condition: {diag2['jacobian_condition']:.3g}, local parameter correlation: {diag2['parameter_correlation']:.3g}.

Holdout RMSE (flow +60%): single {rmse(pred1):.7g}, multi {rmse(pred2):.7g}, weak baseline {rmse(baseline):.7g}.

## Interpretation

A lower local Jacobian condition indicates more balanced information for the two fitted log-parameters, but does not prove structural identifiability. The held-out flow uses the same generative model, so this tests numerical generalization under idealized assumptions, **not** physical predictiveness. Do not claim calibrated Ni/Al2O3 kinetics or published-data validation.

The forward RTD model and the reactive solver share a dispersion operator, but this script does not infer reaction kinetics or validate thermal behavior.
"""
    (folder/'REPORT.md').write_text(summary)
    print(summary)
    return metrics
if __name__=='__main__':run()
