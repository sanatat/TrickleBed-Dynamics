from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rtd_calibration import *

def main():
    out=Path('results/v13_rtd');out.mkdir(parents=True,exist_ok=True)
    truth,d,noise=make_synthetic_data()
    train=d['calibration'];test=d['holdout']; t=train['time_s']
    fitted,opt=fit_rtd(t,train['observed'])
    baseline,_=fit_rtd(t,train['observed'],tanks_fixed=1.0)
    boot=residual_bootstrap(t,train['observed'],fitted,runs=100)
    qlo,qhi=np.percentile(boot,[2.5,97.5],axis=0)
    hold_preds=np.stack([impulse_response(t,RTDParams(*x),1.35) for x in boot]);low,high=np.percentile(hold_preds,[2.5,97.5],axis=0)
    pred_train=impulse_response(t,fitted)
    pred_hold=impulse_response(t,fitted,1.35)
    pred_base=impulse_response(t,baseline,1.35)
    summary={'study_type':'SYNTHETIC proof-of-concept; no experimental data','data_seed':20261026,'relative_flow_calibration':1.0,'relative_flow_holdout':1.35,'truth':vars(truth),'fit':vars(fitted),'baseline_single_tank_fit':vars(baseline),'fit_95pct_bootstrap_intervals':dict(zip(vars(fitted).keys(),zip(qlo.tolist(),qhi.tolist()))),'rmse_calibration':rmse(train['observed'],pred_train),'rmse_holdout':rmse(test['observed'],pred_hold),'rmse_holdout_single_tank':rmse(test['observed'],pred_base),'holdout_fraction_within_bootstrap_band':float(np.mean((test['observed']>=low)&(test['observed']<=high))),'bootstrap_runs':100,'noise_sigma':noise,'caveat':'Bootstrap bands cover parametric estimation uncertainty only and are not observation prediction intervals. Ground truth deliberately follows fitted model family.'}
    (out/'metrics.json').write_text(json.dumps(summary,indent=2))
    with (out/'curves.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s','calibration_synthetic','calibration_fit','holdout_synthetic','holdout_predicted','holdout_single_tank','holdout_param_low','holdout_param_high'])
        w.writerows(zip(t,train['observed'],pred_train,test['observed'],pred_hold,pred_base,low,high))
    with (out/'bootstrap_parameters.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['mean_residence_s','effective_tanks','delay_s']);w.writerows(boot)
    fig,ax=plt.subplots(figsize=(9,5))
    ax.scatter(t,train['observed'],s=8,label='Synthetic calibration observations',alpha=.65)
    ax.plot(t,pred_train,label='RTD model fit',lw=2)
    ax.plot(t,impulse_response(t,truth),':',label='Synthetic truth')
    ax.set(xlabel='Time (s)',ylabel='Outlet impulse RTD E(t) (1/s)',title='Calibration on synthetic RTD at relative flow 1.0')
    ax.legend();fig.tight_layout();fig.savefig(out/'calibration.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    ax.scatter(t,test['observed'],s=9,label='Held-out synthetic observations',alpha=.65)
    ax.plot(t,pred_hold,label='Unrefitted RTD prediction',lw=2)
    ax.plot(t,pred_base,'--',label='Single-CSTR baseline',lw=1.5)
    ax.fill_between(t,low,high,alpha=.25,label='Bootstrap parameter band (not predictive interval)')
    ax.set(xlabel='Time (s)',ylabel='Outlet impulse RTD E(t) (1/s)',title='Holdout at 35% higher flow (no refitting)');ax.legend();fig.tight_layout();fig.savefig(out/'holdout.png',dpi=160);plt.close(fig)
    report=f'''# RTD calibration & uncertainty: synthetic demonstration (v1.3)

**Not experimental calibration, not a validated reactor model, and not reproduced 2012 Alopaeus CFD data.** No experimental RTD dataset was available in the repository. Results must never be presented as measured or as demonstrated generalization to real equipment.

## Protocol
- Generate tracer impulse profiles from a gamma tanks-in-series model, with known synthetic truth `tau=85 s`, `N=5.2`, `delay=7 s` and deterministic random Gaussian measurement noise (seed 20261026).
- Fit `tau`, effective number of tanks `N` and delay on synthetic flow `Q/Qref=1.0` only by constrained nonlinear least squares.
- Predict the **independently generated**, previously unseen dataset at `Q/Qref=1.35` with **no refitting**. Assume tau scales with inverse flow while N and delay stay constant; this extrapolation law is an assumption.
- Compare to a single-CSTR (`N=1`) constrained baseline, fitted only on calibration data.
- Residual bootstrap (`n=100`, deterministic seed) gives approximate confidence bands for parameters. Plot shaded bands are *parameter uncertainty*, not observation intervals and not calibrated predictive uncertainty.

## Results
- RTD model: tau={fitted.mean_residence_s:.3f} s, N={fitted.effective_tanks:.3f}, delay={fitted.delay_s:.3f} s.
- Calibration RMSE={summary['rmse_calibration']:.7f} s^-1.
- Held-out RMSE={summary['rmse_holdout']:.7f} s^-1.
- Held-out single-CSTR RMSE={summary['rmse_holdout_single_tank']:.7f} s^-1.
- 95% residual bootstrap parameter bands (in order: tau, N, delay): {list(zip(qlo.round(3),qhi.round(3)))}.

## Interpretation and limitations
This is a **method/pipeline check** with ground truth intentionally in the fitted model family. Strong held-out performance is expected and is not scientific proof of physical generalization. The gamma RTD is a reduced-order hydraulic surrogate with no reaction, wetting mechanics, explicit catalyst, heat, or gas dynamics. It is NOT yet coupled to the v1.2 four-compartment simulation. Data must be replaced by independent published or experimental tracer observations for external validation, with sample-time and excitation protocol recorded. Parameter identifiability requires further checking on real measurements. No claim of experimental validation is supported.

## Reproduction
`python rtd_experiment_v13.py` then `python -m pytest -q`.
'''
    (out/'RTD_REPORT.md').write_text(report)
    print(json.dumps({k:v for k,v in summary.items() if k.startswith('rmse') or k=='fit'},indent=2))
if __name__=='__main__': main()
