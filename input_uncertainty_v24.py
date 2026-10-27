"""v2.4: Transparent sensitivity audit of Figure 6 input assumptions.
Figure 6 contains a decay curve but no machine-readable inlet trace.
All observed points below are inherited approximate MANUAL readings (v2.0).
No pixel calibration or raw data are claimed.
"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.stats import gamma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/v24_input_uncertainty'

def predict_step(t,shape,mean,transition_start=0.,transition_duration=0.):
    """Outlet fraction for input 1->0, linearly ramping in duration;
    convolve step survival using quadrature when duration>0."""
    t=np.asarray(t,dtype=float)
    if shape<=0 or mean<=0 or transition_duration<0: raise ValueError('invalid parameters')
    if transition_duration==0:
        return gamma.sf(np.maximum(t-transition_start,0),a=shape,scale=mean/shape) * (t>=transition_start)+(t<transition_start)
    # uniform mixture of shifted downward steps
    offsets=np.linspace(transition_start,transition_start+transition_duration,41)
    return np.mean([predict_step(t,shape,mean,float(o),0) for o in offsets],axis=0)

def run():
    OUT.mkdir(parents=True,exist_ok=True)
    d=np.genfromtxt(ROOT/'results/v20_sales_digitization/figure6_approx_experimental.csv',delimiter=',',skip_header=1)
    t,y=d[:,0],d[:,1]
    prev=json.loads((ROOT/'results/v23_cross_flow/metrics.json').read_text())
    shape=prev['figure5_shape']; mean_scaled=prev['predicted_figure6_mean_s']; mean_unscaled=prev['figure5_gamma_mean_s']
    configs=[('scaled_ideal_step',mean_scaled,0.,0.),('scaled_shift_60s',mean_scaled,60.,0.),('scaled_ramp_120s',mean_scaled,0.,120.),('scaled_shift_60_ramp_120',mean_scaled,60.,120.),('unscaled_ideal_step',mean_unscaled,0.,0.)]
    rows=[]; predictions={}
    for name,mean,start,duration in configs:
        pred=predict_step(t,shape,mean,start,duration)
        predictions[name]=pred
        rows.append({'scenario':name,'gamma_mean_s':mean,'input_start_s':start,'input_ramp_s':duration,'rmse_approx':float(np.sqrt(np.mean((pred-y)**2))), 'mae_approx':float(np.mean(abs(pred-y)))})
    # A bounded model sensitivity envelope: ranges are illustrative, not experimental confidence intervals
    envelope=[]
    for start in [0,30,60,90]:
      for dur in [0,60,120]:
        envelope.append(predict_step(t,shape,mean_scaled,start,dur))
    envelope=np.array(envelope)
    with (OUT/'scenario_metrics.csv').open('w',newline='') as f:
      w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with (OUT/'predictions.csv').open('w',newline='') as f:
      w=csv.writer(f);w.writerow(['time_s_manual_approx','figure6_manual_approx']+list(predictions)+['assumption_envelope_min','assumption_envelope_max'])
      for i,ti in enumerate(t): w.writerow([ti,y[i]]+[float(p[i]) for p in predictions.values()]+[float(envelope[:,i].min()),float(envelope[:,i].max())])
    fig,ax=plt.subplots(figsize=(9,5))
    ax.plot(t,y,'k^',label='Published Fig 6: previous manual approx.',ms=5)
    ax.plot(t,predictions['scaled_ideal_step'],lw=2,label='scaled; ideal step')
    ax.plot(t,predictions['scaled_shift_60_ramp_120'],lw=2,label='scaled; 60 s shift + 120 s ramp')
    ax.plot(t,predictions['unscaled_ideal_step'],'--',label='unscaled; ideal step')
    ax.fill_between(t,envelope.min(axis=0),envelope.max(axis=0),alpha=.2,label='assumption-only envelope')
    ax.set(xlabel='Time (s)',ylabel='Normalized outlet fraction',title='Figure 6: sensitivity to undocumented input forcing')
    ax.set_ylim(-.03,1.08);ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(OUT/'input_sensitivity.png',dpi=170);plt.close(fig)
    metrics={'publication':'Sales et al. (2005), DOI 10.1590/S0104-66322005000300012','observations':'manual approximate Fig 6 readings inherited from v2.0','n_points':len(t),'raw_experiment_available':False,'pixel_calibrated':False,'inlet_trace_Fig6_available':False,'source_Fig6_caption': 'QG=300 L/h, QL=10 L/h','model':'Gamma transfer-function from Fig 5; mean adjusted by holdup/velocity ratio from Table 3','shape':shape,'scaled_mean_s':mean_scaled,'ranges':'input offset 0–90 s; ramp duration 0–120 s: exploratory assumed range, not measurement','scenarios':rows,'best_rmse_scenario':min(rows,key=lambda x:x['rmse_approx'])['scenario'],'warning':'Selection based on manually approximated values is not experimental model validation; input alternatives have no measured supporting trace.'}
    (OUT/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    report="""# v2.4 — Published-figure input uncertainty audit

Source: Sales, Maranhão, Pereira & Abreu (2005), *Experimental dynamic evaluation of three-phase reactors*, https://doi.org/10.1590/S0104-66322005000300012.

## Findings grounded in publication
- Figure 5 displays experimental **entrance**, experimental **exit** and theoretical **exit** for liquid flow **5 L/h**, gas flow **300 L/h**.
- Figure 6 displays experimental and modeled **outlet** for liquid flow **10 L/h**, gas flow **300 L/h**; it does not provide a separate sampled inlet curve for this run.
- The plotted Figure 6 response is decaying, consistent with a step-down-type experimental protocol. The exact forcing time course is not recoverable from the published graph.
- Published Table 3 reports fitted RTD hydrodynamic parameters, not a downloadable time-series dataset.

## What is computed here
- Retains the Gamma parameters fitted to **manual approximations of Figure 5** in v2.2. No refit on Figure 6.
- Transfers the mean residence scale using the published holdup/velocity relationship.
- Compares hypothesized input start times and gradual ramps to the Figure 6 **manual** readings.
- Creates an envelope spanning assumed timings; **not** statistical, pixel, or experimental uncertainty.

## Key caveats
- No accurate pixel-calibrated digitization has been performed in this release. The source PDF exists on the public journal site but was not available as a container file.
- There are no independently measured Figure 6 inlet signal samples available in our data files.
- Any improved RMSE obtained by choosing an input assumption from Figure 6 is an in-sample sensitivity diagnostic, **not** held-out validation.
- This is a separate RTD surrogate, **not** a validation of the reactive multiphase CFD equations.

## Next action to enable genuine cross-flow validation
Obtain author-supplied entrance and outlet time series (both flow rates) or a high-resolution journal PDF that permits reproducible pixel calibration. Archive the raw pixel-to-axis map, digitization error, and sampling protocol; fit only the 5 L/h run and freeze all parameters before predicting 10 L/h.
"""
    (OUT/'REPORT.md').write_text(report)
    print(json.dumps({'scenarios':rows,'tests':len(t),'output':str(OUT)},indent=2))
    return metrics
if __name__=='__main__':run()
