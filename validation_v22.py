"""Prospective-style temporal holdout and kernel comparison for approximate Sales Fig 5 points.
This is a numerical methodological demonstration, NOT experimental validation.
"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sales2005_transfer_v21 import transfer
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'/'v22_temporal_holdout'
SRC=ROOT/'results'/'v20_sales_digitization'/'figure5_approx_experimental.csv'

def rmse(x,y): return float(np.sqrt(np.mean((np.asarray(x)-np.asarray(y))**2)))
def predict(t,cin,tau,shape,delay):
    # Gamma convolution applied to entire INPUT trace, but calibrated only against early output observations.
    return transfer(t,t,cin,tau,shape,delay,dt=3.0)
def fit(t,cin,cout,mask,kind):
    if kind=='gamma':
        x0=[160,5,0]; lower=[10,0.7,0];upper=[700,25,350]
        unpack=lambda x:x
    else:
        x0=[160,0];lower=[10,0];upper=[700,350]
        unpack=lambda x:(x[0],1.0,x[1])
    res=least_squares(lambda x:predict(t,cin,*unpack(x))[mask]-cout[mask],x0,bounds=(lower,upper),max_nfev=90)
    return unpack(res.x),res

def run():
    OUT.mkdir(parents=True,exist_ok=True)
    t,cin,cout=np.genfromtxt(SRC,delimiter=',',skip_header=1).T
    # first 9 observations up to 750 s train; last 4 are strictly held out for OUTLET fitting.
    mask=t<=750
    metrics={'data_origin':'Sales et al. 2005 Fig 5; 13 manually approximated points inherited from v2.0, not original measurements or pixel-calibrated digitization',
             'train_until_s':750,'train_n':int(mask.sum()),'holdout_n':int((~mask).sum()),
             'caution':'The complete approximate inlet trace is used as known forcing; only late outlet observations are held out. Same single run, not independent experimental validation.'}
    curves=[]
    for kind in ['gamma','exponential']:
        pars,_=fit(t,cin,cout,mask,kind)
        pred=predict(t,cin,*pars)
        metrics[kind]={'params':{'mean_s':float(pars[0]),'shape':float(pars[1]),'delay_s':float(pars[2])},'train_rmse':rmse(pred[mask],cout[mask]),'holdout_rmse':rmse(pred[~mask],cout[~mask]),'all_rmse':rmse(pred,cout)}
        curves.append(pred)
    metrics['zero_residence_baseline']={'holdout_rmse':rmse(cin[~mask],cout[~mask]),'train_rmse':rmse(cin[mask],cout[mask])}
    # Input-only knowledge does not imply an independent test condition.
    with (OUT/'temporal_holdout.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s_approx','inlet_approx','outlet_approx','train_output_flag','gamma_prediction','exponential_prediction']);
        w.writerows(zip(t,cin,cout,mask.astype(int),*curves))
    fig,ax=plt.subplots(figsize=(9,5));ax.plot(t,cin,'o-',label='Approximate inlet',alpha=.6);ax.scatter(t[mask],cout[mask],marker='s',label='Training outlet');ax.scatter(t[~mask],cout[~mask],marker='D',label='Held-out outlet');ax.plot(t,curves[0],label='Gamma model');ax.plot(t,curves[1],label='Exponential model');ax.axvline(750,linestyle='--',color='gray',alpha=.6);ax.set(xlabel='Time (s)',ylabel='Normalized concentration',title='Sales 2005 Fig. 5 (approximate manual readings)');ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(OUT/'temporal_holdout.png',dpi=160);plt.close(fig)
    (OUT/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    report=f'''# v2.2 — Temporal holdout and model comparison

**Source**: Sales et al. (2005), Figure 5, DOI: 10.1590/S0104-66322005000300012.

**Critical provenance**: Data are 13 *approximate manually read points* inherited from v2.0, not author-supplied raw observations or pixel-calibrated coordinates. The assumed reading errors are not metrological uncertainty estimates.

## What was tested

- Predict outlet from approximate measured inlet by convolution with normalized residence-time kernel.
- Fit only outlet points through 750 s (9 points), and evaluate on four later points **from the same experimental curve**.
- Compare flexible gamma RTD and two-parameter delayed exponential RTD. This is a temporal split, not independent validation on another reactor or flow rate.

## Results

| Model | Training RMSE | Temporal holdout RMSE |
|---|---:|---:|
| Gamma | {metrics['gamma']['train_rmse']:.6f} | {metrics['gamma']['holdout_rmse']:.6f} |
| Exponential | {metrics['exponential']['train_rmse']:.6f} | {metrics['exponential']['holdout_rmse']:.6f} |
| No residence time | {metrics['zero_residence_baseline']['train_rmse']:.6f} | {metrics['zero_residence_baseline']['holdout_rmse']:.6f} |

**Interpretation limitations**: the inlet is known through the held-out interval, and each model sees early values of the same decay curve. With only four test points and uncertain digitization, these scores do not establish generalization or reliable physical parameter values. The models are empirical RTD kernels, not full CFD or the reactive catalyst-pore model.

## Next steps needed for experimental validation

1. Acquire high-resolution source figure and digitize all distinct measurement markers with pixel-to-axis calibration, archiving the selected coordinates and image.
2. Obtain raw data from authors if possible.
3. Measure/obtain separate inlet and outlet traces at a second flow rate; fit one condition and predict the other with a physically supported flow-scaling law.
4. Report RTD normalization, tracer recovery, instrument lag, and identifiability diagnostics.
'''
    (OUT/'REPORT.md').write_text(report)
    print(json.dumps(metrics,indent=2))
    return metrics
if __name__=='__main__':run()
