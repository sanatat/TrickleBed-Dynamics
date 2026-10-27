"""Robust D-optimal flow design under parameter uncertainty.
All measurements are simulated; candidate model and generating model match.
Worst-case log(det(FIM)) is optimized over a pre-registered parameter ensemble.
"""
from dataclasses import replace
from itertools import combinations
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optimal_design_v17 import FLOWS, BASELINE, TRUE, SIGMA, EPS, sensitivities, score
from identifiability_v16 import setup, observed, estimate, SEED

# Encompass roughly 50-200% of nominal exchange & 0.5-2x dispersion.
ENSEMBLE = ((.00125,7e-6),(.00125,2.8e-5),(.0025,1.4e-5),(.005,7e-6),(.005,2.8e-5))
TIMES = np.linspace(0,240,31)
CANDIDATES = list(combinations(FLOWS,2))
HELD_FLOW = 1.75
NOISE_SEEDS = (110, 111, 112, 113, 114)

def build_information():
    bank=[]
    for theta in ENSEMBLE:
        jac={float(f):sensitivities(f,theta=theta,t=TIMES) for f in FLOWS}
        bank.append(jac)
    return bank

def design_records(bank):
    records=[]
    for pair in CANDIDATES:
        vals=np.array([score(pair,j) for j in bank])
        records.append(dict(flows=pair, worst=float(vals.min()),mean=float(vals.mean()),std=float(vals.std()), nominal=float(vals[len(vals)//2])))
    return records

def fit_pair(pair,theta,noise_seed):
    cases=[]
    for i,f in enumerate(pair):
        p=replace(setup(),u_L=setup().u_L*f)
        cases.append((p,TIMES,observed(p,TIMES,theta,noise=SIGMA,seed=SEED+noise_seed*100+i),SIGMA))
    fit,summary=estimate(cases)
    ptest=replace(setup(),u_L=setup().u_L*HELD_FLOW)
    data=observed(ptest,TIMES,theta,noise=SIGMA,seed=SEED+noise_seed*100+90)
    pred=observed(ptest,TIMES,summary['estimate'])
    rmse=float(np.sqrt(np.mean((pred-data)**2)))
    return rmse,summary

def run(out='results/v18_robust_design'):
    folder=Path(out);folder.mkdir(parents=True,exist_ok=True)
    bank=build_information();rec=design_records(bank)
    robust=max(rec,key=lambda r:(r['worst'],r['mean']))
    local=max(rec,key=lambda r:r['nominal'])
    base=next(r for r in rec if r['flows']==BASELINE)
    results=[]
    for theta_idx in (0,4):
        theta=ENSEMBLE[theta_idx]
        for seed in NOISE_SEEDS[:1]:
            for title,pair in [('baseline',BASELINE),('local',local['flows']),('robust',robust['flows'])]:
                rmse,diag=fit_pair(pair,theta,seed)
                results.append({'theta_idx':theta_idx,'theta_exchange':theta[0],'theta_dispersion':theta[1],'replicate':seed,'design':title,'flow1':pair[0],'flow2':pair[1],'holdout_rmse':rmse,'estimated_exchange':diag['estimate'][0],'estimated_dispersion':diag['estimate'][1],'jacobian_condition':diag['jacobian_condition'],'fit_success':diag['success']})
    means={title:float(np.mean([x['holdout_rmse'] for x in results if x['design']==title])) for title in ('baseline','local','robust')}
    wins={title:sum(x['holdout_rmse']<next(y['holdout_rmse'] for y in results if y['theta_idx']==x['theta_idx'] and y['replicate']==x['replicate'] and y['design']=='baseline') for x in results if x['design']==title) for title in ('local','robust')}
    metrics={'data_type':'entirely synthetic; identical forward and generating model','objective':'maximize minimum log-det Fisher information across a fixed 5-point parameter ensemble','parameter_ensemble':[list(t) for t in ENSEMBLE],'sampling_times_s':TIMES.tolist(),'assumed_gaussian_noise_sd':SIGMA,'baseline_design':list(BASELINE),'local_design':list(local['flows']),'robust_design':list(robust['flows']),'scores':{name:{k:v for k,v in r.items() if k!='flows'} for name,r in [('baseline',base),('local',local),('robust',robust)]},'holdout_flow':HELD_FLOW,'n_independent_noise_trials':len(NOISE_SEEDS[:1]),'heldout_synthetic_scenarios':[0,4],'holdout_mean_rmse':means,'wins_vs_baseline':wins,'limitations':['robust over finite preselected ensemble only','no published or real experimental RTD data','local Fisher information linearization within each ensemble member','optimization uses the same model family as synthetic tests','no physical safety/flow operating constraints','no guarantees of improved holdout RMSE']}
    (folder/'metrics.json').write_text(json.dumps(metrics,indent=2))
    with (folder/'scores.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['flow1','flow2','min_logdet','mean_logdet','std_logdet','nominal_logdet']);w.writerows([[*r['flows'],r['worst'],r['mean'],r['std'],r['nominal']] for r in sorted(rec,key=lambda x:x['worst'],reverse=True)])
    with (folder/'holdout_trials.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    fig,ax=plt.subplots(figsize=(8,4.5));ss=sorted(rec,key=lambda x:x['worst'],reverse=True)[:12];ax.barh([str(r['flows']) for r in ss[::-1]],[r['worst'] for r in ss[::-1]]);ax.set(xlabel='Worst-case log det FIM (five parameter scenarios)',ylabel='Flow multipliers',title='Robust design candidates');fig.tight_layout();fig.savefig(folder/'robust_ranking.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4.5));names=list(means);ax.bar(names,[means[n] for n in names]);ax.set(ylabel='Mean synthetic holdout RMSE',title='Stress test across 2 parameter sets × 1 noise seed');fig.tight_layout();fig.savefig(folder/'holdout_comparison.png',dpi=160);plt.close(fig)
    lines=['# v1.8 Robust optimal experimental design','', '**All observations are synthetic.** Neither empirical reactor calibration nor industrial validation is claimed.','',f'Baseline flows: {BASELINE}; local flows: {local["flows"]}; robust flows: {robust["flows"]}.','', '| Design | Worst-case log det FIM | Mean holdout RMSE |','|---|---:|---:|']
    for name,r in [('baseline',base),('local',local),('robust',robust)]:lines.append(f'| {name} | {r["worst"]:.4f} | {means[name]:.8g} |')
    lines += ['', 'Optimization: evaluate the local Fisher matrix on five explicitly selected parameter pairs and maximize the **minimum** log determinant; no test data participates in design selection.', 'Testing: two parameter settings and one Gaussian noise realization each, with a disjoint holdout flow of 1.75×; refit the two parameters for each simulated training pair.', '', 'This demonstrates a computational methodology, not physical predictiveness. Optimized Fisher information is not equivalent to guaranteed lower RMSE.', '', 'Important limitations: the 5-point scenario set is not a probability distribution or confidence region; true parameters and noise process are unmeasured; no operational feasibility constraint; one noise trial per scenario gives highly uncertain error summaries; wetting and chemistry are not newly calibrated.']
    (folder/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'designs':{k:v for k,v in [('baseline',BASELINE),('local',local['flows']),('robust',robust['flows'])]},'worst_logdet':{x:r['worst'] for x,r in [('baseline',base),('local',local),('robust',robust)]},'mean_holdout_rmse':means,'n_tests':len(results)},indent=2))
    return metrics
if __name__=='__main__':run()
