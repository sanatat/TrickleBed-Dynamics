"""Local D-optimal flow design, simulated RTD experiments only.
Uses log-parameters and Fisher information J.T@J; pre-specified baseline.
No external experimental data, no industrial validation.
"""
from dataclasses import replace
from itertools import combinations
from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from identifiability_v16 import setup, observed, estimate, TRUE, SEED

FLOWS=(0.65,0.8,1.0,1.2,1.35,1.6,1.9)
BASELINE=(1.0,1.35)
SIGMA=8e-5
EPS=0.015

def sensitivities(flow, theta=TRUE, t=None):
    if t is None:t=np.linspace(0,240,61)
    p=replace(setup(),u_L=setup().u_L*flow)
    x=np.log10(theta)
    cols=[]
    for i in range(2):
        plus=x.copy();minus=x.copy();plus[i]+=EPS;minus[i]-=EPS
        yp=observed(p,t,10**plus)
        ym=observed(p,t,10**minus)
        cols.append((yp-ym)/(2*EPS*SIGMA))
    return np.column_stack(cols)

def information(design,jacs):
    j=np.concatenate([jacs[float(x)] for x in design],axis=0)
    return j.T@j

def score(design,jacs):
    info=information(design,jacs)
    sign,logdet=np.linalg.slogdet(info)
    return float(logdet) if sign>0 else float('-inf')

def design(out='results/v17_design'):
    folder=Path(out);folder.mkdir(parents=True,exist_ok=True)
    t=np.linspace(0,240,61)
    jacs={float(f):sensitivities(f,t=t) for f in FLOWS}
    candidates=list(combinations(FLOWS,2))
    scores=sorted([(score(c,jacs),c) for c in candidates],reverse=True)
    best=scores[0][1]
    def summary(flows,seed_shift):
        cases=[]
        for i,f in enumerate(flows):
            p=replace(setup(),u_L=setup().u_L*f)
            y=observed(p,t,TRUE,noise=SIGMA,seed=SEED+seed_shift+i)
            cases.append((p,t,y,SIGMA))
        fit,diag=estimate(cases)
        return fit,diag
    fit_base,base=summary(BASELINE,100)
    fit_opt,opt=summary(best,100)
    # Entirely unused test condition and independently generated noise.
    test_flow=1.75
    ptest=replace(setup(),u_L=setup().u_L*test_flow)
    test=observed(ptest,t,TRUE,noise=SIGMA,seed=SEED+10000)
    pred_base=observed(ptest,t,base['estimate'])
    pred_opt=observed(ptest,t,opt['estimate'])
    rmse=lambda v:float(np.sqrt(np.mean((v-test)**2)))
    fi_b=information(BASELINE,jacs); fi_o=information(best,jacs)
    metrics={'data_type':'synthetic','design_type':'local_D_optimal','theta_nominal':list(TRUE),'noise_sd':SIGMA,'flows_candidates':list(FLOWS),'baseline_flows':list(BASELINE),'optimal_flows':list(best),'baseline_logdet_FIM':score(BASELINE,jacs),'optimal_logdet_FIM':score(best,jacs),'logdet_gain':score(best,jacs)-score(BASELINE,jacs),'baseline_FIM_eigenvalues':np.linalg.eigvalsh(fi_b).tolist(),'optimal_FIM_eigenvalues':np.linalg.eigvalsh(fi_o).tolist(),'baseline_parameter_estimates':base['estimate'],'optimal_parameter_estimates':opt['estimate'],'baseline_local_condition':base['jacobian_condition'],'optimal_local_condition':opt['jacobian_condition'],'holdout_flow':test_flow,'holdout_rmse_baseline':rmse(pred_base),'holdout_rmse_optimal':rmse(pred_opt),'warning':'Local nominal optimality with fixed known noise, constant wetting, and same forward/generative model. One random trial; no claim that heldout RMSE must improve.'}
    (folder/'metrics.json').write_text(json.dumps(metrics,indent=2))
    with (folder/'design_scores.csv').open('w',newline='') as h:
        w=csv.writer(h);w.writerow(['flow1','flow2','logdet_fim'])
        for sc,(a,b) in scores:w.writerow([a,b,sc])
    with (folder/'holdout.csv').open('w',newline='') as h:
        w=csv.writer(h);w.writerow(['time_s','holdout_synthetic','baseline_fit','optimal_fit']);w.writerows(zip(t,test,pred_base,pred_opt))
    fig,ax=plt.subplots(figsize=(7.5,4.5))
    ax.plot(t,test,'.',label='Synthetic held-out data',ms=3)
    ax.plot(t,pred_base,label='Baseline flow design')
    ax.plot(t,pred_opt,label='D-optimal flow design')
    ax.set(xlabel='Time (s)',ylabel='Outlet tracer concentration (arbitrary units)',title='Held-out flow response (+75%)');ax.legend();fig.tight_layout();fig.savefig(folder/'holdout.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.5))
    chosen=scores[:12];ax.barh([f'{a:g}, {b:g}' for _,(a,b) in chosen[::-1]],[s for s,_ in chosen[::-1]])
    ax.axvline(metrics['baseline_logdet_FIM'],ls='--',color='gray',label='Baseline')
    ax.set(xlabel='log det (Fisher information)',ylabel='Flow multipliers',title='Local D-optimal flow selection');ax.legend();fig.tight_layout();fig.savefig(folder/'design_ranking.png',dpi=160);plt.close(fig)
    lines=['# v1.7 — Local optimal design for synthetic tracer experiments','','This is a reproducible **in-silico demonstration**, not experimental validation.','',f'Candidate flow multipliers: {FLOWS}',f'Pre-specified baseline: {BASELINE}',f'Optimal pair: {best}',f'Baseline log det FIM: {metrics["baseline_logdet_FIM"]:.4f}',f'Optimal log det FIM: {metrics["optimal_logdet_FIM"]:.4f}',f'Increase: {metrics["logdet_gain"]:.4f}',f'Fitted baseline log10 parameter conditioning: {base["jacobian_condition"]:.3f}',f'Fitted optimal log10 parameter conditioning: {opt["jacobian_condition"]:.3f}',f'Holdout RMSE, baseline: {rmse(pred_base):.8g}',f'Holdout RMSE, optimal: {rmse(pred_opt):.8g}', '', '## Limitations','- The design objective uses sensitivities around assumed true parameters; real-world nominal parameters are unknown.','- Fisher information is a local approximation that depends on sampling schedule and a known constant noise assumption.','- Testing uses one synthetic noise realization and the same model family used to generate the data.','- No external physical RTD data, catalyst kinetic calibration, phase equilibrium, or industrial thermal validation.','- A larger information determinant does not guarantee lower hold-out error in a given noisy experiment.','- The current flow multipliers are relative to a demonstration base superficial velocity; feasibility and safe operation are not evaluated.']
    (folder/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print('Design:',best,'baseline logdet:',round(metrics['baseline_logdet_FIM'],3),'best logdet:',round(metrics['optimal_logdet_FIM'],3),'RMSEs:',rmse(pred_base),rmse(pred_opt))
    return metrics
if __name__=='__main__': design()
