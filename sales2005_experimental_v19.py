"""Rebuild v1.9 published Table 3 experiment-derived regression and plots.
No raw tracer time series are included. See results/v19_sales_experimental/REPORT.md.
"""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path('results/v19_sales_experimental')
def run():
    rows=np.genfromtxt(OUT/'sales2005_table3.csv',delimiter=',',skip_header=1)
    X=np.column_stack([np.ones(6),rows[:,0]-2,(rows[:,1]-4.26)/2.13]); Xs=X[:,:2]
    records=[];plots=[]
    for col,name in [(2,'wetting_efficiency'),(3,'liquid_holdup'),(4,'Peclet_L')]:
        y=rows[:,col];p=[];s=[];b=[]
        for i in range(len(rows)):
            ids=np.arange(len(rows))!=i
            coef=np.linalg.solve(X[ids].T@X[ids]+np.diag([0,.03,.03]),X[ids].T@y[ids])
            coefs=np.linalg.solve(Xs[ids].T@Xs[ids]+np.diag([0,.03]),Xs[ids].T@y[ids])
            p.append(float(X[i]@coef));s.append(float(Xs[i]@coefs));b.append(float(y[ids].mean()))
        for i in range(len(rows)):
            records.append(dict(response=name,row=i+1,observed=y[i],prediction_two_flows=p[i],prediction_liquid_only=s[i],prediction_mean=b[i]))
        plots.append((name,y,np.array(p),np.array(s)))
    with (OUT/'leave_one_out_predictions.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    fig,ax=plt.subplots(1,3,figsize=(12,4))
    for a,(name,y,p,s) in zip(ax,plots):
        a.plot(y,y,label='Ideal');a.scatter(y,p,label='Two-flow LOO');a.scatter(y,s,marker='x',label='Liquid-only LOO');a.set(xlabel='Published fitted value',ylabel='LOO prediction',title=name);a.legend(fontsize=7)
    fig.tight_layout();fig.savefig(OUT/'leave_one_out.png',dpi=170);plt.close(fig)
    print(json.dumps({name:{'LOO_MAE_two_flows':float(np.mean(abs(y-p))),'LOO_MAE_liquid_only':float(np.mean(abs(y-s)))} for name,y,p,s in plots},indent=2))
if __name__=='__main__':run()
