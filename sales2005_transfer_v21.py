"""Exploratory inlet-to-outlet RTD fit from APPROXIMATELY hand-read Sales 2005 Fig 5.
Not raw laboratory data and not pixel-calibrated digitization.
Run python sales2005_transfer_v21.py
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammainc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'/'v21_transfer';OUT.mkdir(parents=True,exist_ok=True)
DATA=ROOT/'results'/'v20_sales_digitization'/'figure5_approx_experimental.csv'
def gamma_cdf(t, mean_s, shape):
    a=np.asarray(t);return gammainc(shape,np.maximum(a,0)*shape/mean_s)
def transfer(t, inlet_t, inlet_y, tau, n, delay, dt=2.):
    """Convolve measured step-down inlet with residence-time distribution.
    Starting from steady state y=1, use convolution of the INPUT DROP derivative.
    A piecewise-linear interpolated inlet is assumed; kernel is normalized on grid.
    """
    time=np.asarray(t,float);fine=np.arange(0,max(time.max()+max(3*tau,200),2400)+dt,dt)
    inlet=np.interp(fine,inlet_t,inlet_y,left=inlet_y[0],right=inlet_y[-1]);
    # numerical time derivative; use gamma bin probabilities to avoid singular density at zero
    drops= -np.diff(inlet,prepend=inlet[0]); drops[0]=0
    bins=gamma_cdf(fine+dt/2-delay,tau,n)-gamma_cdf(fine-dt/2-delay,tau,n)
    bins=np.maximum(bins,0);bins/=bins.sum()
    shifted=np.convolve(drops,bins,mode='full')[:len(fine)]
    outlet=inlet_y[0]-np.cumsum(shifted)
    return np.interp(time,fine,outlet)
def main():
    a=np.genfromtxt(DATA,delimiter=',',skip_header=1)
    t,cin,cout=a.T
    # pointwise reading errors from v2.0 were assumptions, not instrument errors
    fit=least_squares(lambda x: (transfer(t,t,cin,*x)-cout),[120,3,50],bounds=([10,0.6,0],[900,30,500]),max_nfev=160,xtol=1e-8,ftol=1e-8)
    fitted=transfer(t,t,cin,*fit.x)
    baseline=cin.copy()
    metrics={'provenance':'Sales et al. 2005 Figure 5, approx hand-read v2.0 CSV; not raw or calibrated-pixel digitization','conditions':'Figure 5 Q_L=5 L/h, Q_G=300 L/h; do not transfer params to Figure 6 without recalibration','tau_mean_s':float(fit.x[0]),'gamma_shape':float(fit.x[1]),'delay_s':float(fit.x[2]),'rmse_approx_points':float(np.sqrt(np.mean((fitted-cout)**2))),'rmse_zero_residence_baseline':float(np.sqrt(np.mean((baseline-cout)**2))),'fit_npoints':len(t),'estimated_manual_time_uncertainty_s':45,'estimated_manual_y_uncertainty':0.06,'validation_status':'method demonstration against approximate graph readings only'}
    (OUT/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    with (OUT/'fig5_transfer_fit.csv').open('w',newline='') as f:
        w=__import__('csv').writer(f);w.writerow(['time_s_approx','inlet_approx','outlet_approx','predicted_from_inlet','zero_residence_baseline']);w.writerows(zip(t,cin,cout,fitted,baseline))
    tf=np.linspace(0,1800,500)
    fig,ax=plt.subplots(figsize=(9,5));ax.plot(t,cin,'o',label='Approx. experimental inlet');ax.plot(t,cout,'s',label='Approx. experimental outlet');ax.plot(tf,transfer(tf,t,cin,*fit.x),'-',label='Convolution model fitted to outlet');ax.set(xlabel='Time (s)',ylabel='Normalized tracer signal',ylim=(-.03,1.06),title='Sales (2005) Fig. 5 — exploratory inlet-to-outlet fit');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(OUT/'figure5_inlet_to_outlet.png',dpi=165);plt.close(fig)
    # Stress test: reproducibility and sensitivity to manual reading uncertainty, NOT frequentist CI
    rng=np.random.default_rng(2026);params=[]
    for i in range(35):
        perturb_in=np.clip(cin+rng.normal(0,.03,len(cin)),0,1); perturb_out=np.clip(cout+rng.normal(0,.03,len(cout)),0,1)
        perturb_in[0]=1;perturb_out[0]=1
        res=least_squares(lambda x:transfer(t,t,perturb_in,*x)-perturb_out,fit.x,bounds=([10,.6,0],[900,30,500]),max_nfev=90)
        params.append(res.x)
    params=np.array(params)
    metrics['stress_test_35_reading_perturbations_10_90_percentiles']={k:[float(z) for z in np.percentile(params[:,j],[10,90])] for j,k in enumerate(['tau_mean_s','shape','delay_s'])}
    (OUT/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    (OUT/'REPORT.md').write_text(f"""# Sales 2005 Figure 5 transfer-function study (v2.1)\n\n**Source:** Sales et al. (2005), DOI 10.1590/S0104-66322005000300012, Fig. 5 (trickle bed, QL=5 L/h, QG=300 L/h).\n\n**Data status:** CSV inherited from v2.0 consists of {len(t)} manually approximated points from a printed graph. NOT original data, NOT calibrated WebPlotDigitizer points. Figure 6 uses a different flow (QL=10 L/h) and is not a held-out validation of these parameters.\n\n## Method\n\nStarting at normalized steady concentration 1, the measured inlet step-down is treated as an exogenous input. Its negative numerical derivative is convolved with a delayed gamma residence-time kernel to predict the outlet response. Gamma CDF bin differences enforce kernel non-negativity and unit integral. Fit uses scipy.optimize.least_squares.\n\n## Fit to approximate points\n\n- Gamma mean residence time: {fit.x[0]:.3f} s\n- Gamma shape: {fit.x[1]:.3f}\n- Additional delay: {fit.x[2]:.3f} s\n- RMSE on approximate points: {metrics['rmse_approx_points']:.5f}\n- Zero-residence (outlet=inlet) RMSE: {metrics['rmse_zero_residence_baseline']:.5f}\n\n**Not physical validation:** These estimates are highly sensitive to hand-reading, baseline, image quality and model misspecification. The 35 perturbations are a stress test, NOT confidence intervals. This reduced-order residence-time kernel is distinct from the published multi-parameter diffusion/adsorption model and the earlier reactive four-compartment model.\n\n## Next requirements\n\nDigitize original figure using pixel-calibrated x/y axes; report point selection and extraction errors. Fit directly to those points, compare competing residence-time kernels, and test across another condition with physically justified flow scaling. Full reactor validation also requires thermal, kinetic and phase data.\n""")
    print(json.dumps(metrics,indent=2))
if __name__=='__main__':main()
