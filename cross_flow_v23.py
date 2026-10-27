"""Cross-flow, NO-refit diagnostic using published holdup and velocity.
Sales et al. (2005), DOI 10.1590/S0104-66322005000300012.
Figure 6 points are approximate MANUAL readings inherited from v2.0.
This is not pixel-calibrated digitization and not raw experimental data.
"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.stats import gamma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/v23_cross_flow'
def residence_time(length, holdup, superficial_velocity):
    if min(length,holdup,superficial_velocity)<=0: raise ValueError('positive values required')
    return length*holdup/superficial_velocity

def run():
    OUT.mkdir(parents=True,exist_ok=True)
    # Sales Table 3: at QG=300 L/h, corresponding uG=2.13e-2 m/s.
    tau5_hyd=residence_time(0.80,0.17,0.001)
    tau10_hyd=residence_time(0.80,0.21,0.002)
    scale=tau10_hyd/tau5_hyd
    params=json.loads((ROOT/'results/v22_temporal_holdout/metrics.json').read_text())['gamma']['params']
    mean_5=params['mean_s']; shape=params['shape']; delay5=params['delay_s']
    mean10=mean_5*scale
    # For ideal normalized step-down inlet, outlet survival = 1 - CDF of RTD.
    # Published Fig6 inlet shape not available as independently sampled data:
    # the step forcing used here is an EXPLICIT TEST ASSUMPTION.
    data=np.genfromtxt(ROOT/'results/v20_sales_digitization/figure6_approx_experimental.csv',delimiter=',',skip_header=1)
    t,y=data[:,0],data[:,1]
    pred=gamma.sf(np.maximum(t-delay5*scale,0),a=shape,scale=mean10/shape)
    baseline=gamma.sf(t,a=shape,scale=mean_5/shape)
    rmse=lambda p:float(np.sqrt(np.mean((p-y)**2)))
    outcome={'source':'Sales 2005 Figure 6 (manually approximated values inherited v2.0)',
    'figure6_experimental_raw':False,'pixel_calibrated':False,'refit_to_figure6':False,
    'figure6_inlet_assumption':'ideal normalized step-down at t=0 (unverified exact experimental inlet)',
    'length_m':0.8,'Q_gas_lph':300,'table3_Q5':{'uL_m_s':.001,'holdup':.17,'hydraulic_mean_s':tau5_hyd},
    'table3_Q10':{'uL_m_s':.002,'holdup':.21,'hydraulic_mean_s':tau10_hyd},
    'time_ratio':scale,'figure5_gamma_mean_s':mean_5,'figure5_shape':shape,
    'predicted_figure6_mean_s':mean10,'n_approx':len(t),
    'rmse_cross_flow_approx':rmse(pred),'rmse_unscaled_approx':rmse(baseline),
    'assessment':'diagnostic only; independent flow condition but approximate manually read points and ideal step forcing.'}
    (OUT/'metrics.json').write_text(json.dumps(outcome,indent=2)+'\n')
    with (OUT/'cross_flow.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s_approx','figure6_manually_approximated','cross_flow_prediction_no_refit','unscaled_baseline']);w.writerows(zip(t,y,pred,baseline))
    fig,ax=plt.subplots(figsize=(8.3,4.8));ax.plot(t,y,'^',ms=5,label='Sales Fig. 6 approximate readings');ax.plot(t,pred,'-',label='Fig. 5 gamma, holdup-scaled to 10 L/h');ax.plot(t,baseline,'--',label='Fig. 5 gamma, unscaled')
    ax.set(xlabel='Time (s)',ylabel='Normalized tracer fraction',title='Cross-flow diagnostic (not direct experimental validation)');ax.set_ylim(-.03,1.06);ax.grid(alpha=.22);ax.legend();fig.tight_layout();fig.savefig(OUT/'cross_flow.png',dpi=180);plt.close(fig)
    (OUT/'REPORT.md').write_text(f"""# v2.3 — Cross-flow RTD assessment (no refitting)\n\nSource: Sales et al. (2005), Table 3 and Figs 5–6. DOI: 10.1590/S0104-66322005000300012.\n\n## Published values vs assumptions\n- **Published fitted experimental-derived Table 3**: at gas velocity 0.0213 m/s, liquid velocity 0.001/0.002 m/s corresponds to holdup 0.17/0.21.\n- Bed height 0.80 m (published apparatus). Nominal hydraulic times `height × holdup / superficial liquid velocity`: {tau5_hyd:.1f} s (5 L/h) and {tau10_hyd:.1f} s (10 L/h), ratio {scale:.4f}. These are idealized estimates, not measured tracer means.\n- v2.2 Figure 5 Gamma mean {mean_5:.2f} s; shape {shape:.2f}, estimated from 13 **manually approximated** graphical readings. No Figure 6 fitting is performed.\n- Figure 6 input assumed to be an ideal step-down at t=0. This assumption cannot be confirmed from digital inlet measurements for that trial. Figure 6 points were manually estimated in v2.0, not pixel-calibrated or raw.\n\n## Held-out flow comparison (approximate, illustrative only)\n- Scaled model approximate RMSE = {rmse(pred):.5f} normalized fraction.\n- Unscaled baseline approximate RMSE = {rmse(baseline):.5f}.\n- The errors include image-reading and inlet-forcing mismatch, and cannot establish the accuracy of the physical reactor model.\n\n## Next evidence requirements\n1. Obtain the author's original time-concentration data or a high-resolution source plot.\n2. Digitize each experiment with saved pixel coordinates and axis calibration, and assess digitization repeatability.\n3. Get measured inlet signal for Figure 6; do not substitute a step unless documented.\n4. Validate cross-flow without refitting and report experimental measurement uncertainty.\n""")
    print(json.dumps(outcome,indent=2))
    return outcome
if __name__=='__main__':run()
