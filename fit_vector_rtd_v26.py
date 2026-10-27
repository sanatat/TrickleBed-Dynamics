"""Reduced-order input-output RTD fit to vector-extracted published Figure 5.
Figure 6 is summarized but NOT validated: its inlet waveform was not published.
"""
from pathlib import Path
import json,csv
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammaincc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sales2005_transfer_v21 import transfer
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'/'v26_vector_extraction'

def read_curve(name):
 a=np.loadtxt(OUT/f'{name}_vector.csv',delimiter=',',skiprows=1)
 return a[:,0],np.clip(a[:,1],0.,1.)

def quantile_crossing(t,c,level):
 # Piecewise linear first crossing of a falling response.
 j=np.where(c<=level)[0]
 if not len(j):return float('nan')
 i=j[0]
 if i==0:return float(t[0])
 if c[i-1]==c[i]:return float(t[i])
 return float(t[i-1]+(level-c[i-1])*(t[i]-t[i-1])/(c[i]-c[i-1]))

def main():
 OUT.mkdir(parents=True,exist_ok=True)
 t,cin=read_curve('fig5_experimental_entrance');u,cout=read_curve('fig5_experimental_exit')
 if not np.array_equal(t,u): raise ValueError('Figure 5 inlet/outlet must use same grid')
 t6,y6=read_curve('fig6_experimental_exit')
 tm,ym=read_curve('fig6_published_fitted_model')
 paper_model_on_markers=np.interp(t6,tm,ym)
 published_author_graphical_rmse=float(np.sqrt(np.mean((paper_model_on_markers-y6)**2)))
 # Fit a gamma residence-time kernel to actual source VECTOR-rendered series.
 gamma=least_squares(lambda v:transfer(t,t,cin,v[0],v[1],v[2])-cout,
      [140.,3.,0.],bounds=([5.,.6,0],[650.,20.,220.]),
      max_nfev=80,xtol=1e-8,ftol=1e-8)
 exp=least_squares(lambda v:transfer(t,t,cin,v[0],1.0,v[1])-cout,
      [135.,0.],bounds=([5.,0.],[650.,220.]),max_nfev=75)
 pred=transfer(t,t,cin,*gamma.x);pred_exp=transfer(t,t,cin,exp.x[0],1.,exp.x[1]);base=cin
 rmse=lambda y:float(np.sqrt(np.mean((y-cout)**2)))
 metrics={'source':'Sales et al. 2005 Fig. 5 & Fig. 6, directly vector-encoded figure symbols in user PDF',
    'figure5':{'n_resampled_graphical_points':int(len(t)),'gamma_mean_s':float(gamma.x[0]),'gamma_shape':float(gamma.x[1]),'additional_delay_s':float(gamma.x[2]),'gamma_rmse_graphical_trace':rmse(pred),
               'exponential_mean_s':float(exp.x[0]),'exponential_delay_s':float(exp.x[1]),'exponential_rmse_graphical_trace':rmse(pred_exp),'no_residence_rmse':rmse(base),
               'input_t50_s':quantile_crossing(t,cin,.5),'outlet_t50_s':quantile_crossing(t,cout,.5)},
    'figure6':{'n_discrete_published_triangle_markers':int(len(t6)),'outlet_t50_s':quantile_crossing(t6,y6,.5),'inlet_time_series_available':False,'authors_published_model_rmse_on_printed_markers':published_author_graphical_rmse,'note':'authors own model-vs-experiment graphical comparison, NOT new independent validation'},
    'qualifications':[
       'Figure 5 data are resampled symbols on a published drawn curve, NOT raw instrument data nor independent i.i.d. data points.',
       'The RMSE describes Figure 5 fitting to resampled graphical representation, not experimental measurement error.',
       'The PDF contains Figure 6 output triangle markers, but no separately measured Figure 6 inlet trace. Thus no independent inlet-to-outlet transfer validation at the 10 L/h condition.',
       'The model is a phenomenological delayed gamma kernel, not the published heterogeneous axial-dispersion adsorption model and not a reactive industrial CFD simulation.',
       'Symbol location/line thickness uncertainty excludes experimental instrument or printing systematics.',
       'The fitted delay may be correlated with kernel shape, so these estimates are not guaranteed individually identifiable.'
    ]}
 # Plot calibrated figure traces and model
 dense=np.linspace(0,1800,500)
 fig,ax=plt.subplots(figsize=(9,5.2))
 ax.plot(t,cin,'.',ms=3,label='Fig. 5 inlet — vector-trace resampling')
 ax.plot(t,cout,'o',ms=3,label='Fig. 5 outlet — vector-trace resampling')
 ax.plot(dense,transfer(dense,t,cin,*gamma.x),'-',lw=2,label='Fitted gamma transfer kernel')
 ax.plot(dense,transfer(dense,t,cin,exp.x[0],1.,exp.x[1]),'--',label='Fitted exponential kernel')
 ax.set(xlabel='Time (s)',ylabel='Normalized tracer concentration',title='Sales et al. (2005), Fig. 5 — source-vector extraction',xlim=(0,1800),ylim=(-.04,1.06));ax.grid(alpha=.22);ax.legend(loc='upper right');fig.tight_layout();fig.savefig(OUT/'figure5_vector_transfer_fit.png',dpi=170);plt.close(fig)
 fig,ax=plt.subplots(figsize=(8.5,4.8));ax.scatter(t6,y6,marker='^',s=24,label='Fig. 6 published experimental triangle symbols');ax.plot(tm,ym,lw=2,label='Original authors plotted fitted model');ax.set(xlim=(0,1800),ylim=(-.04,1.07),xlabel='Time (s)',ylabel='Normalized tracer signal',title='Sales (2005), Fig. 6 — experimental markers, QL=10 L/h');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(OUT/'figure6_vector_markers.png',dpi=170);plt.close(fig)
 with (OUT/'figure5_model_predictions.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['time_s','inlet_digitized','outlet_digitized','gamma_fitted','exponential_fitted','no_residence_baseline']);w.writerows(zip(t,cin,cout,pred,pred_exp,base))
 (OUT/'MODEL_METRICS.json').write_text(json.dumps(metrics,indent=2)+'\n')
 (OUT/'REPORT.md').write_text(f'''# Sales (2005) PDF-vector RTD extraction and transfer fit — v2.6

## Source and provenance

- Sales, Maranhão, Pereira & Abreu (2005), *Experimental Dynamic Evaluation of Three-Phase Reactors*, Brazilian Journal of Chemical Engineering 22(3), pp. 443–452, printed page **450**, Fig. **5** and **6**; DOI 10.1590/S0104-66322005000300012.
- User provided publisher-style **vector PDF**, not an image screenshot. Script `vector_digitize_v26.py` reads the original vector-triangle and circle marker coordinates with PyMuPDF; no OCR or invented observations. Source PDF SHA-256 is in `PROVENANCE.json`.
- Figure 5: **{len(t)} resampled points per displayed series** reconstructed from high-density vector marker geometry. *Do not interpret them as {len(t)} statistically independent laboratory measurements.* Figure 6: **{len(t6)} printed triangle centers**.
- PDF point axes are calibrated using literal plot tick marks: `PROVENANCE.json` stores full geometry and a conservative ±1pt/±2.2pt positioning scenario (not experimental uncertainties). The raw plotted markers may extend slightly beyond [0,1] due to glyph offset; model fitting clips signals to that interval, while CSV retains direct unmodified coordinates.

## Figure 5: QL=5 L/h, QG=300 L/h — input-output model

We use the actually plotted *negative step inlet* together with a non-negative normalized residence-time transfer kernel (gamma distribution). The output is calculated by convolving the **change in inlet** with the delayed RTD kernel. Fit uses SciPy nonlinear least squares.

| Model | Fitted RTD mean (s) | Delay (s) | Graphical curve RMSE |
|---|---:|---:|---:|
| Gamma | {gamma.x[0]:.3f} | {gamma.x[2]:.3f} | {rmse(pred):.6f} |
| Exponential | {exp.x[0]:.3f} | {exp.x[1]:.3f} | {rmse(pred_exp):.6f} |
| No residence | N/A | N/A | {rmse(base):.6f} |

Gamma shape estimate = {gamma.x[1]:.3f}. Figure 5 half-concentration crossing inlet = {metrics['figure5']['input_t50_s']:.1f}s; outlet = {metrics['figure5']['outlet_t50_s']:.1f}s. RMSE is a **graphical-fit diagnostic** not measurement validation; the dense trace points are correlated.

## Figure 6: QL=10 L/h

- **{len(t6)} observed plotted triangles** extracted from the original PDF. Output half-concentration crossing = {metrics['figure6']['outlet_t50_s']:.1f}s.
- The **authors' own fitted-model curve** is separately reconstructed from the printed vector polyline. Its RMSE compared to the plotted triangle centers is **{published_author_graphical_rmse:.5f}** normalized signal units. This is a comparison of the AUTHORS' in-sample fit to their own figure, **not a validation of our gamma kernel**.
- **No numerical inlet waveform** for this independent run is published in Figure 6, so we do **not** assert a valid cross-flow input-output prediction or calculate a misleading independent transfer-model RMSE.
- An experimentally defensible follow-up requires the actual inlet recording or a documented input forcing protocol. A prompt step cannot be equated to the measured inlet without evidence.

## Scientific interpretation and limitations

1. We have **published graph-derived experimental observations**, not true *raw instrument observations*. This is a substantial provenance improvement over v2.5 hand estimates.
2. This validates the extraction method and enables exploratory calibration of a phenomenological transfer kernel at **one flow condition**; it does **not** validate a fully coupled reactive/thermal industrial reactor model.
3. The principal bottleneck for next-stage publication-quality validation is an independent paired inlet/outlet RTD dataset for a second condition, plus known line width/calibration/instrument errors.
4. The figure in Sales (2005) is a step-down tracer response, **not** a pulse-response probability density.

## Reproduce

```bash
python vector_digitize_v26.py --pdf /path/to/download.pdf --out results/v26_vector_extraction
python fit_vector_rtd_v26.py
pytest -q
```

The source PDF is **not bundled** for copyright reasons; it is user-supplied and must be supplied explicitly. No rounding or inferred pseudo-experiments are added to the published figures.
''')
 print(json.dumps(metrics,indent=2))
if __name__=='__main__':main()
