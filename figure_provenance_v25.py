'''Audit what is supported by the published Sales 2005 RTD figures.'''
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/v25_provenance'

def run():
    OUT.mkdir(parents=True,exist_ok=True)
    prev=np.genfromtxt(ROOT/'results/v20_sales_digitization/figure5_approx_experimental.csv',delimiter=',',skip_header=1)
    t,entrance,exit=prev.T
    fig,ax=plt.subplots(figsize=(9,5))
    ax.plot(t,entrance,'o-',label='Fig. 5 inlet: inherited manual estimates')
    ax.plot(t,exit,'s-',label='Fig. 5 outlet: inherited manual estimates')
    ax.set(xlabel='Time (s)',ylabel='Normalized mass fraction',title='Sales (2005), Figure 5 — legacy approximate readings')
    ax.set(xlim=(0,1800),ylim=(-.05,1.07));ax.legend();ax.grid(alpha=.3);fig.tight_layout();fig.savefig(OUT/'legacy_approx_not_pixel_digitized.png',dpi=160);plt.close(fig)
    checks={'figure5_inlet_fall':bool(entrance[0]>entrance[-1]),'figure5_outlet_fall':bool(exit[0]>exit[-1]),'figure5_outlet_delayed':bool(np.interp(450,t,exit)>np.interp(450,t,entrance)),'figure5_plot_q_liquid_L_h':5,'figure6_plot_q_liquid_L_h':10,'both_q_gas_L_h':300,'measured_fig6_inlet_available':False,'source_pdf_local_available':False,'new_pixel_digitization_completed':False,'raw_experimental_timeseries_available':False,'previous_data_points':len(t)}
    (OUT/'audit.json').write_text(json.dumps(checks,indent=2)+'\n')
    (OUT/'REPORT.md').write_text('''# v2.5 — Source-verification and reproducible digitization workflow

Published paper: Sales et al. (2005), https://doi.org/10.1590/S0104-66322005000300012, Figure 5, Figure 6, Table 3, printed page 450.

## Verified from the actual journal PDF
- Figure 5 gives **experimental entrance**, **experimental exit**, and **theoretical exit**, at QL=5 L/h, QG=300 L/h.
- Figure 6 gives **experimental outlet** and **modelled outlet**, at QL=10 L/h, QG=300 L/h. A separate measured inlet series is not printed.
- Both are step-down-like decaying **concentration** curves, not plotted pulse RTD PDFs. Do not directly normalize the falling curves as residence-time probability density without calculating a response derivative with appropriate noise propagation.
- Table 3 reports six fitted hydrodynamic parameter sets (Pe, holdup, wetting); these are derived parameters, not raw tracer time series.

## What v2.5 does
- Inspects the legacy 13 approximate Figure 5 samples and audits sign/delay.
- Provides `pixel_digitize_v25.py` to manually record pixel coordinates and calibrate measured axes (x: 0–1800 s; y: 0–1), with explicit provenance and user-defined pixel-error propagation.
- Adds regression tests checking calibration and plot protocols.

## What is NOT done
- The journal PDF can be viewed in a browser, but downloading original image bytes into the local code runtime was unsuccessful. No new **calibrated pixels** have been extracted here. The legacy 13-point and Figure 6 manual estimates remain *approximate, non-validated readings*.
- No experimental inlet trace for Figure 6, direct independent experimental RTD validation, or CFD verification is claimed.

## How to finish a fully reproducible extraction on the user’s machine
1. Download the publisher PDF from https://www.scielo.br/j/bjce/a/fN94qbKQVRVXmTqN8R8BPxP/?format=pdf&lang=en and render printed page 450 (PDF page 8) at 300–600 dpi.
2. Crop Figure 5 and Figure 6 separately without distorting the axis aspect ratios.
3. Mark axis pixel references at 0 s, 1800 s, mass fraction 0 and 1.0. The utility expects X0 X1800 Y0 Y1 (pixel y counts downward).
4. Click each observed marker into its own CSV labeled `experimental_inlet` or `experimental_exit` for Figure 5; for Figure 6 use `experimental_exit`. Keep theoretical curves separate.
5. Run e.g. `python pixel_digitize_v25.py --pixels figure5_pixels.csv --axis 120 920 830 210 --output figure5_calibrated.csv --pixel-error 2` (**these axis coordinates are examples only**).
6. Repeat clicking three times or with a second operator; use differences to quantify extraction variability. Do not confuse plotting extraction uncertainty with measurement uncertainty.
7. Freeze calibration on Figure 5 before predicting Figure 6; without Fig.6 measured inlet, compare only under explicitly assumed boundary conditions.

## Scientific decision
The next meaningful improvement is accurate source-image extraction, not another optimised fit to the manually estimated curve. Give a conservative status in CV: paper-grounded RTD analysis with published parameter audit, exploratory inlet–outlet modeling, and documented data limitations.
''')
    return checks
if __name__=='__main__': print(json.dumps(run(),indent=2))
