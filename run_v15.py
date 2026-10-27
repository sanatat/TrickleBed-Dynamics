from dataclasses import replace
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v6_model import Params
from reactive_rtd_v15 import Coupling,integrate

def main():
    out=Path('results/v15_joint');out.mkdir(parents=True,exist_ok=True)
    p=Params(n=12,t_end=360,t_step=180,wet_init=.35,wet_adv=.8,wet_rec=.46,cap_tau=110,u_L=.0015)
    t=np.linspace(0,360,91)
    configs={'zero_dispersion':Coupling(dispersion=0.),'low_dispersion':Coupling(dispersion=5.43e-6), 'high_dispersion':Coupling(dispersion=1.4e-5)}
    result={label:integrate(p,c,t) for label,c in configs.items()}
    metrics={label:dict(final_outlet_octene=float(r['octene'][-1]),peak_temperature_K=float(max(r['peak_T'])),max_outlet_tracer=float(max(r['tracer']))) for label,r in result.items()}
    with (out/'trajectories.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s']+[f'{label}_{key}' for label in result for key in ['octene','peak_T','tracer']]);
        for i in range(len(t)):w.writerow([t[i]]+[float(result[label][key][i]) for label in result for key in ['octene','peak_T','tracer']])
    for key,ylabel,filename in [('octene','Outlet octene (mol/m3)','outlet_octene.png'),('peak_T','Maximum bed temperature (K)','peak_temperature.png'),('tracer','Outlet tracer (arbitrary normalized units)','rtd_tracer.png')]:
        fig,ax=plt.subplots(figsize=(8,4.5))
        for label,r in result.items():ax.plot(r['time'],r[key],label=label.replace('_',' '))
        ax.axvline(p.t_step,ls=':',c='gray',label='Flow step');ax.set(xlabel='Time (s)',ylabel=ylabel);ax.legend();fig.tight_layout();fig.savefig(out/filename,dpi=140);plt.close(fig)
    metadata={'dataset':'synthetic, NO external validation','study':'Concurrent integration of reactive 4-compartment balances and tracer with shared mobile-liquid axial dispersion','parameter_provenance':'all dispersions are illustration; low value inspired by v1.4 synthetic fit; high value synthetic truth; neither empirically verified','metrics':metrics,'limitations':['No experimental calibration','Pore-exchange coefficient for tracer independent of chemistry interface coefficients','Tracer passive; chemical reaction affects neither wettability nor hydraulics','Pressure, phase change and actual pore geometry not solved','Concentration reduction does not equal overall conversion']}
    (out/'metrics.json').write_text(json.dumps(metadata,indent=2))
    (out/'REPORT.md').write_text('# Joint reactive–RTD simulation (v1.5)\n\nNumerically integrate a passive tracer simultaneously with the v6 four-compartment reactive model. Axial dispersion now affects the mobile-liquid reactant/product fields as well as the tracer, providing a shared hydrodynamic transport hypothesis. A step change in liquid flow at 180 s changes residence time and prescribed wetting.\n\n**All inputs and comparison curves are synthetic.** This is *numerical coupling*, not experimental validation or a fully consistent multifluid hydrodynamic closure. The tracer cannot physically change chemical kinetics (it is passive); the shared transport operator links the responses. The calibrated dispersion value from v1.4 is weakly identifiable and not validated. See metrics.json and code.\n')
    print(json.dumps(metrics,indent=2))
if __name__=='__main__':main()
