"""Qualitative checks against 2012 paper: chi=1 vs chi=10, NOT digitized reproduction."""
from pathlib import Path
from dataclasses import replace
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from paper_extension import PaperConfig, simulate, metrics
out=Path(__file__).parent/'results'/'paper_extension';out.mkdir(parents=True,exist_ok=True)
rows=[]
for chi in (1.0,10.0):
    for wet0 in (0.2,0.88):
        c=PaperConfig(chi=chi,wet_init=wet0)
        t,Y=simulate(c);m=metrics(t,Y,c)
        A,H,P,Ap,Hp,Pp,T,f=Y.reshape(8,c.n,-1)
        name=f'chi{int(chi)}_wet{str(wet0).replace(".","p")}'
        with (out/f'{name}.csv').open('w',newline='') as fp:
            writer=csv.writer(fp);writer.writerow(m.keys());writer.writerows(zip(*m.values()))
        rows.append((name,m['conversion_out'][-1],m['T_peak_K'].max(),m['wetting_mean'][-1]))
        for category,data in [('temperature',T[:,-1]),('product',P[:,-1])]:
            plt.figure(category,figsize=(7,4.5));z=(np.arange(c.n)+.5)*c.length/c.n
            plt.plot(z,data,label=f'chi={chi:g}; wet0={wet0:g}')
for category in ('temperature','product'):
    plt.figure(category);plt.xlabel('Axial position (m)');plt.ylabel('Temperature (K)' if category=='temperature' else 'Bulk product (mol/m3)');plt.title('Illustrative reduced model; not Gorshkova Fig. 2')
    plt.legend(fontsize=8);plt.grid(alpha=.25);plt.tight_layout();plt.savefig(out/f'axial_{category}.png',dpi=150);plt.close()
with (out/'summary.csv').open('w',newline='') as fp:
    writer=csv.writer(fp);writer.writerow(('case','outlet_conversion_proxy','max_temperature_K','final_wet_fraction'));writer.writerows(rows)
for r in rows: print(f'{r[0]}: outlet X={r[1]:.5f}; Tpeak={r[2]:.3f} K; wet={r[3]:.4f}')
