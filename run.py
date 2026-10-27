from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import Config, simulate, metrics, unpack

c=Config()
output=Path(__file__).parent/'results';output.mkdir(exist_ok=True)
series={}
for name,wet0 in [('initially_dry',0.20),('prewetted',0.88)]:
    t,y=simulate(c,wet0)
    m=metrics(t,y,c);series[name]=(t,y,m)
    with (output/f'{name}.csv').open('w',newline='') as fh:
        writer=csv.writer(fh);ks=list(m)
        writer.writerow(ks)
        writer.writerows(zip(*[m[k] for k in ks]))
    print(f'{name}: final X={m["conversion_out"][-1]:.4f}, max T={m["temperature_peak_K"].max():.3f} K, final wet={m["wetting_mean"][-1]:.3f}')

for metric,title,file in [
  ('conversion_out','Outlet conversion proxy: 1 - Aout/Ain','outlet_conversion.png'),
  ('temperature_peak_K','Maximum axial temperature (K)','peak_temperature.png'),
  ('wetting_mean','Mean surface wetting fraction','wetting_history.png')]:
    fig,ax=plt.subplots(figsize=(7.5,4.6))
    for name,(t,y,m) in series.items(): ax.plot(t,m[metric],label=name.replace('_',' '))
    ax.axvline(c.step_time,ls='--',color='gray',label='flow step')
    ax.set(xlabel='Time (s)',ylabel=title)
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(output/file,dpi=160);plt.close(fig)
fig,ax=plt.subplots(figsize=(7.5,4.6))
for name,(t,y,m) in series.items():
    A,H,P,Ap,Hp,Pp,T,f=[z.reshape(c.n,-1) for z in unpack(y,c.n)]
    ax.plot((np.arange(c.n)+.5)*c.length/c.n,T[:,-1],label=name.replace('_',' '))
ax.set(xlabel='Axial position (m)',ylabel='Temperature at final time (K)')
ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(output/'axial_final_temperature.png',dpi=160);plt.close(fig)

with (output/'SUMMARY.txt').open('w') as fh:
    fh.write('SYNTHETIC EDUCATIONAL SIMULATION, NOT VALIDATED OR A REPRODUCTION OF 2012 RESULTS.\n')
    for name,(t,y,m) in series.items():
        fh.write(f'{name}: final conversion proxy {m["conversion_out"][-1]:.5f}; max temperature {m["temperature_peak_K"].max():.2f} K; final wetting {m["wetting_mean"][-1]:.4f}\n')
