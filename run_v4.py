from pathlib import Path
from dataclasses import replace
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v4_model import Params,simulate,unpack

out=Path('results/v4');out.mkdir(parents=True,exist_ok=True)
series=[];rows=[]
for wetlabel,wi in [('initially_dry',.20),('prewetted',.88)]:
  for chi in (1.,10.):
    p=replace(Params(),wet_init=wi,chi=chi)
    t,y=simulate(p)
    a=y.reshape(14,p.n,-1)
    out_oct=a[0,-1];max_T=a[12].max(axis=0);mean_w=a[13].mean(axis=0)
    path=out/f'{wetlabel}_chi{int(chi)}.csv'
    with path.open('w',newline='') as f:
      writer=csv.writer(f);writer.writerow(['time_s','outlet_octene_mol_m3_liquid','max_temperature_K','mean_wetting'])
      writer.writerows(zip(t,out_oct,max_T,mean_w))
    rows.append([wetlabel,chi,round(float(1-out_oct[-1]/p.li_A),5),round(float(max_T.max()),4),round(float(mean_w[-1]),5)])
    series.append((f'{wetlabel}, chi={chi:g}',t,out_oct,max_T,mean_w))
with (out/'summary.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['scenario','chi','final_outlet_concentration_proxy','max_temperature_K','final_mean_wet']);w.writerows(rows)
for fname,title,axis,value in [('max_temperature.png','Peak temperature','K',3),('outlet_octene.png','Outlet isooctene','mol/m³',2),('mean_wetting.png','Mean wetting fraction','fraction',4)]:
  fig,ax=plt.subplots(figsize=(7,4.2))
  for entry in series:
    ax.plot(entry[1],entry[value],label=entry[0])
  ax.axvline(200,linestyle='--',alpha=.5)
  ax.set(xlabel='Time (s)',ylabel=axis,title=title);ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/fname,dpi=160);plt.close(fig)
print('Scenario, chi, concentration proxy, peak K, wet final')
for r in rows: print(*r,sep=', ')
