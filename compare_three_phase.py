from dataclasses import replace
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from three_phase import Config,simulate

out=Path(__file__).resolve().parent/'results'/'three_phase';out.mkdir(parents=True,exist_ok=True)
rows=[]
for chi in (1.,10.):
 for wet in (.2,.88):
  c=Config(chi=chi,initial_wet=wet)
  t,Y=simulate(c)
  X=Y.reshape(11,c.n,-1)
  L,G,W,T,f=X[:3],X[3:6],X[6:9],X[9],X[10]
  label=f'chi{int(chi)}_wet{str(wet).replace(".","p")}'
  conversion=1-L[0,-1,-1]/c.inlet_A_liquid
  peak=float(np.max(T))
  rows.append((label,float(conversion),peak,float(f[:,-1].mean()),float(G[1,-1,-1])))
  with open(out/f'{label}.csv','w',newline='') as file:
   w=csv.writer(file);w.writerow(['time_s','outlet_octene_mol_m3liquid','outlet_gas_H2_mol_m3gas','peak_T_K','mean_wet'])
   w.writerows(zip(t,L[0,-1],G[1,-1],T.max(axis=0),f.mean(axis=0)))
  for name,y in [('peak_temperature',T.max(axis=0)),('outlet_octene',L[0,-1]),('wetting',f.mean(axis=0))]:
   plt.figure(name,figsize=(7,4));plt.plot(t,y,label=label)
  z=(np.arange(c.n)+.5)*c.length/c.n
  plt.figure('axial_octene',figsize=(7,4));plt.plot(z,L[0,:,-1],label=label)
for name,ylabel in [('peak_temperature','Peak temperature (K)'),('outlet_octene','Liquid outlet octene (mol/m3-liquid)'),('wetting','Mean dynamic wetting fraction'),('axial_octene','Final liquid octene (mol/m3-liquid)')]:
 plt.figure(name);plt.xlabel('Axial position (m)' if name=='axial_octene' else 'Time (s)');plt.ylabel(ylabel)
 plt.title('Synthetic 1D explicit-gas model; NOT 2012 CFD reconstruction');plt.legend(fontsize=8);plt.grid(alpha=.25);plt.tight_layout();plt.savefig(out/f'{name}.png',dpi=150);plt.close()
with open(out/'summary.csv','w',newline='') as file:
 w=csv.writer(file);w.writerow(['case','outlet_conversion_proxy','peak_T_K','wet_final','outlet_gas_H2']);w.writerows(rows)
for row in rows:print(row)
