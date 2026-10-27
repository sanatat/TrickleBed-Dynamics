from pathlib import Path
from dataclasses import replace
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v5_model import Params,simulate,unpack,volumes

out=Path('results/v5');out.mkdir(parents=True,exist_ok=True)
rows=[]; curves={}
for label,wet in [('initially_dry',.2),('prewetted',.88)]:
 for chi in (1,10):
    p=Params(wet_init=wet,chi=chi)
    t,Y=simulate(p)
    vals=[]
    for j,tt in enumerate(t):
        N,T,s=unpack(Y[:,j],p);C=N/volumes(s,p)[:,None,:]
        vals.append([tt, C[0,0,-1],T.max(),s.mean(),C[2,0,-1],C[3,1,-1]])
    vals=np.array(vals); curves[(label,chi)]=vals
    with (out/f'{label}_chi{chi}.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['time_s','outlet_liquid_A_mol_m3','max_temperature_K','mean_pore_liquid_saturation','outlet_pore_A_mol_m3','outlet_pore_gas_H_mol_m3']);w.writerows(vals)
    rows.append([label,chi,float(vals[-1,1]),float(vals[:,2].max()),float(vals[-1,3])])
with (out/'summary.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['history','chi','final_outlet_liquid_A_mol_m3','peak_temperature_K','final_pore_saturation']);w.writerows(rows)
for index,name,ylab in [(1,'outlet_octene','Liquid outlet isooctene (mol/m3)'),(2,'peak_temperature','Max bed temperature (K)'),(3,'pore_saturation','Pore liquid saturation')]:
 plt.figure(figsize=(8,4.5))
 for (label,chi),v in curves.items():plt.plot(v[:,0],v[:,index],label=f'{label}, chi={chi}')
 plt.axvline(200,ls=':',color='gray');plt.xlabel('Time (s)');plt.ylabel(ylab);plt.legend(fontsize=8);plt.tight_layout();plt.savefig(out/f'{name}.png',dpi=150);plt.close()
print('V5 SYNTHETIC RESULTS')
for r in rows:print(r)
