import numpy as np
from dataclasses import replace
from v6_model import Params,initial,unpack,volumes
from quality_audit_v11 import residual_inventory,metrics

def test_closed_sources_conserve_both_equivalents():
 p=Params(n=4); y=initial(p).reshape(14,p.n);y[:12]*=np.random.default_rng(123).uniform(.6,1.4,(12,p.n))
 assert max(residual_inventory(y.ravel(),p,reaction=True))<1e-9

def test_closed_transfer_only_conserves():
 p=Params(n=5);y=initial(p)
 assert max(residual_inventory(y,p,reaction=False))<1e-9

def test_pore_volumes_sum_to_constant():
 p=Params(n=5);s=np.linspace(.05,.95,p.n);V=volumes(s,p)
 assert np.allclose(V[2]+V[3],p.pore_vol)

def test_varying_assumed_rate_changes_output():
 p=Params(n=6,t_end=150.,t_step=80.)
 a=metrics(replace(p,k_ref=p.k_ref/2))
 b=metrics(replace(p,k_ref=p.k_ref*2))
 assert a['outlet_alkene_mol_m3']!=b['outlet_alkene_mol_m3']

def test_short_simulation_nonnegative():
 m=metrics(Params(n=6,t_end=150.,t_step=80.))
 assert m['min_species_inventory_mol_m3bed']>-1e-7
 assert 0<m['saturation_min']<=m['saturation_max']<1
