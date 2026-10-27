import numpy as np
from dataclasses import replace
from v5_model import Params,initial,unpack,volumes,rhs,invariants,simulate

def test_volumes_change_and_total_constant():
    p=Params(n=3); s=np.array([.2,.5,.88]);v=volumes(s,p)
    assert np.allclose(v[2]+v[3],p.pore_vol)
    assert v[2,0]<v[2,2] and v[3,0]>v[3,2]

def test_closed_invariants_capillary_and_reaction():
    p=Params(n=4);y=initial(p); dy=rhs(0,y,p,flow=False)
    c,h=invariants(dy,p)
    assert abs(c)<1e-9 and abs(h)<1e-9

def test_closed_each_species_without_reaction():
    p=Params(n=3);dy=rhs(0,initial(p),p,flow=False,reaction=False)
    n,_,_=unpack(dy,p)
    assert np.allclose(n.sum(axis=(0,2)),0,atol=1e-10)

def test_bounded_saturation_and_nonnegative_inventories():
    p=Params(n=5,t_end=300)
    t,y=simulate(p); n,_,s=unpack(y[:,-1],p)
    assert s.min()>=p.s_min and s.max()<=p.s_max
    assert n.min()>-1e-6

def test_step_changes_saturation_direction():
    p=Params(n=4,t_end=300)
    t,y=simulate(p); before=unpack(y[:,np.searchsorted(t,200)],p)[2].mean(); end=unpack(y[:,-1],p)[2].mean()
    assert end<before
