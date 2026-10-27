import numpy as np
from dataclasses import replace
from v6_model import Params, sherwood, interface_areas, simulate, unpack, invariants, rhs, initial

def test_sherwood_stagnant():
    assert np.isclose(sherwood(0,10),2)

def test_sherwood_increases():
    assert sherwood(10,100)>sherwood(1,100)>2

def test_sherwood_invalid():
    import pytest
    with pytest.raises(ValueError): sherwood(-1,1)

def test_area_endpoints():
    p=Params()
    a=interface_areas(np.array([0.,1.]),p)
    assert a['GL'][0]==0 and a['GD'][1]==0
    assert all(np.all(v>=0) for v in a.values())

def test_area_partition():
    p=Params()
    s=np.linspace(0,1,11)
    a=interface_areas(s,p)
    base=6*(1-p.void)/p.pellet_diameter
    assert np.allclose(a['GL']/base+a['GD']/(p.pellet_porosity*base),1)

def test_invariant_closed():
    p=Params(n=3,t_end=10,t_step=20)
    t,Y=simulate(p,flow=False)
    q0=invariants(Y[:,0],p); q1=invariants(Y[:,-1],p)
    assert all(abs(a-b)<1e-5 for a,b in zip(q0,q1))

def test_finite():
    p=Params(n=3,t_end=20,t_step=10)
    _,Y=simulate(p)
    assert np.isfinite(Y).all()

def test_changes_when_area_disabled():
    p=Params(n=3,t_end=20,t_step=10)
    _,a=simulate(p)
    _,b=simulate(replace(p,interface_area_mode=False))
    assert not np.allclose(a,b)

def test_initial_rhs_finite():
    p=Params(n=4)
    assert np.isfinite(rhs(0,initial(p),p)).all()
