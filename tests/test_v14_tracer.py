import numpy as np
from v6_model import Params
from tracer_coupled_v14 import TracerConfig,tracer_rhs,wetting_schedule,run_tracer

def test_flux_conservation_closed():
    p=Params(n=6,t_step=1000,wet_init=.4)
    fun=wetting_schedule(p,40)
    y=np.r_[np.linspace(.1,.7,p.n),np.linspace(.3,.9,p.n)]
    dy=tracer_rhs(15,y,p,TracerConfig(),fun,inlet_on=False)
    # Sum changes = advective + diffusive boundary flux only; diffusion end boundaries zero
    outflux=p.u_L*y[p.n-1]/p.V_L
    expected=-outflux/(p.length/p.n)
    assert np.isclose(np.sum(dy), expected,rtol=1e-10,atol=1e-10)

def test_no_tracer_when_inlet_disabled():
    p=Params(n=5,t_step=1000)
    r=run_tracer(p,TracerConfig(),np.linspace(0,100,8),inlet_on=False)
    assert np.allclose(r['outlet_conc'],0)

def test_positive_inventories():
    p=Params(n=8,t_step=1000)
    r=run_tracer(p,TracerConfig(),np.linspace(0,160,20))
    assert np.min(r['outlet_conc'])>=-1e-10
