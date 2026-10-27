import numpy as np
from dataclasses import replace
from v6_model import Params
from reactive_rtd_v15 import Coupling,diffusive_divergence,integrate,full_rhs

def test_diffusion_conservative():
    x=np.array([[1.,3.,2.,4.]])
    assert abs(diffusive_divergence(x,.005,.2).sum())<1e-12

def test_diffusion_constant_zero():
    assert np.allclose(diffusive_divergence(np.ones((3,6)),2e-5,.1),0)

def test_tracer_off_no_source():
    p=Params(n=4,t_step=1000,t_end=10); r=integrate(p,Coupling(),np.linspace(0,10,4),tracer_on=False)
    assert np.max(np.abs(r['tracer']))<1e-12

def test_passive_tracer_no_feedback():
    p=Params(n=4,t_step=1000,t_end=30)
    c=Coupling(dispersion=1e-5)
    a=integrate(p,c,np.linspace(0,30,5),True);b=integrate(p,c,np.linspace(0,30,5),False)
    assert np.allclose(a['octene'],b['octene'],rtol=1e-8,atol=1e-8)
    assert np.allclose(a['peak_T'],b['peak_T'],rtol=1e-8,atol=1e-8)

def test_dispersion_changes_species():
    p=Params(n=6,t_step=1000,t_end=70)
    t=np.linspace(0,70,7)
    a=integrate(p,Coupling(dispersion=0),t);b=integrate(p,Coupling(dispersion=4e-4),t)
    assert abs(a['octene'][-1]-b['octene'][-1])>1e-5
