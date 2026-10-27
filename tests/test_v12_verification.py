from verification_v12 import integrate,budget_derivative,energy_residual
from v6_model import Params,initial
import numpy as np

def test_timestep_event_present():
    p=Params(n=3,t_end=20.,t_step=10.)
    t,Y=integrate(p,1e-6,1e-8,21)
    assert np.any(np.isclose(t,10.)) and np.all(np.diff(t)>0)

def test_equivalents_match_boundaries():
    p=Params(n=4)
    actual,expected=budget_derivative(5.,initial(p),p)
    assert np.allclose(actual,expected,rtol=1e-9,atol=1e-10)

def test_energy_rhs_identity():
    p=Params(n=4)
    assert energy_residual(5.,initial(p),p)<1e-11
