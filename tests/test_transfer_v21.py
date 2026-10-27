import numpy as np
from sales2005_transfer_v21 import transfer,gamma_cdf
def test_constant_input():
 t=np.linspace(0,500,30); y=np.ones_like(t); assert np.allclose(transfer(t,t,y,100,2,20),1)
def test_step_down_monotone_and_bounded():
 t=np.linspace(0,500,101); y=np.where(t<100,1.,0.); p=transfer(t,t,y,60,3,10);assert np.all(np.diff(p)<1e-9);assert np.min(p)>=-1e-8;assert np.max(p)<=1+1e-8
def test_cdf_bounds():
 assert np.isclose(gamma_cdf(0,100,2),0);assert gamma_cdf(1000,100,2)>.99
