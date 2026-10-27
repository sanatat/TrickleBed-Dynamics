import numpy as np
from validation_v22 import rmse,predict,fit

def test_rmse_identity():
    assert rmse([0,1],[0,1])==0

def test_constant_input_stays_constant():
    t=np.arange(0,1200,100.0); y=predict(t,np.ones_like(t),160,3,0)
    assert np.max(np.abs(y-1)) < 1e-10

def test_step_down_is_bounded():
    t=np.linspace(0,1000,21);cin=(t<200).astype(float)
    y=predict(t,cin,150,4,0)
    assert np.all(y>=-1e-7) and np.all(y<=1+1e-7)

def test_fit_returns_finite():
    t=np.linspace(0,900,12); cin=np.exp(-np.maximum(t-100,0)/140);cout=predict(t,cin,130,3,0)
    p,res=fit(t,cin,cout,t<=600,'gamma');assert np.all(np.isfinite(p)) and res.success
