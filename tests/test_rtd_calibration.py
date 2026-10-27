import numpy as np
import pytest
from rtd_calibration import RTDParams, impulse_response, make_synthetic_data, fit_rtd, residual_bootstrap, rmse

def test_normalized():
    t=np.linspace(0,1000,20001)
    e=impulse_response(t,RTDParams(85,5.2,7))
    assert abs(np.trapezoid(e,t)-1)<1e-6

def test_nonnegative():
    assert np.min(impulse_response(np.arange(300),RTDParams(85,5.2,7)))>=0

def test_before_delay_zero():
    assert np.all(impulse_response(np.arange(7),RTDParams(85,5.2,7))==0)

def test_residence_scales_with_flow():
    t=np.linspace(0,1000,20001)
    for q in [1.,1.35]:
        y=impulse_response(t,RTDParams(85,5.2,7),q)
        mean=np.trapezoid(t*y,t)
        assert abs(mean-(7+85/q)) <.01

def test_parameter_recovery_noise_free():
    t=np.arange(0,361,3.)
    y=impulse_response(t,RTDParams(85,5.2,7))
    p,_=fit_rtd(t,y)
    assert abs(p.mean_residence_s-85)<.001
    assert abs(p.effective_tanks-5.2)<.001
    assert abs(p.delay_s-7)<.001

def test_holdout_beats_single_cstr():
    _, d, _=make_synthetic_data()
    t=d['calibration']['time_s']
    p,_=fit_rtd(t,d['calibration']['observed'])
    baseline,_=fit_rtd(t,d['calibration']['observed'],tanks_fixed=1.0)
    held=d['holdout']
    assert rmse(held['observed'],impulse_response(t,p,1.35)) < rmse(held['observed'],impulse_response(t,baseline,1.35))

def test_bootstrap_finite():
    _,d,_=make_synthetic_data()
    t=d['calibration']['time_s']
    p,_=fit_rtd(t,d['calibration']['observed'])
    samples=residual_bootstrap(t,d['calibration']['observed'],p,runs=10)
    assert samples.shape==(10,3)
    assert np.all(np.isfinite(samples))

def test_invalid_parameters():
    with pytest.raises(ValueError):
        impulse_response([1,2],RTDParams(0,5,0))
