import numpy as np
from optimal_design_v17 import information, score, sensitivities, BASELINE

def test_sensitivity_shape_and_finiteness():
    j=sensitivities(1.0,t=np.linspace(0,90,12))
    assert j.shape==(12,2)
    assert np.isfinite(j).all()

def test_information_psd():
    j=sensitivities(1.0,t=np.linspace(0,90,12))
    info=information((1.0,),{1.0:j})
    assert np.allclose(info,info.T)
    assert np.linalg.eigvalsh(info).min()>=-1e-8

def test_score_two_identical_flow_observations():
    j=sensitivities(1.0,t=np.linspace(0,90,12))
    val=score((1.0,1.0),{1.0:j})
    assert np.isfinite(val)
