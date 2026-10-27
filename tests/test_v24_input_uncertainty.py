import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from input_uncertainty_v24 import predict_step,run

def test_step_limit():
    t=np.array([0.,10.,100.]); y=predict_step(t,2,20)
    assert np.all((y>=0)&(y<=1)) and np.all(np.diff(y)<=0)

def test_ramp_monotone():
    t=np.linspace(0,500,100); y=predict_step(t,3,100,30,120)
    assert np.all(np.diff(y)<=1e-12)

def test_invalid():
    import pytest
    with pytest.raises(ValueError): predict_step([0,1],0,10)

def test_run_outputs():
    m=run(); assert m['n_points']>3 and not m['raw_experiment_available']
