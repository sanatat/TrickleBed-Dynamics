import numpy as np
import pytest
from pixel_digitize_v25 import calibrate
from figure_provenance_v25 import run

def test_axes_endpoints():
    v,u=calibrate([[100,500],[900,100]],(100,900,500,100))
    assert np.allclose(v,[[0,0],[1800,1]])

def test_pixel_uncertainty():
    _,u=calibrate([[500,300]],(100,900,500,100),pixel_error=2)
    assert u==(4.5,0.005)

def test_bad_bounds():
    with pytest.raises(ValueError):calibrate([[0,0]],(5,5,200,100))

def test_publication_audit():
    result=run()
    assert result['figure5_inlet_fall'] and result['figure5_outlet_delayed']
    assert not result['new_pixel_digitization_completed']
