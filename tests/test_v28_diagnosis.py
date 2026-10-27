"""Tests for v2.8 closure and diagnostic — structural, not validation of physics."""
from dataclasses import replace
import numpy as np
import pytest
from heterogeneous_rtd_v27 import Params, build_linear_system, check_global_inventory, tracer_response, inventory_weights
from v28_diagnosis import crossing_time, load_vector_curve


def test_invalid_pore_accessibility_rejected():
    with pytest.raises(AssertionError):build_linear_system(replace(Params(n_axial=8),accessible_pore_fraction=0.0))
    with pytest.raises(AssertionError):build_linear_system(replace(Params(n_axial=8),accessible_pore_fraction=1.1))


def test_pore_capacity_scales_with_accessibility_and_mobile_capacity_fixed():
    a=Params(n_axial=8,accessible_pore_fraction=1.0)
    b=replace(a,accessible_pore_fraction=.25)
    w1=inventory_weights(a);w2=inventory_weights(b)
    step=1+2*a.n_radial
    np.testing.assert_allclose(w1[::step],w2[::step],rtol=0,atol=0)
    np.testing.assert_allclose(w2[1:step],.25*w1[1:step],rtol=1e-12)


def test_new_closure_keeps_global_mass_conservation():
    for x in [.05,.25,1.0]:
        p=Params(n_axial=9,accessible_pore_fraction=x)
        r,inlet=check_global_inventory(p)
        assert r<1e-11 and inlet<1e-11


def test_reducing_accessibility_approaches_axial_only():
    p=Params(n_axial=9,n_radial=2)
    t=np.linspace(0,600,31)
    pulse=(t<80).astype(float)
    y0=tracer_response(t,t,pulse,p,enable_pore=False)
    y1=tracer_response(t,t,pulse,replace(p,accessible_pore_fraction=.001))
    assert np.max(np.abs(y0-y1))<.003


def test_crossing_time_monotonic_finite():
    t=np.array([0.,10.,20.]);v=np.array([1.,.6,.2])
    assert np.isclose(crossing_time(t,v),12.5)
    assert np.isnan(crossing_time(t,np.ones(3)))


def test_published_figure5_inlet_and_outlet_have_time_alignment():
    x,cin=load_vector_curve('fig5_experimental_entrance')
    y,cout=load_vector_curve('fig5_experimental_exit')
    assert np.array_equal(x,y)
    assert len(x)>40 and np.all(np.diff(x)>0)
    assert 0<=cin.min()<=1 and 0<=cout.min()<=1
