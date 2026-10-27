"""Structural verification only — cannot substitute for experimental validation."""
from dataclasses import replace
import numpy as np
from heterogeneous_rtd_v27 import (Params, build_linear_system, inventory_weights,
                                    check_global_inventory, tracer_response,
                                    effective_dispersion_per_bed)

def small():
    return Params(n_axial=9,n_radial=3)

def test_reported_Pe_maps_to_positive_dispersion():
    p=small()
    assert np.isclose(effective_dispersion_per_bed(p),
        .8*.001/(4.84*.4))

def test_uniform_equilibrium_stationary():
    p=small()
    A,b=build_linear_system(p)
    eq=np.zeros(A.shape[0]);stride=1+2*p.n_radial
    for j in range(p.n_axial):
        eq[j*stride]=1
        eq[j*stride+1:j*stride+1+p.n_radial]=1
        eq[j*stride+1+p.n_radial:(j+1)*stride]=p.K_ads_m3_kg
    assert np.max(np.abs(A@eq+b))<1e-11

def test_global_discrete_tracer_inventory():
    p=small()
    r,f=check_global_inventory(p)
    assert r<1e-11 and f<1e-11
    w=inventory_weights(p)
    assert np.all(w>0)

def test_inlet_outlet_identity_linear_model():
    p=small()
    time=np.linspace(0,1200,35)
    incoming=np.where(time>=150,0.,1.)
    out=tracer_response(time,time,incoming,p)
    assert 0<=out.min() and out.max()<=1.00001
    assert out[0]>0.999
    assert out[-1]<0.30

def test_zero_surface_exchange_matches_axial_model():
    p=small()
    time=np.linspace(0,900,22);incoming=np.where(time>150.,0.,1.)
    no_pore=tracer_response(time,time,incoming,p,enable_pore=False)
    no_transfer=tracer_response(time,time,incoming,replace(p,interface_scale=0.),enable_pore=True)
    assert np.max(abs(no_pore-no_transfer))<6e-4

def test_adsorption_capacity_changes_when_K_changes():
    p=small()
    A,_=build_linear_system(p)
    A0,_=build_linear_system(replace(p,K_ads_m3_kg=0))
    assert (A-A0).nnz > 0

def test_discrete_conservation_for_random_concentrations():
    p=small()
    A,b=build_linear_system(p)
    w=inventory_weights(p)
    y=np.random.default_rng(15).random(A.shape[0])
    inlet=.73
    flux_balance=float(w@(A@y + b*inlet))
    fluid_out=y[(p.n_axial-1)*(1+2*p.n_radial)]
    expected=p.superficial_u_m_s*(inlet-fluid_out)
    assert abs(flux_balance-expected)<1e-12
