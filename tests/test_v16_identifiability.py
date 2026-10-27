import numpy as np
from identifiability_v16 import setup, observed, estimate
from dataclasses import replace

def test_identifiability_noise_free():
 p=setup(6); t=np.linspace(0,220,40); theta=(.0025,1.4e-5); y=observed(p,t,theta)
 fit,diag=estimate([(p,t,y,1e-4)],guess=(.0025,1.4e-5))
 assert fit.success
 assert np.allclose(diag['estimate'],theta,rtol=1e-4)

def test_tracer_nonnegative():
 p=setup(6);t=np.linspace(0,180,35); y=observed(p,t,(.0025,1.4e-5));assert y.min()>-1e-8
