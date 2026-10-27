from dataclasses import replace
import numpy as np
from v4_model import Params,initial,rhs,simulate,species_inventory,unpack

def test_shape_and_finite():
 p=Params(n=4,t_end=20,t_step=10)
 t,y=simulate(p)
 assert y.shape==(14*p.n,len(t))
 assert np.isfinite(y).all()

def test_closed_mass_stoichiometry():
 p=Params(n=5)
 y=initial(p)
 # vary compartments deliberately to force exchange in all four reservoirs
 arr=y.reshape(14,p.n);arr[0]*=1.2;arr[4]*=.8;arr[10]*=1.3
 d=rhs(0,y,p,flow=False).reshape(14,p.n)
 inv=species_inventory(d.ravel(),p)
 assert np.max(np.abs(inv[0]+inv[2]))<1e-9
 assert np.max(np.abs(inv[1]+inv[2]))<1e-9

def test_individual_species_conserved_without_reaction():
 p=Params(n=6)
 d=rhs(0,initial(p),p,flow=False,reaction=False)
 assert np.max(np.abs(species_inventory(d,p)))<1e-10

def test_wetting_stays_bounded():
 p=Params(n=4,t_end=260)
 t,y=simulate(p)
 w=y.reshape(14,p.n,-1)[13]
 assert (w>=0).all() and (w<=1).all()

def test_step_event_and_hysteresis():
 p=Params(n=4,t_end=400)
 t,y=simulate(p)
 w=y.reshape(14,p.n,-1)[13,0]
 assert w[np.argmin(abs(t-200))]>w[-1]

def test_zero_reaction_no_heat_excursion():
 p=Params(n=4,t_end=100)
 _,y=simulate(p,reaction=False)
 assert np.max(np.abs(y.reshape(14,p.n,-1)[12]-p.T_in))<1e-6
