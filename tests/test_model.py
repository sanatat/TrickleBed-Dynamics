import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from model import Config, simulate, unpack, metrics

def test_all_states_finite_and_positive():
    c=Config(n=10,t_end=400)
    t,y=simulate(c,.2)
    assert np.isfinite(y).all()
    assert y.min()>=-1e-5
    assert t[-1]==400

def test_no_reaction_no_product():
    c=Config(n=8,t_end=350)
    _,y=simulate(c,.2,reaction=False)
    parts=unpack(y,c.n)
    assert max(np.max(np.abs(parts[i])) for i in (2,5)) < 1e-7

def test_initial_wetting_paths_diverge():
    c=Config(n=10,t_end=400)
    _,a=simulate(c,.2)
    _,b=simulate(c,.88)
    ma=metrics(None,a,c);mb=metrics(None,b,c)
    assert abs(ma['temperature_peak_K'].max()-mb['temperature_peak_K'].max())>1e-3

def test_carbon_species_initial_conservation_no_reaction():
    c=Config(n=8,t_end=350)
    _,y=simulate(c,.2,reaction=False)
    A,H,P,Ap,Hp,Pp,T,f=unpack(y,c.n)
    assert np.max(np.abs(P)) < 1e-7
    assert np.max(np.abs(Pp)) < 1e-7
