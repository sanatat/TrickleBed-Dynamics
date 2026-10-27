from dataclasses import replace
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from paper_extension import PaperConfig,simulate,metrics,rate,rhs,initial,unpack

def test_rate_zero_without_reactant():
    c=PaperConfig()
    assert rate(np.array([0.,10.]),np.array([10.,0.]),np.array([400.,400.]),c).max()==0

def test_no_reaction_no_product():
    c=PaperConfig(n=8,t_end=300)
    _,Y=simulate(c,react=False)
    assert abs(Y.reshape(8,c.n,-1)[2]).max()<1e-6
    assert abs(Y.reshape(8,c.n,-1)[5]).max()<1e-6

def test_species_stoichiometry_closed_sources():
    c=PaperConfig(n=8,t_end=300)
    Y=initial(c);d=rhs(0,Y,c)
    dA,dH,dP,dAp,dHp,dPp,dT,df=unpack(d,c.n)
    vb=c.liquid_holdup;vp=(1-c.bed_porosity)*c.pellet_porosity*c.pore_liquid_fraction
    assert np.allclose(vb*(dA+dP)+vp*(dAp+dPp),0,atol=1e-9)

def test_both_chi_stable():
    for chi in (1.,10.):
        c=PaperConfig(n=8,chi=chi,t_end=350)
        t,Y=simulate(c)
        assert np.isfinite(Y).all() and Y.min()>-1e-5
        assert np.all(np.diff(t)>0)
        assert t[-1]==350

def test_hydrogen_and_carbon_balances_in_closed_cell():
    c=PaperConfig(n=1,u_in=0,kLa=0,UA=0,t_end=120,step_time=50)
    _,Y=simulate(c)
    A,H,P,Ap,Hp,Pp,T,f=Y.reshape(8,1,-1)
    vb=c.liquid_holdup;vp=(1-c.bed_porosity)*c.pellet_porosity*c.pore_liquid_fraction
    carb=vb*(A+P)+vp*(Ap+Pp)
    assert np.allclose(carb,carb[0,0],rtol=1e-6,atol=1e-4)
