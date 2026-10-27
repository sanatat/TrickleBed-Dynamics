from dataclasses import replace
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from three_phase import Config,init,rhs,simulate,conserved_totals

def test_mobile_holdups():
 c=Config()
 assert c.liquid_fraction+c.gas_fraction==c.porosity

def test_closed_system_carbon_conserved():
 c=Config(n=4,t_end=100,wet_relax_time=60)
 t,Y=simulate(c,flow=False)
 tot=conserved_totals(Y,c)
 carbon=tot[0]+tot[2]
 assert np.allclose(carbon,carbon[0],atol=1e-5,rtol=1e-7)

def test_closed_system_h2_stoichiometry():
 c=Config(n=4,t_end=100)
 t,Y=simulate(c,flow=False)
 tot=conserved_totals(Y,c)
 assert np.allclose(tot[1]+tot[2],(tot[1]+tot[2])[0],atol=1e-5,rtol=1e-7)

def test_no_reaction_preserves_products():
 c=Config(n=4,t_end=60)
 _,Y=simulate(c,reaction=False)
 a=Y.reshape(11,c.n,-1)
 assert np.abs(a[[2,5,8]]).max()<1e-6

def test_finite_nonnegative_and_monotonic_time():
 for chi in (1,10):
  c=Config(n=4,t_end=80,step_time=40,chi=chi)
  t,Y=simulate(c)
  assert np.isfinite(Y).all()
  assert Y.min()>=-1e-5
  assert np.all(np.diff(t)>0)

def test_closed_no_reaction_individual_species_conserved():
 c=Config(n=4,t_end=80)
 _,Y=simulate(c,flow=False,reaction=False)
 tot=conserved_totals(Y,c)
 assert np.allclose(tot,tot[:,0,None],rtol=1e-7,atol=1e-5)
