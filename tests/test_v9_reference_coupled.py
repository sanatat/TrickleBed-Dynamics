import numpy as np
import pytest
from reference_coupled_v9 import map_reference_feed,run
from reference_audit import audit

def test_molar_flux_at_inlet():
 p,info=map_reference_feed()
 a=audit()
 assert np.isclose(p.li_A*p.u_L, a['liquid_feed_mol_s']['octene'])
 assert np.isclose(p.li_H*p.u_L, a['liquid_feed_mol_s']['H2'])
 assert np.isclose(p.ga_A*p.u_G, a['gas_feed_mol_s']['octene'])
 assert np.isclose(p.ga_H*p.u_G, a['gas_feed_mol_s']['H2'])

def test_density_changes_liquid_velocity_only():
 p,_=map_reference_feed(600);q,_=map_reference_feed(800)
 assert p.u_L>q.u_L and np.isclose(p.u_G,q.u_G)

def test_bad_density_rejected():
 with pytest.raises(ValueError):map_reference_feed(0)

def test_source_feed_temperature_pressure():
 p,_=map_reference_feed()
 assert p.T_in==400 and p.pressure_bar==20
