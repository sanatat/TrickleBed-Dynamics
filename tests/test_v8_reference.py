from reference_audit import audit, complete_mass_fractions, component_molar_flows, DATA, MW_KG_MOL
import math

def test_closes_mass_fractions():
    for phase in ('liquid','gas'):
        assert abs(sum(complete_mass_fractions(DATA[f'{phase}_mass_fractions']).values())-1)<1e-12

def test_flow_mass_recovered():
    for phase in ('liquid','gas'):
        flow=DATA[f'{phase}_mass_flow_kg_s']
        mol=component_molar_flows(flow,DATA[f'{phase}_mass_fractions'])
        assert math.isclose(sum(mol[k]*MW_KG_MOL[k] for k in mol),flow,rel_tol=1e-12)

def test_total_flow_additivity():
    d=audit()
    for k in MW_KG_MOL:
        assert math.isclose(d['total_feed_mol_s'][k],d['liquid_feed_mol_s'][k]+d['gas_feed_mol_s'][k])

def test_reaction_limit():
    d=audit();assert d['stoichiometric_upper_bound_reaction_mol_s']<=d['total_feed_mol_s']['octene']
    assert d['stoichiometric_upper_bound_reaction_mol_s']<=d['total_feed_mol_s']['H2']
    assert 0<=d['theoretical_max_octene_conversion_at_inlet_stoichiometry']<=1

def test_bad_mass_fractions_rejected():
    import pytest
    with pytest.raises(ValueError):complete_mass_fractions({'octene':.9,'H2':.2})
