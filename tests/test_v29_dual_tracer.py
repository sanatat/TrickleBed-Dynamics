"""v2.9 non-experimental two-tracer verification tests."""
import numpy as np
from dataclasses import replace
from dual_tracer_v29 import parameterized, fisher, fiducial_diagnostics, sensitivity, predicted
from heterogeneous_rtd_v27 import Params, check_global_inventory


def test_kcl_nonadsorbing_explicit_assumption():
    x=np.log([0.2,8.2e-5])
    p=parameterized(Params(),x,'KCl_assumed_nonadsorbing')
    assert p.K_ads_m3_kg==0
    assert p.accessible_pore_fraction==0.2


def test_kcl_output_does_not_depend_on_nagl_adsorption():
    p=replace(Params(),n_axial=9,n_radial=2)
    t=np.linspace(0,700,17);ci=np.where(t<=140,1.,0.)
    a=predicted(t,ci,p,np.log([.3,2e-5]),'KCl_assumed_nonadsorbing')
    b=predicted(t,ci,p,np.log([.3,3e-4]),'KCl_assumed_nonadsorbing')
    assert np.allclose(a,b,atol=1e-9)


def test_fisher_information_gains_are_nonnegative():
    p=replace(Params(),n_axial=9,n_radial=2)
    t=np.linspace(0,600,15);ci=np.where(t<100,1.,0.)
    theta=np.log([.2,8.2e-5])
    Jn=sensitivity(t,ci,p,theta,'NaGl')
    Jk=sensitivity(t,ci,p,theta,'KCl_assumed_nonadsorbing')
    F=fisher(Jn);Fdual=F+fisher(Jk)
    assert np.linalg.eigvalsh(Fdual-F).min()>-1e-10
    assert fiducial_diagnostics(Fdual)['smallest_eigenvalue']>=fiducial_diagnostics(F)['smallest_eigenvalue']-1e-9


def test_mass_conservation_at_low_accessibility():
    p=replace(Params(),accessible_pore_fraction=.2)
    stationary,inlet=check_global_inventory(p)
    assert stationary<1e-11 and inlet<1e-11


def test_real_data_not_labeled_as_kcl():
    from pathlib import Path
    module=Path(__file__).resolve().parents[1]/'dual_tracer_v29.py'
    s=module.read_text()
    assert 'KCl trajectories SYNTHETIC' in s
