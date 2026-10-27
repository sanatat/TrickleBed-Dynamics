import numpy as np
from robust_design_v18 import ENSEMBLE,CANDIDATES,FLOWS,design_records

def test_ensemble_positive():
    assert len(ENSEMBLE)==5 and all(all(x>0 for x in t) for t in ENSEMBLE)

def test_number_of_designs():
    assert len(CANDIDATES)==21

def test_worst_case_score():
    j={float(f):np.eye(2) for f in FLOWS}
    rec=design_records([j,j])
    assert len(rec)==21 and all(np.isclose(r['worst'],np.log(4.0)) for r in rec)
