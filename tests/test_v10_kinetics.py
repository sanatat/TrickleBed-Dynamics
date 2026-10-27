import numpy as np
import pytest
from kinetics_provenance_v10 import arrhenius_ratio,EA

def test_energies(): assert EA=={'TMP-1':34000.,'TMP-2':49000.}
def test_reference():
 for Ea in EA.values(): assert np.isclose(arrhenius_ratio(Ea,400.),1)
def test_temp_increase():
 for Ea in EA.values(): assert arrhenius_ratio(Ea,420)>1

def test_higher_ea_greater_temp_sensitivity(): assert arrhenius_ratio(49000,420)>arrhenius_ratio(34000,420)

def test_invalid():
 with pytest.raises(ValueError):arrhenius_ratio(34000,0)
