import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from cross_flow_v23 import residence_time,run
def test_times():
 assert abs(residence_time(.8,.17,.001)-136)<1e-8
 assert abs(residence_time(.8,.21,.002)-84)<1e-8
def test_invalid():
 import pytest
 with pytest.raises(ValueError):residence_time(.8,-.1,.001)
def test_run():
 m=run()
 assert m['refit_to_figure6'] is False and m['n_approx']>10
 assert m['rmse_cross_flow_approx']>=0
