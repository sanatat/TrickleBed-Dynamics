import json
from pathlib import Path

def test_paper_reference_not_claimed_validated():
    root=Path(__file__).resolve().parents[1]
    txt=(root/'results/v7_validation/VALIDATION_REPORT.md').read_text()
    assert 'NOT quantitatively validated' in txt
    assert 'NOT the published mass-fraction input' in txt

def test_reference_values_are_from_text_not_curve_data():
    root=Path(__file__).resolve().parents[1]
    data=json.loads((root/'results/v7_validation/paper_text_targets.json').read_text())
    assert data['reacted_axial_extent_chi1_approx_m']==1.0
    assert data['reacted_axial_extent_chi10_approx_m']==.25
    assert 'NOT digitized' in data['notes']

def test_reported_results_are_finite_and_lower_than_paper_temperature():
    root=Path(__file__).resolve().parents[1]
    arr=json.loads((root/'results/v7_validation/model_metrics.json').read_text())
    assert len(arr)==2
    assert arr[0]['model_peak_T_K']<500 and arr[1]['model_peak_T_K']<500
    assert arr[1]['model_outlet_liquid_proxy_conversion'] > arr[0]['model_outlet_liquid_proxy_conversion']
