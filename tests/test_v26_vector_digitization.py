"""Deterministic verification of figure-calibration and exported point provenance."""
from pathlib import Path
import os
import numpy as np
import pytest
from vector_digitize_v26 import AXES,tick_calibrate,extract
from fit_vector_rtd_v26 import quantile_crossing
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'v26_vector_extraction'

def test_axis_endpoints():
    for a in AXES.values():
        t,c=tick_calibrate([a['x0'],a['x1800']],[a['y0'],a['y1']],a)
        assert np.allclose(t,[0,1800],rtol=0,atol=1e-10)
        assert np.allclose(c,[0,1],rtol=0,atol=1e-10)

def test_printed_markers_saved():
    for name in ['fig5_experimental_entrance','fig5_experimental_exit','fig6_experimental_exit']:
        a=np.loadtxt(OUT/f'{name}_vector.csv',delimiter=',',skiprows=1)
        assert len(a)>=50
        assert np.all(np.diff(a[:,0])>0)
        assert (a[:,0].min()>=0) and a[:,0].max()<=1800.5
        assert a[:,1].min()>-.02 and a[:,1].max()<1.02

def test_half_concentration_levels():
    def load(name):
        a=np.loadtxt(OUT/f'{name}_vector.csv',delimiter=',',skiprows=1)
        return a[:,0],np.clip(a[:,1],0,1)
    ti,ci=load('fig5_experimental_entrance');to,co=load('fig5_experimental_exit'); t6,c6=load('fig6_experimental_exit')
    assert quantile_crossing(ti,ci,.5)<quantile_crossing(to,co,.5)
    assert 150<quantile_crossing(t6,c6,.5)<350

def test_pdf_geometry_when_source_available():
    source = os.environ.get('SALES2005_PDF', '')
    if not source: pytest.skip('set SALES2005_PDF to the locally downloaded original PDF to test source extraction')
    pdf=Path(source)
    if not pdf.is_file(): pytest.skip('SALES2005_PDF does not exist; optional source extraction test skipped')
    res=extract(pdf)
    assert len(res['fig6_experimental_exit']['xy'])==61
    assert res['fig5_experimental_entrance']['source_mark_count']>900

def test_data_is_declared_graph_derived():
    import json
    p=json.loads((OUT/'PROVENANCE.json').read_text())
    assert p['not_raw_instrument_data'] is True
    assert 'inlet waveform' in p['limits']
