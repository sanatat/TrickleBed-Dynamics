"""Trace the marker geometry already embedded in Sales et al. (2005), PDF p.8.

No OCR, image classification, or invented raw measurements. Figure 5 dense
vector symbols are resampled as a DRAWN TRACE, not independent experiments.
Figure 6 triangle centers are measured observations as depicted by authors.
Reproduction requires the user's downloaded source PDF via --pdf.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
import fitz

SOURCE='Sales et al. (2005), Brazilian Journal of Chemical Engineering 22(3), 443-452, Fig. 5-6, printed p.450'
DOI='10.1590/S0104-66322005000300012'
# Axis limits derive from PDF vector tick-mark paths, NOT fitted to output data.
AXES={
  'figure5': {'x0':101.52,'x1800':258.48,'y1':463.04,'y0':642.80},
  'figure6': {'x0':351.60,'x1800':510.72,'y1':462.56,'y0':638.96},
}

def tick_calibrate(x,y,axes):
    """PDF point (top-left origin) -> seconds and normalized concentration."""
    a=axes
    return 1800*(np.asarray(x)-a['x0'])/(a['x1800']-a['x0']), (a['y0']-np.asarray(y))/(a['y0']-a['y1'])

def extract(pdf_path):
    doc=fitz.open(str(pdf_path))
    if len(doc)<8: raise ValueError('source PDF needs page 8 (printed page 450)')
    drawings=doc[7].get_drawings()
    series={}
    # The first two drawn datasets are extremely dense vector marker sequences
    # (same trace rendered many times); the source contains no raw instrument file.
    for name,k in [('fig5_experimental_exit',3),('fig5_experimental_entrance',4)]:
        xx=[];yy=[]
        for d in drawings:
            b=d['rect'];
            if d['type']!='fs' or len(d['items'])!=k:continue
            if not (101<=b.x0<258.5 and 457<b.y0<650):continue
            # Keep only small markers; oversized triangle in legend is excluded.
            if (k==3 and b.width>1.70) or (k==4 and b.width>1.75):continue
            xx.append((b.x0+b.x1)/2);yy.append((b.y0+b.y1)/2)
        raw=np.column_stack((xx,yy));raw=raw[np.argsort(raw[:,0])]
        # Quantize to original PDF's x resolution, merge duplicates/overprints.
        ux=np.unique(np.round(raw[:,0],2)); uy=np.array([np.median(raw[np.abs(raw[:,0]-x)<0.011,1]) for x in ux]);
        t,y=tick_calibrate(ux,uy,AXES['figure5'])
        mask=(t>=-0.5)&(t<=1800.5)
        series[name]={'pdf_xy':np.column_stack((ux[mask],uy[mask])), 'xy':np.column_stack((t[mask],y[mask])), 'source_mark_count':len(raw), 'kind':'resampled graphical trace (not independent point observations)'}
    # Figure 6 contains discrete large triangular EXPERIMENTAL markers.
    pts=[]
    for d in drawings:
        b=d['rect']
        if d['type']=='fs' and len(d['items'])==3 and 351.3<=b.x0<511 and 450<b.y0<644 and b.width>3.0:
            x=(b.x0+b.x1)/2;y=(b.y0+b.y1)/2
            # The label-legend's triangle sits at x about 429, y~571.
            if 407<x<478 and 541<y<600:continue
            pts.append((x,y))
    raw=np.unique(np.round(np.array(pts),2),axis=0);raw=raw[np.argsort(raw[:,0])]
    t,y=tick_calibrate(raw[:,0],raw[:,1],AXES['figure6']);mask=(t>=0)&(t<=1800.5)
    series['fig6_experimental_exit']={'pdf_xy':raw[mask],'xy':np.column_stack((t[mask],y[mask])),'source_mark_count':len(raw),'kind':'printed discrete experimental triangle symbols'}
    if len(series['fig6_experimental_exit']['xy'])<30:raise ValueError('Insufficient triangle markers; verify this is Sales 2005 original PDF')
    # The original authors' MODEL curve is a stroked PDF polyline, not observations.
    # Reconstruct it separately for a source-figure consistency check only.
    linepts=[]
    for d in drawings:
        b=d['rect']
        if d['type']!='s' or len(d['items'])<10:continue
        if not (351.0<b.x0<510.7 and 460<b.y0<645):continue
        if not (d['items'][0][0]=='l'): continue
        for item in d['items']:
            if item[0]=='l':linepts.append((item[1].x,item[1].y));linepts.append((item[2].x,item[2].y))
    if linepts:
        pts=np.unique(np.round(np.array(linepts),3),axis=0)
        pts=pts[np.argsort(pts[:,0])]
        # Source polyline contains duplicate x due to vertical drawing steps.
        px=np.unique(pts[:,0]);py=np.array([np.median(pts[pts[:,0]==x,1]) for x in px]);
        t,c=tick_calibrate(px,py,AXES['figure6']);valid=(t>=0)&(t<=1800.5)
        series['fig6_published_fitted_model']={'pdf_xy':np.column_stack((px[valid],py[valid])),
           'xy':np.column_stack((t[valid],c[valid])),'source_mark_count':len(pts),
           'kind':'authors original plotted model polyline (NOT new observations)'}
    return series

def save_data(pdf_path, output_dir, sample_step=30.):
    output_dir=Path(output_dir);output_dir.mkdir(parents=True,exist_ok=True)
    series=extract(pdf_path)
    info={
      'source':SOURCE,'doi':DOI,'page_index_0_based':7,
      'sha256_source_pdf':hashlib.sha256(Path(pdf_path).read_bytes()).hexdigest(),
      'pdf_coordinate_unit':'PDF points with origin top-left',
      'axes_pdf_points':AXES,
      'method':'read embedded vector drawing marker shapes with PyMuPDF; calibrate using x-axis ticks 0,1800 s and y-axis ticks 0,1',
      'not_raw_instrument_data':True,
      'limits':'Figure 5 dense vector markers are a plotted continuous trace. Sampling them yields correlated reconstructions, not independent experimental points. Figure 6 triangle symbols carry finite plotting-marker uncertainty; inlet waveform for Figure 6 is not separately published.',
      'uncertainty_conservative':'Use ±1.0 PDF pt for Fig 5 and ±2.2 PDF pt for Fig 6 symbol positioning, PLUS unquantified digitization/author plotting systematic errors. These are not instrument uncertainty or confidence intervals.',
      'series':{},
    }
    for name,dat in series.items():
        coords=dat['pdf_xy'];a=dat['xy'];
        if name.startswith('fig5'):
            # Create 61 evenly spaced points as a graphical *resampling*.
            grid=np.arange(0,1800+0.001,sample_step)
            keep=(grid>=min(a[:,0]))&(grid<=max(a[:,0])); grid=grid[keep]
            y=np.interp(grid,a[:,0],a[:,1]);
            xp=np.interp(grid,a[:,0],coords[:,0]);yp=np.interp(grid,a[:,0],coords[:,1]);
            table=np.column_stack((grid,y,xp,yp)); uncertainty_pt=1.0
        else:
            table=np.column_stack((a,coords));uncertainty_pt=2.2
        with (output_dir/f'{name}_vector.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['time_s','normalized_signal','pdf_x_pt','pdf_y_pt']);w.writerows(table)
        key='figure5' if name.startswith('fig5') else 'figure6'; bounds=AXES[key]
        info['series'][name]={'embedded_vector_markers':dat['source_mark_count'],'exported_rows':len(table),'kind':dat['kind'],
          'time_point_approx_uncertainty_s':1800*uncertainty_pt/(bounds['x1800']-bounds['x0']),
          'fraction_approx_uncertainty':uncertainty_pt/(bounds['y0']-bounds['y1'])}
    (output_dir/'PROVENANCE.json').write_text(json.dumps(info,indent=2)+'\n')
    return info

def main():
    p=argparse.ArgumentParser();p.add_argument('--pdf',required=True);p.add_argument('--out',default='results/v26_vector_extraction');a=p.parse_args()
    print(json.dumps(save_data(a.pdf,a.out),indent=2))
if __name__=='__main__':main()
