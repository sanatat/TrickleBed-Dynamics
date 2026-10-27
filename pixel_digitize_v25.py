"""Reproducible coordinate-calibration utility for published Figure 5/6.

Inputs are manually *clicked pixel coordinates*, on a local plot image.
This script never manufactures published observations or performs OCR.
Use axes reference pixels based on x=0,1800 sec and y=0,1 mass fraction.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SOURCE_URL='https://doi.org/10.1590/S0104-66322005000300012'

def calibrate(points, bounds, pixel_error=2.):
    """bounds=(pixel_x_at_0,pixel_x_at_1800,pixel_y_at_0,pixel_y_at_1).
    Follows plot coordinate conventions (y screen grows downward).
    """
    a,b,c,d=map(float,bounds)
    if not (b>a and d<c):raise ValueError('bad image axis calibration')
    p=np.asarray(points,dtype=float)
    if p.ndim!=2 or p.shape[1]!=2:raise ValueError('points must be Nx2')
    result=np.column_stack((1800*(p[:,0]-a)/(b-a),(p[:,1]-c)/(d-c)))
    u_t=1800*pixel_error/(b-a);u_y=pixel_error/abs(d-c)
    return result,(u_t,u_y)

def read_pixels(path):
    with open(path,newline='') as f:
        r=csv.DictReader(f)
        if not {'series','pixel_x','pixel_y'}.issubset(r.fieldnames):raise ValueError('required series,pixel_x,pixel_y')
        return [(row['series'],float(row['pixel_x']),float(row['pixel_y'])) for row in r]

def export(input_csv,output_csv,bounds,pixel_error=2):
    rows=read_pixels(input_csv)
    calibrated,unc=calibrate([(x,y) for _,x,y in rows],bounds,pixel_error)
    p=Path(output_csv);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['series','time_s','mass_fraction','pixel_x','pixel_y','uncertainty_time_s','uncertainty_fraction'])
        for (name,x,y),(t,frac) in zip(rows,calibrated):w.writerow([name,float(t),float(frac),x,y,*unc])
    meta={'data_class':'pixel-digitized publication figure; not raw measurements','source':SOURCE_URL,'bounds_px':list(map(float,bounds)),'pixel_error_assumed_px':pixel_error,'approx_coordinate_uncertainty':{'time_s':unc[0],'mass_fraction':unc[1]},'rows':len(rows),'input_path':str(input_csv),'limitations':'manual click + pixel-error only; does not cover instrument or line thickness ambiguity'}
    p.with_suffix('.provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta

def click_image(image_path, output_csv):
    """Runs only in a local desktop Matplotlib interactive backend."""
    import matplotlib
    matplotlib.use('TkAgg',force=True)
    import matplotlib.pyplot as plt
    im=plt.imread(image_path)
    fig,ax=plt.subplots(figsize=(12,8));ax.imshow(im);ax.set_title('Click curve points left-to-right; press ENTER to finish')
    pts=plt.ginput(n=-1,timeout=0);plt.close(fig)
    with open(output_csv,'w',newline='') as f:
        w=csv.writer(f);w.writerow(['series','pixel_x','pixel_y']);w.writerows([('experimental_exit',x,y) for x,y in pts])
    return len(pts)

def main():
    pa=argparse.ArgumentParser(description=__doc__)
    pa.add_argument('--image',help='Locally saved/cropped source plot PNG')
    pa.add_argument('--click',action='store_true',help='Manually click a plot and save points; requires Tk')
    pa.add_argument('--pixels',required=True,help='CSV of labeled pixel coordinates')
    pa.add_argument('--output',required=True)
    pa.add_argument('--axis',nargs=4,type=float,metavar=('X0','X1800','Y0','Y1'),help='Four measured axis reference positions in pixels')
    pa.add_argument('--pixel-error',type=float,default=2.)
    args=pa.parse_args()
    if args.click:
        if not args.image:pa.error('--image is required with --click')
        n=click_image(args.image,args.pixels)
        print(f'Saved {n} clicked points to {args.pixels}')
    if args.axis:
        print(json.dumps(export(args.pixels,args.output,args.axis,args.pixel_error),indent=2))
    elif not args.click:pa.error('Provide --axis for calibration or --click for point entry')
if __name__=='__main__':main()
