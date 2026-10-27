from pathlib import Path
import numpy as np

def test_fig6_curve_bounds():
 p=Path(__file__).resolve().parents[1]/'results'/'v20_sales_digitization'/'figure6_approx_experimental.csv'
 a=np.genfromtxt(p,delimiter=',',skip_header=1)
 assert len(a)==19 and np.all(np.diff(a[:,0])>0) and np.all((a[:,1]>=0)&(a[:,1]<=1))

def test_fig5_two_curves_and_times():
 p=Path(__file__).resolve().parents[1]/'results'/'v20_sales_digitization'/'figure5_approx_experimental.csv'
 a=np.genfromtxt(p,delimiter=',',skip_header=1)
 assert len(a)==13 and a.shape[1]==3 and np.all(np.diff(a[:,0])>0)
