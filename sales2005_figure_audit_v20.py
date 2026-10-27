"""Refit approximate manually read published-figure points; NOT raw experimental data."""
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammaincc
from matplotlib import pyplot as plt
p=Path(__file__).parent/'results'/'v20_sales_digitization'
a=np.genfromtxt(p/'figure6_approx_experimental.csv',delimiter=',',skip_header=1)
t,y=a[:,0],a[:,1]
def predicted(v):return gammaincc(v[1],np.maximum(t-v[0],0)/v[2])
r=least_squares(lambda v: predicted(v)-y,[110,3,60],bounds=([0,.3,5],[450,30,500]))
print('Estimated phenomenological params:',r.x,'approx curve RMSE',np.mean((predicted(r.x)-y)**2)**.5)
