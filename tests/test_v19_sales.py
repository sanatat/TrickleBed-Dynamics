import numpy as np
from pathlib import Path

def test_published_rows():
 a=np.genfromtxt(Path(__file__).resolve().parents[1]/'results/v19_sales_experimental/sales2005_table3.csv',delimiter=',',skip_header=1)
 assert a.shape==(6,5)
 assert np.allclose(a[:,0],[1,2,3,1,2,3])
 assert np.all((a[:,2]>0)&(a[:,2]<1))
 assert np.all(a[:,4]>0)
