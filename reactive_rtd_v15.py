"""Exploratory SIMULTANEOUS tracer/reaction integration on the v6 grid.
Same axial dispersion affects all dissolved mobile reactive species and tracer.
No measured transport coefficients; this is a synthetic coupling experiment.
"""
from dataclasses import dataclass
from scipy.integrate import solve_ivp
import numpy as np
from v6_model import Params, initial, rhs, unpack, volumes
from tracer_coupled_v14 import pulse, TracerConfig

@dataclass(frozen=True)
class Coupling:
    dispersion: float = 1.4e-5  # m2/s, assumed/illustrative; NOT validated
    pore_exchange: float = 0.0025  # 1/s approximately, assumed
    pulse_time: float = 35.0
    pulse_sigma: float = 7.0

def diffusive_divergence(C,D,dz):
    """Conservative no-flux end boundaries, per cross sectional area; shape (..., n)."""
    F=np.zeros((*C.shape[:-1],C.shape[-1]+1))
    F[...,1:-1]=-D*(C[...,1:]-C[...,:-1])/dz
    return (F[...,:-1]-F[...,1:])/dz

def full_rhs(t,y,p,c,tracer_on=True):
    m=14*p.n
    chem=y[:m]; mobile=y[m:m+p.n]; pore=y[m+p.n:]
    dy=rhs(t,chem,p)
    N,T,s=unpack(chem,p)
    V=volumes(s,p)
    C=N/np.maximum(V[:,None,:],1e-12)
    dz=p.length/p.n
    # Diffusion is volume-scaled into molar inventories in MOBILE LIQUID.
    diffusion=p.V_L*diffusive_divergence(C[0],c.dispersion,dz)
    dN=dy[:12*p.n].reshape(4,3,p.n)
    dN[0]+=diffusion
    L=mobile/p.V_L
    W=pore/(p.pore_vol*np.clip(s,p.s_min,p.s_max))
    # Use same volume-exchange closure as v6, with an additional mobile/pore transfer.
    flow_ratio=1. if t<p.t_step else p.flow_ratio
    u=p.u_L*flow_ratio
    inlet=pulse(t,TracerConfig(pulse_time=c.pulse_time,pulse_sigma=c.pulse_sigma)) if tracer_on else 0.0
    face=np.empty(p.n+1)
    face[0]=u*inlet
    face[1:]=u*L
    face[1:-1]+=-c.dispersion*(L[1:]-L[:-1])/dz
    flux=c.pore_exchange*np.clip(s,p.s_min,p.s_max)*(L-W)
    target=p.wet_adv if flow_ratio>=1 else p.wet_rec
    ds=(np.clip(target,p.s_min,p.s_max)-s)/p.cap_tau
    dvol=p.pore_vol*ds
    fill=np.maximum(dvol,0)*L+np.minimum(dvol,0)*W
    dmobile=(face[:-1]-face[1:])/dz-flux-fill
    dpore=flux+fill
    return np.r_[dy,dmobile,dpore]

def integrate(p,c,times,tracer_on=True):
    times=np.asarray(times,dtype=float)
    if times[0]!=0 or np.any(np.diff(times)<=0): raise ValueError('times must be increasing from 0')
    y0=np.r_[initial(p),np.zeros(2*p.n)]
    def solve(a,b,y,tg):
        z=solve_ivp(full_rhs,(a,b),y,t_eval=tg,args=(p,c,tracer_on),method='BDF',rtol=2e-7,atol=1e-9)
        if not z.success:raise RuntimeError(z.message)
        return z
    if 0<p.t_step<times[-1]:
        t1=np.unique(np.r_[times[times<p.t_step],p.t_step]); t2=np.unique(np.r_[p.t_step,times[times>p.t_step]])
        a=solve(0,p.t_step,y0,t1);b=solve(p.t_step,times[-1],a.y[:,-1],t2)
        ts=np.r_[a.t,b.t[1:]];Y=np.c_[a.y,b.y[:,1:]]
    else:
        z=solve(0,times[-1],y0,times);ts=z.t;Y=z.y
    m=14*p.n
    octene=[];prod=[]; temp=[]; tracer=[]; total=[]
    for j in range(len(ts)):
        N,T,s=unpack(Y[:m,j],p);V=volumes(s,p);C=N/np.maximum(V[:,None,:],1e-12)
        octene.append(C[0,0,-1]);prod.append(C[0,2,-1]);temp.append(T.max())
        tracer.append(Y[m+p.n-1,j]/p.V_L)
        total.append((Y[m:,j].sum())*p.length/p.n)
    return dict(time=ts,octene=np.array(octene),product=np.array(prod),peak_T=np.array(temp),tracer=np.array(tracer),tracer_inventory=np.array(total),raw=Y)
