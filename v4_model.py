"""Paper-inspired, CONSERVATIVE 1D 4-compartment prototype, NOT 2012 CFD reproduction.
Compartments: mobile liquid L, mobile gas G, pore liquid W, pore gas D.
Species A (isooctene proxy), H (hydrogen), P (isooctane proxy).
ALL transport, phase equilibrium, reaction, and wetting constants illustrative.
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp

@dataclass(frozen=True)
class Params:
    n: int=12
    length: float=1.5             # 2012 Table 2
    void: float=.37              # 2012 Table 2
    pellet_porosity: float=.5     # 2012 Table 2
    pellet_diameter: float=.003   # 2012 Table 2
    T_in: float=400.0            # 2012 Table 2
    pressure_bar: float=20.0     # 2012 Table 2; NOT solved
    V_L: float=.24               # illustration, volume per bed volume
    V_G: float=.13               # illustration
    V_W: float=.15               # illustration; FIXED pore liquid inventory
    V_D: float=.10               # illustration; FIXED pore gas inventory
    u_L: float=.0008             # superficial m/s, illustration
    u_G: float=.0015
    li_A: float=300.0
    li_H: float=8.0
    ga_H: float=500.0
    ga_A: float=6.0
    K_A: float=.02               # C_g / C_l, illustrative
    K_H: float=62.5
    K_P: float=.01
    k_GL: float=.0006            # mol/m3 bed/s per conc unit
    k_GD: float=.0003
    k_LW: float=.0012
    k_DW: float=.0007
    chi: float=1.
    k_ref: float=.002           # illustrative kinetic scale
    ads_A: float=.002
    ads_H: float=.04
    Ea: float=30000.
    dHr: float=-90000.
    heatcap: float=1.5e6
    heatloss: float=180.
    wet_init: float=.20
    # Deliberate wetting hysteresis: advancing and receding branches
    wet_adv: float=.84
    wet_rec: float=.46
    tau_wet: float=160.
    t_step: float=200.
    t_end: float=600.
    flow_ratio: float=.7
    
SPECIES=('A','H','P')

def inlets(p):
    return np.array([p.li_A,p.li_H,0.]),np.array([p.ga_A,p.ga_H,0.])

def initial(p):
    L,G=inlets(p)
    return np.concatenate([np.repeat(v[:,None],p.n,axis=1) for v in (L,G,L,G)]+[np.full((1,p.n),p.T_in),np.full((1,p.n),p.wet_init)]).ravel()

def unpack(y,p):
    a=np.asarray(y).reshape(14,p.n)
    return a[0:3],a[3:6],a[6:9],a[9:12],a[12],a[13]

def rate(A,H,T,p):
    A=np.maximum(A,0.);H=np.maximum(H,0.)
    k=p.k_ref*np.exp(np.clip(-p.Ea/8.314462618*(1/np.maximum(T,250.)-1/p.T_in),-20.,20.))
    return k*A*H/(1+p.ads_A*A+p.ads_H*H)**2

def rhs(t,y,p,flow=True,reaction=True,transfer=True):
    L,G,W,D,T,w=unpack(y,p)
    dz=p.length/p.n
    ratio=(1. if t<p.t_step else p.flow_ratio) if flow else 0.
    uL=p.u_L*ratio;uG=p.u_G if flow else 0.
    li,ga=inlets(p)
    def adv(c,inc,u):
        return u/dz*(np.column_stack((inc,c[:,:-1]))-c)
    dL=adv(L,li,uL)/p.V_L
    dG=adv(G,ga,uG)/p.V_G
    dW=np.zeros_like(W);dD=np.zeros_like(D)
    wet=np.clip(w,0.,1.)
    if transfer:
        K=np.array([p.K_A,p.K_H,p.K_P])[:,None]
        # each J is in mol/(m3 bed s); all pairwise terms cancel identically
        jGL=p.k_GL*(G-K*L)                  # G->L
        jGD=p.k_GD*(1.-wet)[None,:]*(G-D) # G->D
        jLW=p.k_LW*p.chi*wet[None,:]*(L-W) # L->W
        jDW=p.k_DW*wet[None,:]*(D-K*W)    # D->W
        dL+=(jGL-jLW)/p.V_L
        dG+=(-jGL-jGD)/p.V_G
        dW+=(jLW+jDW)/p.V_W
        dD+=(jGD-jDW)/p.V_D
    r=p.V_W*wet*rate(W[0],W[1],T,p) if reaction else np.zeros(p.n)
    dW[0]-=r/p.V_W;dW[1]-=r/p.V_W;dW[2]+=r/p.V_W
    adv_heat=(uL+uG)*.5/dz*(np.r_[p.T_in,T[:-1]]-T)
    dT=adv_heat+(-p.dHr*r-p.heatloss*(T-p.T_in))/p.heatcap
    # phenomenological wetting memory: different flow targets when increasing/reducing inlet flow
    target=p.wet_adv if ratio>=1. else p.wet_rec
    dw=(target-w)/p.tau_wet
    return np.concatenate((dL,dG,dW,dD,dT[None,:],dw[None,:])).ravel()

def simulate(p,flow=True,reaction=True,transfer=True):
    if not np.isclose(p.V_G+p.V_L,p.void):
        raise ValueError('Mobile holdups must add up to bed porosity')
    if min(p.V_G,p.V_L,p.V_W,p.V_D)<=0: raise ValueError('Holdups must be positive')
    y0=initial(p)
    def solve(t0,t1,y,tt):
        out=solve_ivp(rhs,(t0,t1),y,args=(p,flow,reaction,transfer),method='BDF',t_eval=tt,rtol=1e-7,atol=1e-9)
        if not out.success: raise RuntimeError(out.message)
        return out
    grid=np.linspace(0,p.t_end,121)
    if flow and 0<p.t_step<p.t_end:
        ts=np.r_[grid[grid<p.t_step],p.t_step]
        e=solve(0,p.t_step,y0,ts)
        ts2=np.r_[p.t_step,grid[grid>p.t_step]]
        f=solve(p.t_step,p.t_end,e.y[:,-1],ts2)
        return np.r_[e.t,f.t[1:]],np.c_[e.y,f.y[:,1:]]
    s=solve(0,p.t_end,y0,grid)
    return s.t,s.y

def species_inventory(y,p):
    a=np.asarray(y).reshape(14,p.n)
    return p.V_L*a[0:3]+p.V_G*a[3:6]+p.V_W*a[6:9]+p.V_D*a[9:12]
