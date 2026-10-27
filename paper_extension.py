"""Traceable 1D reduced-order extensions motivated by Gorshkova et al. (2012).
Not the author's Fluent CFD implementation. ALL default fitted values are illustrative.
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp

@dataclass(frozen=True)
class PaperConfig:
    n: int = 24
    length: float = 1.5            # m, 2012 Table 2
    bed_porosity: float = 0.37     # 2012 Table 2
    pellet_porosity: float = 0.50  # 2012 Table 2
    particle_diameter: float = 0.003  # m, 2012 Table 2
    pore_diameter: float = 1e-5    # m, 2012 Table 2
    pressure_bar: float = 20.0     # 2012 Table 2; no momentum calculation
    T_in: float = 400.0           # K, 2012 Table 2
    A_in: float = 1500.0          # mol/m3, ILLUSTRATIVE
    H_in: float = 8.0             # mol/m3, ILLUSTRATIVE
    H_star: float = 12.0          # mol/m3, ILLUSTRATIVE Henry-equilibrium proxy
    u_in: float = 0.001           # m/s, ILLUSTRATIVE
    liquid_holdup: float = 0.25   # m3 per m3-bed, ILLUSTRATIVE
    pore_liquid_fraction: float = 0.65  # ILLUSTRATIVE
    kLa: float = 0.022            # 1/s, ILLUSTRATIVE
    k_pore: float = 0.008         # 1/s, ILLUSTRATIVE at chi=1
    k_ref: float = 0.014          # 1/s, ILLUSTRATIVE
    K_A: float = 1/1500          # m3/mol, ILLUSTRATIVE adsorption
    K_H: float = 0.05            # m3/mol, ILLUSTRATIVE adsorption
    E_a: float = 30000.0         # J/mol, ILLUSTRATIVE
    delta_H: float = -90000.0     # J/mol, ILLUSTRATIVE
    C_v: float = 1.45e6          # J/(m3-bed K), ILLUSTRATIVE
    UA: float = 220.0            # W/(m3-bed K), ILLUSTRATIVE
    tau_wet: float = 260.0       # s, ILLUSTRATIVE memory closure
    chi: float = 1.0             # dimensionless, 2012 comparison chi=1,10
    wet_init: float = 0.2        # ILLUSTRATIVE initial condition
    step_time: float = 200.0     # s, ILLUSTRATIVE
    flow_after: float = 0.7      # ILLUSTRATIVE
    t_end: float = 1200.0


def operating(t,c):
    ratio=1.0 if t<c.step_time else c.flow_after
    return ratio, np.clip(0.72*ratio**0.40, 0.05,0.95) # ILLUSTRATIVE, not Lappalainen correlation


def unpack(y,n):
    a=y.reshape(8,n)
    return tuple(a[i] for i in range(8))


def initial(c):
    n=c.n
    return np.r_[np.full(n,c.A_in),np.full(n,c.H_in),np.zeros(n),
                 np.full(n,c.A_in),np.full(n,c.H_in),np.zeros(n),
                 np.full(n,c.T_in),np.full(n,c.wet_init)]


def rate(A,H,T,c):
    """Single-olefin lumped analogue of Gorshkova Eq.24/25; NOT original fitted kinetics.
    Units: rate in mol/(m3-pore-liquid s) with effective k_ref adjusted accordingly.
    """
    A=np.maximum(A,0);H=np.maximum(H,0)
    k=c.k_ref*np.exp(-c.E_a/8.314462618*(1/T-1/c.T_in))
    return k*A*H/(1+c.K_A*A+c.K_H*H)**2


def rhs(t,y,c,react=True,supply=True):
    n=c.n;dz=c.length/n
    A,H,P,Ap,Hp,Pp,T,f=unpack(y,n)
    ratio,ftarget=operating(t,c)
    u=c.u_in*ratio
    vb=c.liquid_holdup
    vp=(1-c.bed_porosity)*c.pellet_porosity*c.pore_liquid_fraction
    fw=np.clip(f,0.0,1.0)
    # Eq. 11-like concentration-difference transfer; chi represents effective diffusion length scaling.
    # Same chi on all pore exchange is a simplification, not the complete five-interface CFD closure.
    kp=c.k_pore*c.chi
    JA=kp*vp*fw*(A-Ap)
    JH=kp*vp*fw*(H-Hp)
    JP=kp*vp*fw*(P-Pp)
    r=vp*fw*rate(Ap,Hp,T,c) if react else np.zeros(n)
    def adv(C,Cin):
        up=np.r_[Cin,C[:-1]]
        return u/dz*(up-C)
    gas=(c.kLa*vb*fw*(c.H_star-H)) if supply else np.zeros(n)
    dA=(adv(A,c.A_in)-JA)/vb
    dH=(adv(H,c.H_in)-JH+gas)/vb
    dP=(adv(P,0)-JP)/vb
    dAp=(JA-r)/vp
    dHp=(JH-r)/vp
    dPp=(JP+r)/vp
    upstreamT=np.r_[c.T_in,T[:-1]]
    dT=(c.C_v*u/dz*(upstreamT-T)-c.delta_H*r-c.UA*(T-c.T_in))/c.C_v
    df=(ftarget-f)/c.tau_wet
    return np.r_[dA,dH,dP,dAp,dHp,dPp,dT,df]


def simulate(c,react=True,supply=True):
    y0=initial(c)
    eval_all=np.linspace(0,c.t_end,241)
    def solve(a,b,y,eval_times):
        sol=solve_ivp(rhs,(a,b),y,args=(c,react,supply),method='BDF',t_eval=eval_times,
                      rtol=1e-7,atol=1e-9)
        if not sol.success: raise RuntimeError(sol.message)
        return sol
    if 0<c.step_time<c.t_end:
        t1=np.r_[eval_all[eval_all<c.step_time],c.step_time]
        t2=np.r_[c.step_time,eval_all[eval_all>c.step_time]]
        a=solve(0,c.step_time,y0,t1)
        b=solve(c.step_time,c.t_end,a.y[:,-1],t2)
        return np.r_[a.t,b.t[1:]],np.c_[a.y,b.y[:,1:]]
    a=solve(0,c.t_end,y0,eval_all)
    return a.t,a.y


def metrics(t,Y,c):
    A,H,P,Ap,Hp,Pp,T,f=Y.reshape(8,c.n,-1)
    return {'time_s':t,'conversion_out':1-A[-1]/c.A_in,
            'T_peak_K':T.max(axis=0),'T_out_K':T[-1],
            'wetting_mean':f.mean(axis=0),'product_out_molm3':P[-1]}
