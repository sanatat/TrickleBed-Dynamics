"""Explicit mobile gas + liquid + stagnant intrapellet liquid, fixed catalyst.
A paper-inspired 1-D prototype, NOT a reconstruction of the 2012 Fluent model.
Concentrations: mol/m^3 of own phase; fluxes: mol/(m^3-bed s).
Species: A=octene, H=hydrogen, P=octane; gas H and vapor A/P.
All transfer, kinetic, thermodynamic and wetting constants are ILLUSTRATIVE.
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp

@dataclass(frozen=True)
class Config:
    n: int = 16
    length: float = 1.5               # m: Gorshkova Table 2
    porosity: float = 0.37           # bed, Table 2
    pellet_porosity: float = 0.5     # Table 2
    pellet_diameter: float = .003    # m, Table 2
    pressure_bar: float = 20.0       # input reference, not solved
    temp_in: float = 400.0           # K, Table 2
    liquid_fraction: float = .24     # ILLUSTRATIVE fraction of bed volume
    gas_fraction: float = .13        # ILLUSTRATIVE, should equal porosity-liquid_fraction
    pore_fraction: float = .18       # ILLUSTRATIVE, fraction of bed occupied by pore liquid
    liquid_velocity: float = .001    # m/s superficial ILLUSTRATIVE
    gas_velocity: float = .002       # m/s superficial ILLUSTRATIVE
    inlet_A_liquid: float = 300.0    # mol/m3-liquid ILLUSTRATIVE
    inlet_H_liquid: float = 8.0      # mol/m3-liquid ILLUSTRATIVE
    inlet_P_liquid: float = 0.0
    inlet_H_gas: float = 500.0       # mol/m3-gas ILLUSTRATIVE
    inlet_A_gas: float = 6.0         # mol/m3-gas ILLUSTRATIVE
    inlet_P_gas: float = 0.0
    # Equilibrium Cgas=K*C_liquid; proxies, not real flash or validated Henry/Antoine
    K_A: float = .02
    K_H: float = 62.5
    K_P: float = .01
    k_gl: float = .0008              # 1/s effective over bed vol ILLUSTRATIVE
    k_pore: float = .0015           # 1/s ILLUSTRATIVE
    chi: float = 1.0                # scaled diffusion length proxy
    k_ref: float = .002             # kinetic coefficient calibrated nowhere
    ads_A: float = .002             # m3/mol ILLUSTRATIVE
    ads_H: float = .04              # m3/mol ILLUSTRATIVE
    activation_energy: float = 30000. # J/mol ILLUSTRATIVE
    reaction_enthalpy: float = -90000. # J/mol ILLUSTRATIVE
    volumetric_heat_capacity: float = 1.5e6 # J/(m3-bed K) ILLUSTRATIVE
    heat_loss: float = 200.         # W/(m3-bed K) ILLUSTRATIVE
    wet_relax_time: float = 200.    # s ILLUSTRATIVE; no genuine hysteresis
    initial_wet: float = .2
    wet_equilibrium: float = .72
    step_time: float = 200.         # s
    flow_after: float = .7
    t_end: float = 600.

SPECIES=('A','H','P')

def inlet(c):
    return np.array([c.inlet_A_liquid,c.inlet_H_liquid,c.inlet_P_liquid]), np.array([c.inlet_A_gas,c.inlet_H_gas,c.inlet_P_gas])

def init(c):
    li,ga=inlet(c)
    L=np.repeat(li[:,None],c.n,axis=1)
    G=np.repeat(ga[:,None],c.n,axis=1)
    W=L.copy()
    return np.concatenate([L,G,W,np.full((1,c.n),c.temp_in),np.full((1,c.n),c.initial_wet)]).ravel()

def split(y,c):
    a=np.asarray(y).reshape(11,c.n)
    return a[:3],a[3:6],a[6:9],a[9],a[10]

def kinetic(A,H,T,c):
    A=np.maximum(A,0);H=np.maximum(H,0)
    k=c.k_ref*np.exp(np.clip(-c.activation_energy/8.314462618*(1/np.maximum(T,250)-1/c.temp_in),-20,20))
    return k*A*H/(1+c.ads_A*A+c.ads_H*H)**2

def rhs(t,y,c,reaction=True,exchange=True,flow=True):
    L,G,W,T,f=split(y,c)
    dz=c.length/c.n
    li,ga=inlet(c)
    ratio=(1.0 if t<c.step_time else c.flow_after) if flow else 0.
    ul=c.liquid_velocity*ratio; ug=c.gas_velocity if flow else 0.
    wet=np.clip(f,0,1)
    target=c.wet_equilibrium*(ratio**.4 if flow else 1.)
    # Forward finite-volume advection of conserved phase inventories
    def adv(C,Cin,velocity):
        return velocity/dz*(np.r_[Cin,C[:-1]]-C)
    dL=np.array([adv(L[i],li[i],ul) for i in range(3)])/c.liquid_fraction
    dG=np.array([adv(G[i],ga[i],ug) for i in range(3)])/c.gas_fraction
    dW=np.zeros_like(W)
    if exchange:
        K=np.array([c.K_A,c.K_H,c.K_P])[:,None]
        # positive Jgl = gas -> liquid, obeys exact species conservation
        Jgl=c.k_gl*(G-K*L)
        # positive Jlw = bulk liquid -> pore liquid
        Jlw=c.k_pore*c.chi*wet[None,:]*(L-W)
        dL+=(Jgl-Jlw)/c.liquid_fraction
        dG-=Jgl/c.gas_fraction
        dW+=Jlw/c.pore_fraction
    rate=c.pore_fraction*wet*kinetic(W[0],W[1],T,c) if reaction else np.zeros(c.n)
    dW[0]-=rate/c.pore_fraction
    dW[1]-=rate/c.pore_fraction
    dW[2]+=rate/c.pore_fraction
    # pseudo-mixture thermal closure with advective speed; NOT individual phase temperatures
    heat_adv=(ul+ug)*.5/dz*(np.r_[c.temp_in,T[:-1]]-T)
    dT=heat_adv+(-c.reaction_enthalpy*rate-c.heat_loss*(T-c.temp_in))/c.volumetric_heat_capacity
    df=(target-f)/c.wet_relax_time
    return np.concatenate([dL,dG,dW,dT[None,:],df[None,:]]).ravel()

def simulate(c,reaction=True,exchange=True,flow=True):
    if abs(c.liquid_fraction+c.gas_fraction-c.porosity)>1e-9: raise ValueError('Mobile phase fractions must sum to bed porosity')
    x0=init(c)
    grid=np.linspace(0,c.t_end,121)
    if flow and 0<c.step_time<c.t_end:
        g1=np.r_[grid[grid<c.step_time],c.step_time]
        g2=np.r_[c.step_time,grid[grid>c.step_time]]
        s1=solve_ivp(rhs,(0,c.step_time),x0,args=(c,reaction,exchange,flow),method='BDF',t_eval=g1,rtol=2e-6,atol=1e-8)
        if not s1.success:raise RuntimeError(s1.message)
        s2=solve_ivp(rhs,(c.step_time,c.t_end),s1.y[:,-1],args=(c,reaction,exchange,flow),method='BDF',t_eval=g2,rtol=2e-6,atol=1e-8)
        if not s2.success:raise RuntimeError(s2.message)
        return np.r_[s1.t,s2.t[1:]],np.c_[s1.y,s2.y[:,1:]]
    s=solve_ivp(rhs,(0,c.t_end),x0,args=(c,reaction,exchange,flow),method='BDF',t_eval=grid,rtol=2e-6,atol=1e-8)
    if not s.success:raise RuntimeError(s.message)
    return s.t,s.y

def conserved_totals(Y,c):
    V=np.array([c.liquid_fraction,c.gas_fraction,c.pore_fraction])
    a=Y.reshape(11,c.n,-1)
    return np.einsum('p,pcnt->cnt',V,a[:9].reshape(3,3,c.n,-1)).sum(axis=1)
