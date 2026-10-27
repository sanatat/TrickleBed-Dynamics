"""Transparent educational 1D reduced-order trickle-bed model.
Not a reproduction of Gorshkova et al. (2012), nor a validated industrial CFD model.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp

@dataclass(frozen=True)
class Config:
    n: int = 24
    length: float = 1.5             # m, paper test geometry
    bed_porosity: float = 0.37      # paper Table 2
    pellet_porosity: float = 0.5    # paper Table 2
    particle_diameter: float = 0.003 # m, paper Table 2
    pressure_bar: float = 20.0      # paper Table 2; not dynamically solved
    Tin: float = 400.0             # K, paper Table 2
    cAin: float = 1500.0           # mol/m3, scenario assumption
    cHin: float = 8.0              # mol/m3, scenario assumption
    hstar: float = 12.0            # mol/m3, assumed dissolved H2 saturation
    superficial_u: float = 0.0010  # m/s, assumption (NOT paper flow converted)
    kLa: float = 0.022             # 1/s, assumed gas-to-liquid transfer
    k_pore: float = 0.008          # 1/s, effective bulk/pore exchange timescale
    k_rxn: float = 0.014          # 1/s, normalized, assumed reaction kinetics
    H_half: float = 4.0           # mol/m3, assumed H saturation
    dH: float = -90000.0           # J/mol, assumed heat of reaction
    heat_capacity: float = 1.45e6 # J/(m3-bed K), effective heat capacity
    heat_loss: float = 220.0       # W/(m3-bed K), assumed
    tau_wet: float = 260.0        # s, assumed memory relaxation
    step_time: float = 200.0       # s
    t_end: float = 1200.0         # s
    flow_after: float = 0.7        # post-step ratio of liquid flow


def operation(t, c):
    flow_ratio = 1.0 if t < c.step_time else c.flow_after
    # An illustrative equilibrium-wetting correlation; NOT a published correlation.
    target_wet = np.clip(0.72 * flow_ratio**0.40, 0.05, 0.95)
    return flow_ratio, target_wet


def unpack(y, n):
    return [y[i*n:(i+1)*n] for i in range(8)]


def initial_state(c, wet0):
    n = c.n
    # bulk A,H,P; pore A,H,P; temperature; wetting fraction
    return np.concatenate([np.full(n,c.cAin), np.full(n,c.cHin), np.zeros(n),
                           np.full(n,c.cAin), np.full(n,c.cHin), np.zeros(n),
                           np.full(n,c.Tin), np.full(n,wet0)])


def rhs(t, y, c, reaction=True, gas_supply=True):
    n = c.n
    A,H,P,Ap,Hp,Pp,T,f = unpack(y,n)
    dz = c.length/n
    ratio, ftarget = operation(t,c)
    u = c.superficial_u * ratio
    # Fixed liquid holdup and effective pore volume; simplifying assumptions.
    vb = 0.25                 # m3 liquid per m3 bed, assumed
    vp = (1-c.bed_porosity)*c.pellet_porosity*0.65 # m3 wetted pore capacity / m3 bed
    # exchange mol/(m3-bed s); f affects exchange surface and reaction contact
    wet = np.clip(f,0.01,0.99)
    exchange = c.k_pore * vp * wet
    JA = exchange * (A-Ap)
    JH = exchange * (H-Hp)
    JP = exchange * (P-Pp)
    # Source mol/(m3-bed s) for A consumption and P formation.
    r = (c.k_rxn * vp * wet * np.maximum(Ap,0) *
         np.maximum(Hp,0)/(c.H_half+np.maximum(Hp,0))) if reaction else np.zeros(n)
    def adv(C, Cin):
        upstream = np.concatenate([[Cin],C[:-1]])
        return (u/dz)*(upstream-C)
    gas = c.kLa*vb*wet*(c.hstar-H) if gas_supply else np.zeros(n)
    dA = (adv(A,c.cAin)-JA)/vb
    dH = (adv(H,c.cHin)-JH+gas)/vb
    dP = (adv(P,0)+JP*(-1))/vb
    dAp = (JA-r)/vp
    dHp = (JH-r)/vp
    dPp = (JP+r)/vp
    # Thermally mixed pseudo-homogeneous energy balance, not full gas/liquid/solid temperatures.
    upstreamT = np.concatenate([[c.Tin],T[:-1]])
    adv_heat = c.heat_capacity * (u/dz) * (upstreamT-T)
    dT = (adv_heat + (-c.dH)*r - c.heat_loss*(T-c.Tin))/c.heat_capacity
    df = (ftarget-f)/c.tau_wet
    return np.concatenate([dA,dH,dP,dAp,dHp,dPp,dT,df])


def simulate(c=Config(), wet0=0.2, reaction=True, gas_supply=True):
    y0 = initial_state(c,wet0)
    # Force integration to resolve discontinuity at step_time.
    ts=np.linspace(0,c.t_end,241)
    pre = ts[ts<=c.step_time]
    post = ts[ts>=c.step_time]
    a=solve_ivp(rhs,(0,c.step_time),y0,args=(c,reaction,gas_supply),
                method='BDF',t_eval=pre,rtol=1e-6,atol=1e-8)
    if not a.success: raise RuntimeError(a.message)
    b=solve_ivp(rhs,(c.step_time,c.t_end),a.y[:,-1],args=(c,reaction,gas_supply),
                method='BDF',t_eval=post,rtol=1e-6,atol=1e-8)
    if not b.success: raise RuntimeError(b.message)
    return np.concatenate([a.t,b.t[1:]]),np.concatenate([a.y,b.y[:,1:]],axis=1)


def metrics(t,y,c):
    A,H,P,Ap,Hp,Pp,T,f = [a.reshape(c.n,-1) for a in unpack(y,c.n)]
    return {'time_s':t,'conversion_out':1-A[-1]/c.cAin,
            'temperature_peak_K':T.max(axis=0),'temperature_out_K':T[-1],
            'wetting_mean':f.mean(axis=0),'outlet_A_molm3':A[-1],
            'outlet_P_molm3':P[-1]}
