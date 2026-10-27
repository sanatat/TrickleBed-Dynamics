"""Passive dissolved tracer coupled to the v6 axial wetting/holdup trajectory.
Tracer is NOT coupled back to the chemistry: one-way transport coupling only.
All data generated here are synthetic, not published tracer measurements.
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from v6_model import Params, simulate, unpack, volumes

@dataclass(frozen=True)
class TracerConfig:
    k_exchange: float = 0.0018  # m3_liquid / m3_bed /s parameterized exchange
    dispersion: float = 2e-5   # m2/s, effective axial dispersion
    pulse_time: float = 25.0
    pulse_sigma: float = 5.0

def pulse(t, cfg):
    return np.exp(-.5*((t-cfg.pulse_time)/cfg.pulse_sigma)**2)/(cfg.pulse_sigma*np.sqrt(2*np.pi))

def wetting_schedule(p, duration):
    # Same phenomenological wetting state equation as v6_model.rhs.
    # Treated as one-way hydrodynamic driver; reaction and tracer do not modify wetting.
    def dw(t, w):
        target = p.wet_adv if t < p.t_step else p.wet_rec
        return (target-w)/p.cap_tau
    sol=solve_ivp(dw,(0,duration),np.full(p.n,p.wet_init),dense_output=True,
                  method='RK45',rtol=1e-10,atol=1e-12)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.sol

def tracer_rhs(t,y,p,cfg,wet_fun, inlet_on=True):
    n=p.n; dz=p.length/n
    ML=y[:n]; MW=y[n:]
    s=np.clip(wet_fun(t),p.s_min,p.s_max)
    vl=p.V_L; vw=p.pore_vol*s
    L=ML/vl; W=MW/vw
    u=p.u_L*(1. if t<p.t_step else p.flow_ratio)
    ci=pulse(t,cfg) if inlet_on else 0.
    # Positive interface flux into pore; conserved across L and W.
    exch=cfg.k_exchange*s*(L-W)
    ds=(p.wet_adv-s)/p.cap_tau if t<p.t_step else (p.wet_rec-s)/p.cap_tau
    dv=p.pore_vol*ds
    pore_fill=np.maximum(dv,0)*L + np.minimum(dv,0)*W
    # Conservative finite-volume axial face flux, Danckwerts-like inlet with pure advective input.
    F=np.empty(n+1); F[0]=u*ci
    F[1:]=u*L
    if cfg.dispersion>0:
        F[1:-1]-=cfg.dispersion*(L[1:]-L[:-1])/dz
    dL=(F[:-1]-F[1:])/dz-exch-pore_fill
    dW=exch+pore_fill
    return np.r_[dL,dW]

def run_tracer(p,cfg,times,inlet_on=True):
    duration=float(max(times)); fun=wetting_schedule(p,duration)
    sol=solve_ivp(tracer_rhs,(0,duration),np.zeros(2*p.n),args=(p,cfg,fun,inlet_on),
                  t_eval=np.asarray(times),method='BDF',rtol=2e-7,atol=1e-10)
    if not sol.success: raise RuntimeError(sol.message)
    L=sol.y[:p.n]/p.V_L
    u=np.where(sol.t < p.t_step,p.u_L,p.u_L*p.flow_ratio)
    return {'time':sol.t,'outlet_conc':L[-1], 'outlet_flux':u*L[-1],
            'mobile_inventory':sol.y[:p.n].sum(axis=0)*p.length/p.n,
            'pore_inventory':sol.y[p.n:].sum(axis=0)*p.length/p.n}

def fit_tracer(p,t,measured,initial=(.001,1e-5)):
    def residual(x):
        cfg=TracerConfig(k_exchange=x[0],dispersion=x[1]);
        return (run_tracer(p,cfg,t)['outlet_conc']-measured)
    result=least_squares(residual, initial,bounds=([1e-6,1e-7],[.012,2e-4]),max_nfev=30,
                         xtol=2e-5,ftol=2e-5,gtol=2e-5)
    return result
