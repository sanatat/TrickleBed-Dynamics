"""Synthetic RTD proof-of-concept: gamma tanks-in-series impulse and hold-out test.

NOT experimental calibration and NOT 2012 CFD validation.
"""
from dataclasses import dataclass
import numpy as np
from scipy.stats import gamma
from scipy.optimize import least_squares

@dataclass(frozen=True)
class RTDParams:
    mean_residence_s: float
    effective_tanks: float
    delay_s: float

def impulse_response(t_s, p: RTDParams, relative_flow: float = 1.0):
    """E(t) [1/s] for an impulse; integrate over 0..infinity gives one.

    Flow-only extrapolation assumption: mean residence varies inversely with
    flow while effective tanks and dead-time stay constant.
    """
    t = np.asarray(t_s, dtype=float)
    if p.mean_residence_s <= 0 or p.effective_tanks <= 0 or p.delay_s < 0 or relative_flow <= 0:
        raise ValueError('Positive residence time, tanks, flow and nonnegative delay required')
    tau = p.mean_residence_s / relative_flow
    dt = t - p.delay_s
    return np.where(dt > 0, gamma.pdf(np.maximum(dt,0), a=p.effective_tanks, scale=tau/p.effective_tanks), 0.0)

def fit_rtd(time, observed, relative_flow=1.0, tanks_fixed=None, initial=None):
    t=np.asarray(time,float)
    y=np.asarray(observed,float)
    if initial is None: initial=[82.0, 4.0, 5.0]
    if tanks_fixed is None:
        def fun(x): return impulse_response(t, RTDParams(*x), relative_flow)-y
        res=least_squares(fun, x0=initial, bounds=([15.,1.,0.],[220.,20.,30.]),xtol=1e-11,ftol=1e-11,gtol=1e-11,max_nfev=250)
        p=RTDParams(*res.x)
    else:
        def fun(x): return impulse_response(t, RTDParams(x[0],tanks_fixed,x[1]), relative_flow)-y
        res=least_squares(fun,x0=[initial[0],initial[2]],bounds=([15.,0.],[220.,30.]),max_nfev=250)
        p=RTDParams(res.x[0],tanks_fixed,res.x[1])
    return p,res

def rmse(y,pred):
    return float(np.sqrt(np.mean((np.asarray(y)-np.asarray(pred))**2)))

def make_synthetic_data(seed=20261026):
    rng=np.random.default_rng(seed)
    t=np.arange(0.,361.,3.)
    truth=RTDParams(85.,5.2,7.)
    noise_sigma=0.00023
    data={}
    for label,q in [('calibration',1.),('holdout',1.35)]:
        clean=impulse_response(t,truth,q)
        measured=np.clip(clean+rng.normal(0,noise_sigma,len(t)),0,None)
        data[label]={'time_s':t,'relative_flow':q,'clean':clean,'observed':measured}
    return truth, data, noise_sigma

def residual_bootstrap(t, observed, fitted, runs=80, seed=20261027):
    rng=np.random.default_rng(seed)
    base=impulse_response(t,fitted)
    residual=np.asarray(observed)-base
    params=[]
    for _ in range(runs):
        y=np.clip(base+rng.choice(residual,size=len(residual),replace=True),0,None)
        p,_=fit_rtd(t,y, initial=[fitted.mean_residence_s,fitted.effective_tanks,fitted.delay_s])
        params.append([p.mean_residence_s,p.effective_tanks,p.delay_s])
    return np.asarray(params)
