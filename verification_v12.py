"""Reproducible integrator-tolerance, conservation, and energy-closure diagnostics.
Not external validation. All kinetic and transport parameters remain illustrative.
"""
from dataclasses import replace
from pathlib import Path
import json,csv
import numpy as np
from scipy.integrate import solve_ivp
from v6_model import Params, initial, rhs, unpack, volumes


def integrate(p, rtol=1e-7, atol=1e-9, n_output=101):
    y=initial(p)
    times=np.linspace(0,p.t_end,n_output)
    def segment(a,b,y,ts):
        sol=solve_ivp(lambda t,x:rhs(t,x,p), (a,b), y, method='BDF',
                      t_eval=ts, rtol=rtol,atol=atol)
        if not sol.success: raise RuntimeError(sol.message)
        return sol
    if 0<p.t_step<p.t_end:
        ta=np.unique(np.r_[times[times<p.t_step],p.t_step]);
        a=segment(0,p.t_step,y,ta)
        tb=np.unique(np.r_[p.t_step,times[times>p.t_step]])
        b=segment(p.t_step,p.t_end,a.y[:,-1],tb)
        return np.r_[a.t,b.t[1:]],np.c_[a.y,b.y[:,1:]]
    s=segment(0,p.t_end,y,times)
    return s.t,s.y

def total_equivalents(y,p):
    N,_,_=unpack(y,p)
    return np.array([np.sum(N[:,0,:]+N[:,2,:]),np.sum(N[:,1,:]+N[:,2,:])])

def budget_derivative(t,y,p):
    """Compare equivalent inventory derivatives against explicit inlet/outlet fluxes.
    Each equivalent is summed in all four compartments; spatial integration is dz.
    """
    N,T,S=unpack(y,p)
    C=N/volumes(S,p)[:,None,:]
    li=np.array([p.li_A,p.li_H,0.]);ga=np.array([p.ga_A,p.ga_H,0.])
    f=1. if t<p.t_step else p.flow_ratio
    # Advection driven boundary fluxes only (other sources conserve equivalents)
    flux_L=p.u_L*f*(li-C[0,:,-1]);flux_G=p.u_G*(ga-C[1,:,-1]);flux=flux_L+flux_G
    net=np.array([flux[0]+flux[2],flux[1]+flux[2]])
    dN,_,_=unpack(rhs(t,y,p),p)
    actual=np.array([np.sum(dN[:,0,:]+dN[:,2,:]),np.sum(dN[:,1,:]+dN[:,2,:])])*p.length/p.n
    return actual,net

def energy_residual(t,y,p):
    N,T,S=unpack(y,p)
    W=N[2]/volumes(S,p)[2][None,:]
    from v4_model import rate
    r=volumes(S,p)[2]*S*rate(W[0],W[1],T,p)
    f=1. if t<p.t_step else p.flow_ratio
    dz=p.length/p.n
    adv=(p.u_L*f+p.u_G)*.5/dz*(np.r_[p.T_in,T[:-1]]-T)
    predicted=adv+(-p.dHr*r-p.heatloss*(T-p.T_in))/p.heatcap
    _,actual,_=unpack(rhs(t,y,p),p)
    return float(np.max(np.abs(actual-predicted)))

def run(output=Path('results/v12_verification')):
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    p=Params(n=8,t_end=130.,t_step=60.,wet_init=.2)
    rows=[];solutions={}
    for rtol,atol in ((1e-5,1e-7),(1e-7,1e-9),(1e-9,1e-11)):
        t,Y=integrate(p,rtol,atol,81)
        N,T,S=unpack(Y[:,-1],p)
        outC=float((N/volumes(S,p)[:,None,:])[0,0,-1])
        record={'rtol':rtol,'atol':atol,'outlet_alkene_mol_m3':outC,
                'peak_T_K':float(T.max()),'min_inventory':float(N.min()),
                'sat_min':float(S.min()),'sat_max':float(S.max())}
        rows.append(record);solutions[str(rtol)]=Y
    tight=rows[-1]
    for x in rows:x['relative_outlet_vs_tight']=abs(x['outlet_alkene_mol_m3']-tight['outlet_alkene_mol_m3'])/max(abs(tight['outlet_alkene_mol_m3']),1e-12)
    # Perform local conservation checks for representative state before and after flow switch.
    checkpoints=[]
    for tm in (30.,90.):
        tt,_=integrate(p,1e-9,1e-11,81)
        idx=int(np.argmin(np.abs(tt-tm)))
        tm=float(tt[idx])
        y=solutions['1e-09'][:,idx]
        derivative,boundary=budget_derivative(tm,y,p)
        checkpoints.append({'t_s':tm,'carbon_derivative':float(derivative[0]),'carbon_boundary':float(boundary[0]),
          'hydrogen_derivative':float(derivative[1]),'hydrogen_boundary':float(boundary[1]),
          'max_inventory_flux_residual':float(np.max(np.abs(derivative-boundary))),
          'energy_rhs_residual_K_s':energy_residual(tm,y,p)})
    with (out/'tolerance_study.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    audit={'status':'numerical verification ONLY; synthetic uncalibrated 1D model',
           'tolerances':rows,'budget_checkpoints':checkpoints,
           'max_budget_residual':max(c['max_inventory_flux_residual'] for c in checkpoints),
           'warnings':['Energy RHS identity is a code consistency check, not a physical energy conservation proof.',
           'Wetting closure, pseudo-mixture heat equation, pressure/phase changes remain incomplete.',
           'Parameter values do not reproduce published 2012 CFD outputs.']}
    (out/'audit.json').write_text(json.dumps(audit,indent=2))
    (out/'VERIFICATION_REPORT.md').write_text(f'''# v1.2 targeted numerical verification\n\n**Model:** Uncalibrated one-dimensional four-compartment finite-volume prototype. NOT experimentally validated.\n\n## Improvements\n\n- Three independent BDF relative/absolute tolerance settings, explicitly splitting integration at the imposed feed step.\n- All-compartment carbon- and hydrogen-equivalent differential balances checked against external boundary flux.\n- Explicit energy RHS identity and diagnostics for positivity and pore-saturation bounds.\n\n## Results\n\n- Maximum equivalent flux budget residual: **{audit['max_budget_residual']:.3e}** (concentration×velocity units; per unit frontal area).\n- Max relative outlet variation vs tightest tolerance: **{max(r['relative_outlet_vs_tight'] for r in rows):.3e}**.\n- Lowest species inventory encountered at final state among runs: **{min(r['min_inventory'] for r in rows):.4g}** mol per m³ bed.\n\nThe finite-volume inventory flux identity checks implementation consistency; it is not a fit to physical experiment. Energy RHS comparison only checks the written equation against its derivative, not the accuracy of thermal constitutive laws. The physical model still lacks heat capacity variation, latent heats, pressure closure, genuine wetting correlations, two-isomer kinetics, measured calibration, and 2012 Figure 2 quantitative match.\n''')
    return audit
if __name__=='__main__':
    result=run()
    print(json.dumps({'max_budget_residual':result['max_budget_residual'], 'tolerances':result['tolerances']},indent=2))
