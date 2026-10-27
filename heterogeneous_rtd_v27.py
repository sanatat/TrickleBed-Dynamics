"""One-dimensional dynamic gas-liquid-solid TRICKLE BED TRACER comparison.
Sales et al. 2005 equations 7--9 motivate axial dispersion, sphere pore diffusion,
partial surface wetting and linear reversible adsorption; not their exact Laplace
transfer function (equations 13--15). Explicitly labelled assumptions below.

The numerical calculation is linear nonreactive; it is NOT the Alopaeus thermal model.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from pathlib import Path
import csv, json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.sparse import lil_matrix, csr_matrix
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'/'v27_heterogeneous_rtd'
OLD=ROOT/'results'/'v26_vector_extraction'

@dataclass(frozen=True)
class Params:
    length_m: float = .8                  # experimental reactor, Sales p447
    superficial_u_m_s: float = .001       # Sales Table 3, QL=5 L/h row
    liquid_holdup: float = .17            # Sales Table 3
    wetting_efficiency: float = .38       # Sales Table 3
    peclet_reported: float = 4.84        # Sales Table 3
    bed_porosity_assumed: float = .4      # not supplied explicitly
    particle_radius_m: float = .0025     # midpoint diameter 4--6 mm
    internal_porosity: float = .33        # Sales p447
    particle_bulk_density: float = 1100. # kg/m3, Sales p447, used for adsorption
    surface_kt_m_s: float = 8e-5         # bracketed by Sales Table 2 6.11--11.81 e-5
    effective_pore_D_m2_s: float = 7.65e-10  # Sales p450 text; Table 2 discrepant
    K_ads_m3_kg: float = 8.2e-5          # Sales p450 text
    adsorption_relaxation_s_inv: float = .03 # illustrative, UNPUBLISHED
    axial_dispersion_scale: float = 1.   # adjustable scale; reference=1
    interface_scale: float = 1.          # adjust interface exchange; reference=1
    accessible_pore_fraction: float = 1. # v2.8: uncertain share of wetted porous solid that participates in tracer exchange
    n_axial: int = 24
    n_radial: int = 3


def effective_dispersion_per_bed(p: Params):
    # Sales nomenclature Pe=Z*u/(D_ax*epsilon); this is the paper's
    # definition. Our FV uses the resulting D_ax directly in bed flux.
    return p.length_m*p.superficial_u_m_s / (p.peclet_reported*p.bed_porosity_assumed) * p.axial_dispersion_scale


def build_linear_system(p: Params, enable_pore=True):
    """d(x-xeq)/dt = A (x-xeq) + b*(Cin-1).
    States in each axial cell: fluid c, radial shell c[0..Nr-1], q[0..Nr-1].
    Surface exchange is one flux with opposite signs -> strict inventory balance.
    """
    nx, nr = p.n_axial, p.n_radial
    assert 0<p.liquid_holdup<p.bed_porosity_assumed<1
    assert 0<=p.wetting_efficiency<=1 and 0<p.accessible_pore_fraction<=1 and nr>=1 and nx>=2
    stride=1+2*nr
    size=nx*stride
    A=lil_matrix((size,size),dtype=float)
    b=np.zeros(size)
    dz=p.length_m/nx
    u=p.superficial_u_m_s
    h=p.liquid_holdup
    D=effective_dispersion_per_bed(p)
    f_solid=(1-p.bed_porosity_assumed)*p.wetting_efficiency*p.accessible_pore_fraction
    R=p.particle_radius_m
    edges=np.linspace(0.,R,nr+1)
    centers=((edges[:-1]**3+edges[1:]**3)/2)**(1/3)
    volfrac=np.diff(edges**3)/R**3
    A_outer=f_solid*3./R
    def liq(j): return j*stride
    def pore(j,k): return j*stride+1+k
    def ads(j,k): return j*stride+1+nr+k
    for j in range(nx):
        ix=liq(j)
        # conservative upwind advection: (Fleft-Fright)/(h*dz)
        A[ix,ix]-=u/(h*dz)
        if j>0: A[ix,liq(j-1)]+=u/(h*dz)
        else: b[ix]+=u/(h*dz)
        # non-equilibrium Danckwerts inlet; zero-gradient outlet.
        if j>0:
            A[ix,ix]-=D/(h*dz**2)
            A[ix,liq(j-1)]+=D/(h*dz**2)
        if j<nx-1:
            A[ix,ix]-=D/(h*dz**2)
            A[ix,liq(j+1)]+=D/(h*dz**2)
        if not enable_pore:
            continue
        k_transfer=A_outer * p.surface_kt_m_s*p.interface_scale
        shell_last=pore(j,nr-1)
        A[ix,ix]-=k_transfer/h
        A[ix,shell_last]+=k_transfer/h
        # shell accumulation: active pellet volume per bed * shell volume * eps_in
        for k in range(nr):
            ic=pore(j,k); iq=ads(j,k)
            vol=f_solid*volfrac[k]
            capacity=vol*p.internal_porosity
            if capacity<=0: continue
            # adsorption: pore inventory +/- adsorbate, q (kg tracer/kg catalyst)
            kads=p.adsorption_relaxation_s_inv
            A[iq,ic]+=kads*p.K_ads_m3_kg
            A[iq,iq]-=kads
            A[ic,ic]-=vol*p.particle_bulk_density*kads*p.K_ads_m3_kg/capacity
            A[ic,iq]+=vol*p.particle_bulk_density*kads/capacity
            if k==nr-1:
                A[ic,ix]+=k_transfer/capacity
                A[ic,ic]-=k_transfer/capacity
        for k in range(nr-1):
            a=pore(j,k);bidx=pore(j,k+1)
            rface=edges[k+1]
            conductance=f_solid*(3*rface**2/R**3)*p.effective_pore_D_m2_s/(centers[k+1]-centers[k])
            ca=f_solid*volfrac[k]*p.internal_porosity
            cb=f_solid*volfrac[k+1]*p.internal_porosity
            A[a,a]-=conductance/ca; A[a,bidx]+=conductance/ca
            A[bidx,bidx]-=conductance/cb; A[bidx,a]+=conductance/cb
    return csr_matrix(A),b


def inventory_weights(p: Params):
    """Inventory weights per unit reactor cross-sectional area per axial slice."""
    nx,nr=p.n_axial,p.n_radial
    stride=1+2*nr
    w=np.zeros(nx*stride)
    dz=p.length_m/nx
    edges=np.linspace(0,p.particle_radius_m,nr+1)
    shells=np.diff(edges**3)/p.particle_radius_m**3
    active=(1-p.bed_porosity_assumed)*p.wetting_efficiency*p.accessible_pore_fraction
    for j in range(nx):
        k=j*stride
        w[k]=p.liquid_holdup*dz
        for i in range(nr):
            vol=active*shells[i]*dz
            w[k+1+i]=vol*p.internal_porosity
            w[k+1+nr+i]=vol*p.particle_bulk_density
    return w


def check_global_inventory(p: Params):
    """No hidden source/sink except external inlet/outlet for any state."""
    A,b=build_linear_system(p)
    w=inventory_weights(p)
    outlet=np.zeros(A.shape[0]);outlet[(p.n_axial-1)*(1+2*p.n_radial)]=p.superficial_u_m_s
    residual=np.asarray(w@A).ravel()+outlet
    return float(np.max(np.abs(residual))),float(abs(w@b-p.superficial_u_m_s))


def tracer_response(times, inlet_t, inlet_c, p, enable_pore=True, rtol=2e-6, atol=1e-8):
    t=np.asarray(times,float)
    if np.any(np.diff(t)<0) or min(t)<0: raise ValueError('times must increase and be nonnegative')
    A,b=build_linear_system(p,enable_pore)
    inp_t=np.asarray(inlet_t,float); inp_c=np.clip(inlet_c,0.,1.)
    def fun(now,y):
        return A@y + b*(float(np.interp(now,inp_t,inp_c,left=1.,right=inp_c[-1]))-1.)
    s=solve_ivp(fun,(0.,max(t)),np.zeros(A.shape[0]),method='BDF',jac=A,
         rtol=rtol,atol=atol,t_eval=t)
    if not s.success: raise RuntimeError(s.message)
    outlet=1.+s.y[(p.n_axial-1)*(1+2*p.n_radial)]
    return outlet


def load_vector_curve(name):
    data=np.genfromtxt(OLD/f'{name}_vector.csv',delimiter=',',names=True)
    return np.asarray(data['time_s']),np.clip(data['normalized_signal'],0,1)


def metrics(x,y): return float(np.sqrt(np.mean((np.asarray(x)-np.asarray(y))**2)))


def simulate_and_report():
    OUT.mkdir(parents=True,exist_ok=True)
    t,cin=load_vector_curve('fig5_experimental_entrance')
    ty,obs=load_vector_curve('fig5_experimental_exit')
    if not np.allclose(t,ty):raise ValueError('misaligned Fig 5 vector series')
    p=Params()
    # Physical *conditional* predictions: all chosen constants printed and sourced.
    pred_no_pore=tracer_response(t,t,cin,p,enable_pore=False)
    pred_pore=tracer_response(t,t,cin,p,enable_pore=True)
    # Calibrate only one global dispersion scale, then two parameters to diagnose identifiability.
    def objective(logscale):
        pp=replace(p,axial_dispersion_scale=float(np.exp(logscale)))
        return np.mean((tracer_response(t,t,cin,pp)-obs)**2)
    # A 1-D bounded scalar minimization avoids tiny finite-difference gradients
    # swamped by ODE tolerances (a failure mode in previous v2.7 draft).
    opt=minimize_scalar(objective,bounds=(-np.log(12.),np.log(12.)),method='bounded',
         options={'xatol':.03,'maxiter':19})
    fitted_p=replace(p,axial_dispersion_scale=float(np.exp(opt.x)))
    pred_cal=tracer_response(t,t,cin,fitted_p)
    # Gamma curve already fitted in v2.6 to the exact same digitized figure.
    from sales2005_transfer_v21 import transfer
    gamma_metrics=json.loads((OLD/'MODEL_METRICS.json').read_text())['figure5']
    pred_gamma=transfer(t,t,cin,gamma_metrics['gamma_mean_s'],gamma_metrics['gamma_shape'],gamma_metrics['additional_delay_s'])
    # Variation in assumed epsilon and uncertain 200x discrepancy in published D_eff.
    variation=[]
    for bed_eps in (.35,.4,.45):
        for diff in (7.65e-10,1.8e-7):
            pp=replace(p,bed_porosity_assumed=bed_eps,effective_pore_D_m2_s=diff)
            yy=tracer_response(t,t,cin,pp)
            variation.append(dict(bed_porosity_assumed=bed_eps,pore_D_m2_s=diff,figure5_graphical_rmse=metrics(yy,obs)))
    # Diagnostic: mass conservation for independent closed fluid + surface submodel via A property.
    # Uniform concentration should be invariant if Cin uniform, i.e. A@ones+b=0.
    A,b=build_linear_system(p)
    cols=np.zeros(A.shape[0])
    stride=1+2*p.n_radial
    for j in range(p.n_axial):
        cols[j*stride]=1.0
        cols[j*stride+1:j*stride+1+p.n_radial]=1.0
        cols[j*stride+1+p.n_radial:(j+1)*stride]=p.K_ads_m3_kg
    fixed_point_max=float(np.max(np.abs(A@cols+b)))
    inventory_residual, inlet_flux_error=check_global_inventory(p)
    num={
      'source':'Sales et al. (2005), published Figure 5 vector-derived curves, not instrument data',
      'physical_model':'conservative 1-D FV axial dispersion + wetted spherical intrapellet radial Fick diffusion + reversible linear adsorption',
      'comparison':{
        'no_pore_reference_rmse':metrics(pred_no_pore,obs),
        'published_parameter_heterogeneous_rmse':metrics(pred_pore,obs),
        'one_parameter_dispersion_fit_rmse':metrics(pred_cal,obs),
        'gamma_three_parameter_fit_rmse':metrics(pred_gamma,obs),
        'dispersion_scale_fit':fitted_p.axial_dispersion_scale,
        'dispersion_scale_at_bound':bool(np.isclose(abs(opt.x),np.log(12),atol=.05)),
        'fit_success':bool(opt.success),
      },
      'numerics':{'nx':p.n_axial,'nr':p.n_radial,'dimension':A.shape[0], 'fixed_point_residual_max_s_inverse':fixed_point_max,
                  'global_inventory_balance_residual_m_s':inventory_residual,'global_inlet_flux_error_m_s':inlet_flux_error},
      'literature_and_assumptions':{
        'liquid_holdup_Table3':p.liquid_holdup,'wetting_Table3':p.wetting_efficiency,
        'Pe_Table3':p.peclet_reported,'surface_kt_Table2_bracketed_m_s':p.surface_kt_m_s,
        'Deff_p450_prose_m2_s':p.effective_pore_D_m2_s,
        'Deff_table2_contradiction_m2_s':1.8e-7,
        'Kads_p450_prose_m3_per_kg':p.K_ads_m3_kg,
        'bed_porosity_assumed':p.bed_porosity_assumed,'adsorption_rate_assumed_s_inv':p.adsorption_relaxation_s_inv,
        'particle_radius_midrange_m':p.particle_radius_m,
        'Dax_bed_flux_m2_s':effective_dispersion_per_bed(p),
      },'sensitivity':variation,
      'major_limitations':['Model motivated by Sales Eqs. 7-9, not exact analytical transfer function Eqs. 13-15.',
        'Rate of adsorption k_ads and bed porosity unreported in available paper and are assumptions.',
        'Figure 5 is one flow condition; fitting RMSE is not experimental generalization.',
        'Figure 6 inlet series missing; no independent external-test prediction.',
        'No reaction, gas-phase hydrodynamics, heat transfer or Alopaeus CFD validation.',
        'Dense vector-resampled figure samples are strongly correlated.']}
    with (OUT/'metrics.json').open('w') as f:json.dump(num,f,indent=2)
    with (OUT/'figure5_models.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s','inlet_vector','outlet_vector','axial_only','heterogeneous_fixed','heterogeneous_one_parameter_calibrated','gamma_three_parameter']);w.writerows(zip(t,cin,obs,pred_no_pore,pred_pore,pred_cal,pred_gamma))
    with (OUT/'parameter_sensitivity.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=variation[0].keys());w.writeheader();w.writerows(variation)
    fig,ax=plt.subplots(figsize=(10,5.8))
    ax.plot(t,cin,':',c='gray',lw=1.8,label='Published Fig. 5 inlet (digitized)')
    ax.plot(t,obs,'o',ms=3,c='black',label='Published Fig. 5 outlet (digitized)')
    ax.plot(t,pred_no_pore,'--',label='Axial dispersion alone (fixed)')
    ax.plot(t,pred_pore,'-',label='Heterogeneous: literature/assumed params')
    ax.plot(t,pred_cal,'-',label='Heterogeneous: 1 parameter fitted')
    ax.plot(t,pred_gamma,'-.',label='Gamma: 3 parameters fitted (v2.6)')
    ax.set(xlabel='Time (s)',ylabel='Normalized tracer concentration',xlim=(0,1800),ylim=(-.04,1.05),title='Sales (2005) Fig. 5 — models vs graph-derived observations');ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'figure5_heterogeneous_comparison.png',dpi=190);plt.close(fig)
    report=f'''# v2.7 — Heterogeneous RTD inspired by Sales et al. (2005)

## Scientific basis

User-supplied Sales et al. (2005), *Experimental Dynamic Evaluation of Three-Phase Reactors*, **printed pages 445–450**, esp. Eqs. **(7)–(9)** for external liquid transport, spherical intrapellet tracer diffusion, adsorption; Eq. (11)/(12) Danckwerts boundaries and internal sphere flux; Tables **2–3** and Figs. **5–6**. We do **not** claim a literal implementation of analytical transfer function Eqs. (13)–(15).

## Computation and comparison

Using Figure 5 vector-extracted graphical traces from v2.6 (QL = **5 L/h**, QG = **300 L/h**), simulate an initially equilibrated negative tracer step with measured inlet forcing, conservative upwind finite volume advection, axial diffusion, wetted external-sphere exchange, radial spherical diffusion and reversible linear adsorption. Analytical adsorption relation is replaced by a first-order relaxation dynamic qdot = k_ads*(K_ads*C-q), whose kinetic coefficient is **assumed**, not claimed from the paper.

| Variant | Number tuned parameters | Graphical RMSE |
|---|---:|---:|
| Axial-only fixed parameters | 0 | {metrics(pred_no_pore,obs):.6f} |
| Heterogeneous fixed literature/assumed | 0 | {metrics(pred_pore,obs):.6f} |
| Heterogeneous + one fitted axial dispersion multiplier | 1 | {metrics(pred_cal,obs):.6f} |
| Gamma RTD (v2.6) | 3 | {metrics(pred_gamma,obs):.6f} |

One-parameter fit: D_ax multiplier = **{fitted_p.axial_dispersion_scale:.5g}**; reached optimization bound: **{num['comparison']['dispersion_scale_at_bound']}**. Lower graph RMSE is not alone proof of improved physics. These are in-sample graphical fits and have unequal numbers of calibrated parameters.

Conservative stationary check max |A*x_equilibrium+b| = **{fixed_point_max:.3g} s^-1**. The **global inventory identity** w@A + outlet_flux = 0 yields maximum residual **{inventory_residual:.3g} m/s** and inlet residual **{inlet_flux_error:.3g} m/s** for any state. This certifies discrete tracer inventory bookkeeping, not the physical realism of the closures.

## Explicit provenance of parameters

From Table 3: hL=0.17, fe=0.38, Pe=4.84; uL=0.001 m/s, for the low gas-flow row. Table 2 gives k_LS **(6.11–11.81)e-5 m/s**, choose **8e-5 m/s** *within the range*, not explicitly for one Figure 5 run. Page 450 text states Deff = **7.65e-10 m²/s**, whereas **Table 2 displays a different range in 10^-7 units**; these should not be reconciled without original author clarification. We therefore perform a two-order-of-magnitude sensitivity check, not claim either is universally correct. Text gives Kads = 8.20e-5 m³/kg. The adsorption kinetics rate 0.03/s and packed-bed porosity 0.40 are **assumptions**. Other reported dimensions: 0.8m length, particle diameter 4–6mm; choose midpoint 5mm.

The paper's Pe definition Pe=Z*u/(D_ax*epsilon) is applied as an effective bulk flux coefficient using **assumed epsilon**, not as universal dispersion calibration. This discretization and boundary implementation differ from the authors' inverse Laplace method. Film area and interior liquid volume are taken proportional to fe, which is a simplifying closure rather than the article's exact microstructure treatment.

## Why Figure 6 is not claimed as independent validation

Published Figure 6 (QL = 10 L/h) has extracted **experimental outlet markers** but **no paired numerical inlet signal**. Choosing an unmeasured inlet waveform to maximize output agreement would be circular. Thus our v2.7 report is a rigorous **model comparison for Figure 5**, not a fully validated second-flow predictive test. Likewise it is **not** validation of the reactive or exothermic 2012 Alopaeus model.

## Sensitivity

See `parameter_sensitivity.csv`, which varies assumed bed porosity and the disputed Deff scale. The resulting changes are diagnostic *model uncertainty*, not confidence intervals.

## Files and run

- `heterogeneous_rtd_v27.py` — sparse FV physical model, figures/metrics generator.
- `figure5_models.csv` — vector-derived measured figure curves and predictions.
- `figure5_heterogeneous_comparison.png` — graphic.
- `metrics.json` — reproducible results, explicit assumptions and limitations.
- `parameter_sensitivity.csv` — physical assumption stress test.

```bash
python heterogeneous_rtd_v27.py
python -m pytest -q
```

No copied source PDF included in archive; raw source PDF is held separately by the user.
'''
    (OUT/'REPORT.md').write_text(report)
    print(json.dumps(num['comparison'],indent=2))
    print('fixed_point_residual:',fixed_point_max)
    return num

if __name__=='__main__':simulate_and_report()
