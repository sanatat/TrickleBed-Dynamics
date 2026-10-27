"""Conservative variable-pore-saturation 1-D study, synthetic parameters.

State: species MOLES per m³ of total bed for L,G,W,D (3 each), shared T,
wetting saturation s. The 2012 paper motivated separate inner inventories.
This is NOT its capillary/evaporation/CFD reproduction.
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.integrate import solve_ivp
from v4_model import Params as Base, inlets, rate

@dataclass(frozen=True)
class Params(Base):
    interface_area_mode: bool = True  # paper Eq. interface geometry; NOT hydrodynamic closure
    wet_surface_exponent: float = 1.0  # illustrative mapping pore saturation -> outer wet fraction
    sherwood_mode: bool = True
    Re_p: float = 2.0  # representative not evaluated from coupled fluid hydraulics
    Sc_p: float = 100.0  # representative not evaluated from local properties
    pore_vol: float = .315  # (1-void)*pellet_porosity; paper-supported geometry
    s_min: float = .04
    s_max: float = .96
    cap_tau: float = 140.
    # No phase-change model; constant idealized phase equilibrium ratios


def sherwood(Re, Sc):
    """Gorshkova et al. (2012), Eq. 7 (Ranz-Marshall form)."""
    if np.any(np.asarray(Re) < 0) or np.any(np.asarray(Sc) <= 0):
        raise ValueError('Re must be >=0 and Sc >0')
    return 2.0 + 0.6*np.sqrt(Re)*np.cbrt(Sc)

def interface_areas(s,p):
    """Eqs in Gorshkova 2012 p398, per bed volume m2/m3.

    fb and fs are ONLY proxies from s; published Lappalainen wetting closure
    is unavailable here. Interfacial areas are structural, not fitted fluxes.
    """
    s=np.clip(np.asarray(s),0,1)
    fb=np.clip(s**p.wet_surface_exponent,0,1)
    fs=np.clip(s,0,1)
    a=6*(1-p.void)/p.pellet_diameter
    # area GL=fb*a, GD=(1-fb)*epsilon_c*a, LW=fb*epsilon_c*a
    return {'GL':fb*a,'GD':(1-fb)*p.pellet_porosity*a,
            'LW':fb*p.pellet_porosity*a,
            'GW':6*p.pellet_porosity*(1-p.void)*np.maximum(fs-fb,0)/p.pellet_diameter,
            'DW':6*p.pellet_porosity*(1-p.void)*fs*(1-fs)/p.pellet_diameter}

def volumes(s,p):
    return np.array([np.full_like(s,p.V_L),np.full_like(s,p.V_G),p.pore_vol*s,p.pore_vol*(1-s)])

def unpack(y,p):
    a=np.asarray(y).reshape(14,p.n)
    return a[:12].reshape(4,3,p.n),a[12],a[13]

def initial(p):
    l,g=inlets(p); s=np.full(p.n,p.wet_init)
    V=volumes(s,p)
    C=np.stack([np.repeat(v[:,None],p.n,axis=1) for v in (l,g,l,g)])
    return np.concatenate([(V[:,None,:]*C).ravel(),np.full(p.n,p.T_in),s]).astype(float)

def rhs(t,y,p,flow=True,reaction=True,transfer=True):
    N,T,s=unpack(y,p)
    V=volumes(s,p)
    C=N/np.maximum(V[:,None,:],1e-12)
    L,G,W,D=C
    dN=np.zeros_like(N)
    dz=p.length/p.n
    f=(1. if t<p.t_step else p.flow_ratio) if flow else 0.
    li,ga=inlets(p)
    for idx,inc,u in [(0,li,p.u_L*f),(1,ga,p.u_G if flow else 0.)]:
        dN[idx]+=u/dz*(np.column_stack((inc,C[idx,:,:-1]))-C[idx])
    K=np.array([p.K_A,p.K_H,p.K_P])[:,None]
    if transfer:
        # Flux baseline parameters remain synthetic; area factor is relative to
        # a fully wetted pellet. Published geometric interface areas replace
        # the previous ad-hoc scalar s factors (when enabled).
        if p.interface_area_mode:
            a=interface_areas(s,p)
            abase=6*(1-p.void)/p.pellet_diameter
            GL=a['GL']/abase; GD=a['GD']/(p.pellet_porosity*abase)
            LW=a['LW']/(p.pellet_porosity*abase)
            DW=a['DW']/(p.pellet_porosity*abase)
        else:
            GL=np.ones_like(s);GD=1-s;LW=s;DW=s
        # Sherwood scaling is dimensionless, referenced to input Re,Sc.
        # Kinetic constants still NOT independently measured.
        sh_factor=(sherwood(p.Re_p,p.Sc_p)/sherwood(2.0,100.0)) if p.sherwood_mode else 1.0
        J=[(0,1,p.k_GL*sh_factor*GL[None,:]*(G-K*L)),
           (1,3,p.k_GD*GD[None,:]*(G-D)),
           (0,2,p.k_LW*p.chi*sh_factor*LW[None,:]*(L-W)),
           (3,2,p.k_DW*DW[None,:]*(D-K*W))]
        dN[0]+=J[0][2];dN[1]-=J[0][2]
        for source,dest,j in J[1:]:
            dN[source]-=j;dN[dest]+=j
    # Prescribed wetting trajectory, not a measured or published hysteresis closure
    target=p.wet_adv if f>=1 else p.wet_rec
    target=np.clip(target,p.s_min,p.s_max)
    ds=(target-s)/p.cap_tau
    # Volume carried into (or out of) pores: species flux equals volume flux * donor concentration.
    # Compensation of gas displacement D->G; external liquid and gas holdups
    # remain fixed: this is inventory-conservative, NOT volume-consistent hydrodynamics.
    dv=p.pore_vol*ds
    ingress=np.maximum(dv,0.)
    egress=np.maximum(-dv,0.)
    JL=ingress[None,:]*L-egress[None,:]*W   # L->W (+)
    JD=ingress[None,:]*D-egress[None,:]*G   # D->G (+)
    dN[0]-=JL;dN[2]+=JL
    dN[3]-=JD;dN[1]+=JD
    r=V[2]*s*rate(W[0],W[1],T,p) if reaction else np.zeros(p.n)
    dN[2,0]-=r;dN[2,1]-=r;dN[2,2]+=r
    # Heat equation: single pseudo-mixture, capillary latent/sensible effects omitted
    dh=(p.u_L*f+(p.u_G if flow else 0.))*.5/dz*(np.r_[p.T_in,T[:-1]]-T)
    dT=dh+(-p.dHr*r-p.heatloss*(T-p.T_in))/p.heatcap
    return np.r_[dN.ravel(),dT,ds]

def simulate(p,flow=True,reaction=True,transfer=True):
    if not np.isclose(p.V_L+p.V_G,p.void):raise ValueError('bulk void mismatch')
    if not np.isclose(p.pore_vol,(1-p.void)*p.pellet_porosity):raise ValueError('pore volume mismatch')
    if not (p.s_min<=p.wet_init<=p.s_max and p.s_min<p.wet_rec<p.s_max and p.s_min<p.wet_adv<p.s_max):raise ValueError('saturation bounds')
    y0=initial(p)
    def sol(t0,t1,y,grid):
        z=solve_ivp(rhs,(t0,t1),y,args=(p,flow,reaction,transfer),method='BDF',t_eval=grid,rtol=1e-7,atol=1e-9)
        if not z.success:raise RuntimeError(z.message)
        return z
    grid=np.linspace(0,p.t_end,121)
    if flow and 0<p.t_step<p.t_end:
        a=sol(0,p.t_step,y0,np.r_[grid[grid<p.t_step],p.t_step])
        b=sol(p.t_step,p.t_end,a.y[:,-1],np.r_[p.t_step,grid[grid>p.t_step]])
        return np.r_[a.t,b.t[1:]],np.c_[a.y,b.y[:,1:]]
    z=sol(0,p.t_end,y0,grid)
    return z.t,z.y

def invariants(y,p):
    N,_,_=unpack(y,p)
    q=N.sum(axis=(0,2))
    return q[0]+q[2],q[1]+q[2]
