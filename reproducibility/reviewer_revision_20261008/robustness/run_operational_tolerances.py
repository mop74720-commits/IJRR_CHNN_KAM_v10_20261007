"""Reproduce robust full-FE initialization, port-error and work-governor results.

All inputs are relative to this file. Requirements: numpy, scipy, numba,
threadpoolctl and python-flint. Frozen nodal M,K define the physical plant;
the archived CHNN is used only to request force/field, so governor safety
does not require exact modal diagonalization or exact torus initialization.
"""
from pathlib import Path
from fractions import Fraction
import argparse, hashlib, importlib.metadata, json, os, sys, time
HERE = Path(__file__).resolve().parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / 'numba_cache'))
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE / 'inputs'))
import numpy as np
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits
from flint import arb, arb_mat, fmpq, ctx
from feature_model import energy_gradient, energy_hessian
threadpool_limits(limits=1)
ctx.prec = 256
MU, E0, ALPHA, LENGTH = 2e-5, 6e-7, .55558524, .08
CLEARANCE, DURATION, INITIAL_ERROR_Q = .004, 4., 5e-12
DELTA_THETA, DELTA_OMEGA, DELTA_V = .001, .01, .0001
RELATIVE_F, RELATIVE_B = .01, .01
ADDITIVE_F, ADDITIVE_B = .5e-9, .05e-6
FMAX = (20e-6-ADDITIVE_F)/(1+RELATIVE_F)
BMAX = (.02-ADDITIVE_B)/(1+RELATIVE_B)
GAIN = 5.
DATA = np.load(HERE/'inputs/model_matrices.npz')
FIT = np.load(HERE/'inputs/integral_fit.npz')
M, K, PHI = [DATA[k] for k in ['M', 'K', 'all_phi']]
OMEGA = DATA['all_omega']; W0 = float(DATA['omega0_rad_s'])
S = np.sqrt(E0/(W0*OMEGA)); N = 64
BT, BW = PHI[-1], PHI[-2]
KB = K.copy(); KB[-1,-1] += E0
MT = PHI.T@M@PHI; KBT = PHI.T@KB@PHI
MINV = np.linalg.inv(MT); DYNAMIC_K = MINV@KBT
FB, TB = MINV@BW, MINV@BT
AUTHORITY_INV = np.linalg.inv(np.column_stack([BW[:2], BT[:2]]))
COEFF, COSINE, DIREC, NU = (FIT['polynomial_coefficients'],
 float(FIT['cosine_multiplier']), FIT['d'], FIT['nu'])

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def exact(x):
    if isinstance(x, str): return arb(fmpq(str(Fraction(x))))
    return arb(fmpq(str(Fraction.from_float(float(x)))))
def matrix(a): return arb_mat([[exact(v) for v in row] for row in a])
def upper(x): return float(np.nextafter(float(x.abs_upper()), np.inf))
def lower(x): return float(np.nextafter(float(x.lower()), -np.inf))
def save_json(path, obj):
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\r\n')

def learned(z):
    h, g = energy_gradient(z, COEFF, COSINE, DIREC)
    return h, g, energy_hessian(z, COEFF, COSINE, DIREC)

def reference_initial():
    z = np.load(HERE/'inputs/finite_fourier_initial_state.npz')['state']
    _, grad, _ = learned(z)
    a, v = np.zeros(N), np.zeros(N)
    a[:2] = S[:2]*z[:2]; v[:2] = S[:2]*W0*grad[2:]
    return a, v

def bernstein_rows():
    rows=[]; le=LENGTH/32
    for el in range(32):
        b=np.zeros((4,66))
        b[0,2*el]=1; b[1,2*el]=1; b[1,2*el+1]=le/3
        b[2,2*el+2]=1; b[2,2*el+3]=-le/3; b[3,2*el+2]=1
        rows.extend(b[:,2:])
    return np.array(rows)
ROWS=bernstein_rows()

def certify():
    """Outward ball calculations for the exact exported nodal matrices."""
    am, ak = matrix(M), matrix(K)
    assert np.array_equal(M,M.T) and np.array_equal(K,K.T)
    mi, ki = am.inv(), ak.inv()
    # Ball LDL elimination proves positive definiteness of the frozen forms.
    for form in [am, ak]:
        work=arb_mat(form)
        for pivot in range(N):
            assert work[pivot,pivot] > 0
            for i in range(pivot+1,N):
                for j in range(i,N):
                    value=work[j,i]-work[j,pivot]*work[i,pivot]/work[pivot,pivot]
                    work[j,i]=value;work[i,j]=value
    # Exact Hermite Bernstein slope coefficients use length/32/3=1/1200.
    rb=arb_mat(len(ROWS),N)
    for i,row in enumerate(ROWS):
        for j,val in enumerate(row):
            if val: rb[i,j]=exact('1' if val==1 else ('1/1200' if val>0 else '-1/1200'))
    bc=rb*ki*rb.transpose()
    cw=max(bc[j,j].abs_upper() for j in range(len(ROWS)))
    ctheta=ki[N-1,N-1].abs_upper()
    ec=(exact('1/250')**2/(2*cw)).lower(); cap=exact('0.00000002266')
    assert cap<ec
    mf, mt = mi[N-2,N-2].abs_upper(), mi[N-1,N-1].abs_upper()
    residual_force=exact(str(ADDITIVE_F))*mf.sqrt()+exact(str(MU))*exact(str(ADDITIVE_B))*mt.sqrt()
    power_reserve=(residual_force**2/(4*exact(str(ALPHA)))).upper()
    reserve=(exact(str(DURATION))*power_reserve).upper()
    reference=np.load(HERE/'inputs/reference_nodal_state.npz')
    qb=matrix(reference['position_m_rad'].reshape(-1,1)); vb=matrix(reference['velocity_m_rad_per_s'].reshape(-1,1))
    kb=matrix(K);kb[N-1,N-1]+=exact(str(E0))
    eq=(vb.transpose()*am*vb+qb.transpose()*kb*qb)[0,0]/2
    errorq=exact(str(INITIAL_ERROR_Q))
    initial_bound=((eq.sqrt()+errorq.sqrt())**2).upper()
    command_cap=(cap-reserve).lower()
    assert initial_bound < command_cap
    kb_inv=kb.inv()
    error_tip=(2*kb_inv[N-2,N-2]*errorq).sqrt().abs_upper()
    error_slope=(2*kb_inv[N-1,N-1]*errorq).sqrt().abs_upper()
    rms_velocity=(2*errorq/(exact('1100')*arb.pi()*exact('0.0006')**2*exact('0.08'))).sqrt().abs_upper()
    out={'precision_bits':ctx.prec,'duration_s':DURATION,'initial_error_Q_J':INITIAL_ERROR_Q,
         'compliance_centerline_m_per_N_upper':upper(cw),'slope_compliance_upper':upper(ctheta),
         'critical_energy_J_lower':lower(ec),'energy_cap_J':lower(cap),
         'reference_quadratic_energy_J':upper(eq),'initial_energy_bound_J':upper(initial_bound),
         'command_work_cap_J':lower(command_cap),'additive_power_reserve_W':upper(power_reserve),
         'additive_work_reserve_J':upper(reserve),'mass_inverse_force_port_upper':upper(mf),
         'mass_inverse_rotation_port_upper':upper(mt),'shape_bound_m':upper((2*cw*cap).sqrt()),
         'tip_slope_bound_rad':upper((2*ctheta*cap).sqrt()),
         'error_set_implied_tip_m_upper':upper(error_tip),'error_set_implied_tip_angle_rad_upper':upper(error_slope),
         'error_set_implied_mass_weighted_rms_velocity_m_s_upper':upper(rms_velocity),
         'sensor_angle_error_rad':DELTA_THETA,'sensor_angular_rate_error_rad_s':DELTA_OMEGA,
         'sensor_tip_velocity_error_m_s':DELTA_V,'relative_force_error':RELATIVE_F,
         'relative_field_error':RELATIVE_B,'additive_force_error_N':ADDITIVE_F,
         'additive_field_error_T':ADDITIVE_B,
         'strict_positive_initial_work_margin_J':lower(command_cap-initial_bound),
         'balls':{k:str(v)for k,v in [('Cw_upper',cw),('Ec_lower',ec),('Ecap',cap),('Q_reference',eq),
          ('W0_upper',initial_bound),('additive_reserve_upper',reserve),('command_cap_lower',command_cap)]},
         'initialization_scope':'Full-state energy ellipsoid about finite exported Fourier reference; no exact-torus membership required.'}
    save_json(HERE/'analytic_certificate.json',out)
    return out

def invert_legendre(x,v):
    goal=v/(S[:2]*W0); y=goal/NU
    for _ in range(8):
        z=np.r_[x,y];_,g,h=learned(z); dy=np.linalg.solve(h[2:,2:],g[2:]-goal);y-=dy
        if max(abs(dy))<2e-14: return np.r_[x,y]
    raise RuntimeError('Legendre inverse failed')

# Eight scalar arguments control all nonlinear Jacobian dependence.
HMAP=np.zeros((8,2*N+3))
HMAP[0,0]=1/S[0];HMAP[1,1]=1/S[1]
HMAP[2,N]=1;HMAP[3,N+1]=1
HMAP[4,:N]=BT;HMAP[5,N:2*N]=BT;HMAP[6,N:2*N]=BW;HMAP[7,2*N]=1

def port_control(t,s,cap_normalized,governed=True):
    x,vl=s[:2],s[2:4]
    theta,angular,velocity,wn=s[4:]
    z=invert_legendre(x,vl)
    # Bounded errors in the nominal low-mode request are allowed; the work
    # theorem places no accuracy requirement on this requesting controller.
    z=z+2e-4*np.sin(np.array([3.,5.,7.,11.])*t+.3)
    _,g,h=learned(z)
    acc=S[:2]*W0**2*(h[2:,:2]@g[2:]-h[2:,2:]@g[:2])
    F,U=AUTHORITY_INV@(acc+ALPHA*vl+OMEGA[:2]**2*(S[:2]*x))
    B=(U-E0*(theta-np.sin(theta)))/(MU*np.cos(theta))
    request=np.clip([F,B],[-FMAX,-BMAX],[FMAX,BMAX])
    th=theta+DELTA_THETA*np.sin(13*t)
    om=angular+DELTA_OMEGA*np.cos(17*t)
    ve=velocity+DELTA_V*np.sin(11*t)
    ahat=np.array([ve,MU*np.cos(th)*om])
    delta=np.array([DELTA_V,MU*(DELTA_OMEGA+abs(om)*DELTA_THETA)])
    robust=delta+np.array([RELATIVE_F,RELATIVE_B])*(abs(ahat)+delta)
    preq=float(ahat@request+robust@abs(request))
    limit=GAIN*max(cap_normalized-wn,0)*E0
    factor=min(1.,limit/preq) if governed and preq>limit else 1.
    command=request*factor
    pbound=float(ahat@command+robust@abs(command))
    actual=command*(1+np.array([RELATIVE_F*np.sin(7*t),RELATIVE_B*np.cos(5*t)]))
    actual+=np.array([ADDITIVE_F*np.cos(19*t),ADDITIVE_B*np.sin(23*t)])
    moment=E0*(theta-np.sin(theta))+MU*actual[1]*np.cos(theta)
    actual_power=actual[0]*velocity+MU*actual[1]*np.cos(theta)*angular
    return np.array([actual[0],moment,pbound/E0,actual_power/E0]),command,actual,factor

def run(cert,refined=False):
    cap=cert['command_work_cap_J']/E0
    reference=np.load(HERE/'inputs/reference_nodal_state.npz')
    qref=reference['position_m_rad'];vref=reference['velocity_m_rad_per_s']
    dq=np.linalg.solve(KB,np.eye(N)[:,-2]);dq*=20e-6/dq[-2]
    dv=PHI[:,0]*(.00015/PHI[-2,0])
    qe=.5*(dq@KB@dq+dv@M@dv)
    dq*=np.sqrt(.99*INITIAL_ERROR_Q/qe);dv*=np.sqrt(.99*INITIAL_ERROR_Q/qe)
    a0=np.linalg.solve(PHI,qref+dq);v0=np.linalg.solve(PHI,vref+dv)
    # Verify the actual initialized error, including coordinate-solve rounding.
    eq=matrix(PHI)*matrix(a0.reshape(-1,1))-matrix(qref.reshape(-1,1))
    ev=matrix(PHI)*matrix(v0.reshape(-1,1))-matrix(vref.reshape(-1,1))
    kk=matrix(K);kk[N-1,N-1]+=exact(str(E0))
    qe_ball=(eq.transpose()*kk*eq+ev.transpose()*matrix(M)*ev)[0,0]/2
    assert qe_ball<exact(str(INITIAL_ERROR_Q))
    y0=np.r_[a0,v0,cert['initial_energy_bound_J']/E0,0.,0.]
    def rhs(t,y):
        a,v=y[:N],y[N:2*N]
        port,_,_,_=port_control(t,HMAP@y,cap)
        acceleration=-DYNAMIC_K@a-ALPHA*v+FB*port[0]+TB*port[1]
        return np.r_[v,acceleration,port[2],port[3],ALPHA*(v@MT@v)/E0]
    basejac=np.zeros((2*N+3,2*N+3));basejac[:N,N:2*N]=np.eye(N)
    basejac[N:2*N,:N]=-DYNAMIC_K;basejac[N:2*N,N:2*N]=-ALPHA*np.eye(N)
    def jac(t,y):
        s=HMAP@y;deriv=np.empty((4,8))
        for k in range(8):
            step=1e-6*max(abs(s[k]),1e-4 if k in [2,3,6] else .01)
            sp=s.copy();sm=s.copy();sp[k]+=step;sm[k]-=step
            deriv[:,k]=(port_control(t,sp,cap)[0]-port_control(t,sm,cap)[0])/(2*step)
        d=deriv@HMAP;j=basejac.copy()
        j[N:2*N]+=FB[:,None]*d[0]+TB[:,None]*d[1]
        j[2*N]=d[2];j[2*N+1]=d[3];j[2*N+2,N:2*N]=2*ALPHA*(MT@y[N:2*N])/E0
        return j
    times=np.linspace(0,DURATION,2001)
    atol=np.r_[np.full(N,2e-13),np.full(N,2e-11),np.full(3,2e-13)]
    scale=.1 if refined else 1.
    start=time.perf_counter()
    sol=solve_ivp(rhs,(0,DURATION),y0,method='Radau',jac=jac,rtol=2e-8*scale,
                  atol=atol*scale,t_eval=times,max_step=.008)
    if not sol.success: raise RuntimeError(sol.message)
    y=sol.y.T;a=y[:,:N];v=y[:,N:2*N];q=a@PHI.T;qd=v@PHI.T;theta=q[:,-1]
    energy=.5*np.einsum('ij,jk,ik->i',qd,M,qd)+.5*np.einsum('ij,jk,ik->i',q,K,q)+E0*2*np.sin(theta/2)**2
    wc=y[:,2*N]*E0;work=y[:,2*N+1]*E0;loss=y[:,2*N+2]*E0
    records=[port_control(t,HMAP@state,cap)for t,state in zip(times,y)]
    commands=np.array([r[1]for r in records]);actual=np.array([r[2]for r in records]);factors=np.array([r[3]for r in records])
    extent=np.max(abs(q@ROWS.T),axis=1)
    enclosure=wc+times*cert['additive_power_reserve_W']
    ledger=energy-energy[0]-work+loss
    result={'refined':refined,'duration_s':DURATION,'nfev':sol.nfev,'elapsed_s':time.perf_counter()-start,
     'actual_initial_energy_J':float(energy[0]),'actual_initial_tip_error_m':float(dq[-2]),
     'actual_initial_tip_velocity_error_m_s':float(dv[-2]),'initial_error_Q_J':upper(qe_ball),
     'maximum_energy_J':float(max(energy)),'maximum_work_observer_J':float(max(wc)),
     'maximum_energy_minus_certified_observer_J':float(max(energy-enclosure)),
     'maximum_work_cap_excess_J':float(max(wc-cert['command_work_cap_J'])),
     'maximum_centerline_Bernstein_bound_at_saved_times_m':float(max(extent)),
     'minimum_clearance_from_saved_Bernstein_bounds_m':float(CLEARANCE-max(extent)),
     'maximum_tip_displacement_m':float(max(abs(q[:,-2]))),'governor_active_fraction':float(np.mean(factors<1-1e-10)),
     'minimum_request_scale':float(min(factors)),'maximum_applied_field_T':float(max(abs(actual[:,1]))),
     'maximum_applied_tip_force_N':float(max(abs(actual[:,0]))),'maximum_energy_ledger_residual_J':float(max(abs(ledger))),
     'maximum_sensor_angle_error_rad':DELTA_THETA,'maximum_sensor_angular_rate_error_rad_s':DELTA_OMEGA,
     'maximum_sensor_tip_velocity_error_m_s':DELTA_V,
     'claim':'Finite-time robustness and whole-beam clearance under the proved work envelope; no noisy-trajectory KAM invariance claimed.'}
    assert energy[0] <= cert['initial_energy_bound_J']*(1+1e-12)
    assert max(wc-cert['command_work_cap_J'])<1e-13
    assert max(energy-enclosure)<1e-13
    assert max(extent)<CLEARANCE
    label='refined' if refined else 'primary'
    np.savez_compressed(HERE/(label+'_trajectory.npz'),time_s=times,modal_position=a,modal_velocity=v,
       tip_m=q[:,-2],tip_angle_rad=theta,energy_J=energy,work_observer_J=wc,
       energy_enclosure_J=enclosure,applied_force_N=actual[:,0],applied_field_T=actual[:,1],
       requested_scale=factors,command=commands,bernstein_centerline_bound_m=extent,
       input_work_J=work,dissipation_J=loss,power_balance_J=ledger,initial_nodal_position_error=dq,
       initial_nodal_velocity_error=dv)
    save_json(HERE/(label+'_results.json'),result)
    print(json.dumps(result,indent=2),flush=True)
    return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['certificate','primary','refined','compare','all'],default='all')
    args=parser.parse_args()
    cert=certify()
    if args.stage in ['primary','all']:run(cert)
    if args.stage in ['refined','all']:run(cert,True)
    if args.stage in ['compare','all']:
        a=np.load(HERE/'primary_trajectory.npz');b=np.load(HERE/'refined_trajectory.npz')
        comparison={'maximum_tip_difference_m':float(max(abs(a['tip_m']-b['tip_m']))),
         'maximum_energy_difference_J':float(max(abs(a['energy_J']-b['energy_J']))),
         'maximum_work_observer_difference_J':float(max(abs(a['work_observer_J']-b['work_observer_J']))),
         'relative_and_absolute_tolerances_refined_by':10}
        save_json(HERE/'refinement_comparison.json',comparison);print(comparison)
    deps={k:importlib.metadata.version(k)for k in ['numpy','scipy','numba','threadpoolctl','python-flint']}
    save_json(HERE/'environment.json',{'python':sys.version,'dependencies':deps,'inputs':{p.name:sha(p)for p in (HERE/'inputs').iterdir()if p.is_file()}})

if __name__=='__main__': main()
