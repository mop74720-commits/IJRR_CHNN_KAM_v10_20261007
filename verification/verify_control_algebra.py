"""Independent exact checks of the condensed control and clearance identities.

This is an algebraic consistency check, not a new KAM or continuum existence proof.
"""
from pathlib import Path
import json
import sympy as s

th,mu,F,U,E0,By=s.symbols('theta mu F U E0 By', real=True, nonzero=True)
G=F/mu*s.Matrix([[-s.sin(th),s.cos(th)],[s.cos(th),s.sin(th)]])
m=mu*s.Matrix([s.cos(th),s.sin(th)])
assert s.simplify(G*m-s.Matrix([0,F]))==s.zeros(2,1)
assert s.simplify(G.T-G)==s.zeros(2)
assert s.simplify(s.trace(G))==0
assert s.simplify(G.T*G-(F/mu)**2*s.eye(2))==s.zeros(2)
cmd=(U-E0*(th-s.sin(th)))/(mu*s.cos(th))
assert s.simplify(E0*(th-s.sin(th))+mu*cmd*s.cos(th)-U)==0
x1,x2,y1,y2=s.symbols('x1 x2 y1 y2')
xs=s.Matrix([x1,x2]);ys=s.Matrix([y1,y2])
# Nonseparable polynomial includes mixed momenta; identity must retain h_YY.
h=(x1*x1+2*x2*x2+y1*y1+2*y2*y2)/2+x1*x1*y2*y2/7+x2*y1**3/11+y1*y1*y2*y2/13
hx=s.Matrix([s.diff(h,x) for x in xs]);hy=s.Matrix([s.diff(h,y) for y in ys])
flow=xs.col_join(ys);field=hy.col_join(-hx)
direct=hy.jacobian(flow)*field
claimed=hy.jacobian(xs)*hy-hy.jacobian(ys)*hx
assert s.simplify(direct-claimed)==s.zeros(2,1)
ell,xi=s.symbols('ell xi',positive=True)
assert s.integrate((ell-xi)**2,(xi,0,ell))==ell**3/3
EI=s.Rational(5_000_000)*s.pi*s.Rational(6,10000)**4/4
L=s.Rational(8,100);c=s.Rational(4,1000)
Ec=3*EI*c*c/(2*L**3)
bound=s.sqrt(2*L**3*(s.Rational(95,100)*Ec)/(3*EI))
report={'status':'PASS','checks':['magnetic_force','symmetric_traceless_gradient',
  'gradient_spectral_norm','magnetic_torque_command','nonseparable_Hamiltonian_acceleration',
  'clamped_beam_Cauchy_Schwarz_kernel'],
  'clearance_energy_nJ':float(Ec*10**9),'energy_cap_centerline_mm':float(bound*1000),
  'scope':'Exact algebra and parameter substitution for main application equations; model assumptions and proof hypotheses remain necessary.'}
out=Path(__file__).with_name('core_control_algebra_validation.json')
out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
