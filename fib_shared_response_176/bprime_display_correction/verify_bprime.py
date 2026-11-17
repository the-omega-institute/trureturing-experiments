"""Uniform Bprime theorem 9 certificate and literal-source regression.
Reuses the prior exact outward integer arithmetic; no floating point, sampling,
Lean verification, independent review, or all-response recovery is asserted.
"""
import hashlib, json, re, sys
from pathlib import Path
from fractions import Fraction as F
from math import comb
import sympy as sp
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'dependencies' / 'fib_x_section8_display_correction'))
from verify_section8 import Box, centered_poly as centered_multi, polybox, display, intervals, rational
assert sys.version_info >= (3, 12)
def centered_poly(poly, variables, ranges):
    if len(variables) != 1:
        return centered_multi(poly, variables, ranges)
    var=variables[0]
    center=sum(ranges[var].endpoints())/2
    shifted=sp.expand(poly.subs(var,var+sp.Rational(center)))
    return polybox(shifted,variables,{var:ranges[var]-center})
THETA, XI = F(3,7), F(-35,29)
path = Path(sys.argv[1])
data = path.read_bytes()
section = data.decode().split('## 9. ', 1)[1].split('## 10. ', 1)[0]
checks=[]
measurements={}
def check(name, ok, **details):
    checks.append(dict(name=name, passed=bool(ok), **details))
    if not ok: print('FAIL', name, flush=True)
u, b, n, a = sp.symbols('u b n a')
h, ell = u*(u + 1), u*u + 3*u - 4
al, be = 3*h - 10, 8*u*u + 18*u - 30
phi = al*a + be*n - 23*h
delta = n*h + a*n*ell - 4*a*h
A = (10 - 3*b)*u*u + (8 - 3*b)*u + 10*b - 32
B = (10 - 8*b)*u*u + (19 - 18*b)*u + 30*b - 36
C = (23*b - 29)*h
Q = sp.Matrix([[1, -sp.Rational(1, 2), -u - 1],
               [-sp.Rational(1, 2), -1, 2*u + 1], [-u - 1, 2*u + 1, -u*u + u + 1]])
T = sp.Matrix([[-al, (THETA*al - be - 5*XI*ell)/2, (23 + 20*XI)*h/2],
               [(THETA*al - be - 5*XI*ell)/2, THETA*be, -(23*THETA + 5*XI)*h/2],
               [(23 + 20*XI)*h/2, -(23*THETA + 5*XI)*h/2, 0]])
vec = sp.Matrix([A, B, C])
rho = sp.expand((vec.T*Q.adjugate()*vec)[0])
sigma = sp.expand((vec.T*T.adjugate()*vec)[0])
tau = sp.expand((vec.T*((Q + T).adjugate() - Q.adjugate() - T.adjugate())*vec)[0])
R = sp.Poly(tau*tau - 4*rho*sigma, u, b)
p = sp.Matrix([[n + a - 3*u - 1, n + 2*a - 4*u - 3],
               [1 + 2*u - a, 2 + 3*u - n - a]])
x = sp.Matrix([a, n, 1])
check('VI 24.3 source determinant', sp.expand((x.T*Q*x)[0] - p.det()) == 0)
check('VI 24.4 parent polynomial', sp.expand((x.T*T*x)[0] - ((THETA*n - a)*phi - 5*XI*delta)) == 0)
psi6 = (n + 2*a)*(10*u*u + 19*u - 36) - 10*a*ell - 29*h
check('VI 24.2 actual M response', sp.expand(A*a + B*n + C - (psi6 - b*phi)) == 0)
# Independent literal action matrices and four-node interpolation establish
# the formulas for the same parent history and the actual sixth action.
M, J = sp.Matrix([[0, 1], [1, 1]]), sp.Matrix([[0, 1], [1, 0]])
acts = [sp.eye(2)]
for action in (M, M, J, M, M):
    acts.append(action*acts[-1])
trace = lambda action: sp.expand(sp.trace(action*p))
check('8.1 source traces on literal MMJMM parent',
      [trace(action) for action in acts] == [1, u, u + 1, n, a, n + a])
check('8.3/8.4 actual sixth and seventh trace factors',
      [trace(M*acts[-1]), trace(M*M*acts[-1]), trace(J*M*acts[-1]),
       trace(J*acts[-1]), trace(J*J*acts[-1]), trace(M*J*acts[-1])] ==
      [n + 2*a, 2*n + 3*a, 4*n - u + 1, 1 + 2*n, n + a, u + 2*a])
drift = sp.Symbol('drift')
nodes = [-drift, -drift/u, -drift/(u + 1), (6 - drift)/n]

def interpolate(clock):
    return sum(value*sp.prod(sp.Rational(clock - j, i - j)
                            for j in range(4) if j != i)
               for i, value in enumerate(nodes))

V6, X7 = 10*u*u + 19*u - 36, 20*u*u + 34*u - 70
g6, g7 = 120/n + drift*(V6/h - 20/n), 210/n + drift*(X7/h - 35/n)
check('8.2 gain from whole cubic interpolation',
      sp.cancel(interpolate(6) - g6) == 0 and sp.cancel(interpolate(7) - g7) == 0)
theta_local, xi_local = a/n + drift*delta/(24*h*n), drift*phi/(120*h)
q_src = 5*delta/(n*phi)
psi7 = 4*(2*n + 3*a)*X7 - 105*a*ell - 381*h
gamma = 120*h/phi*(1 + (1 + 2*n)*(V6/h - 20/n))
zeta = (120*h/phi*(1 + (n + a)*(X7/h - 35/n)) - 210*q_src)/30
identities = [
    (n + 2*a)*g6 + drift - (120*(1 + 2*theta_local) + 120*xi_local*psi6/phi),
    (2*n + 3*a)*g7 + drift - (210*(2 + 3*theta_local) + 30*xi_local*psi7/phi),
    (1 + 2*n)*g6 + drift - (240 + 120/n + xi_local*gamma),
    (n + a)*g7 + drift - (210*(1 + theta_local) + 30*xi_local*zeta),
    a*interpolate(4) + drift - 24*theta_local,
    (n + a)*interpolate(5) + drift - 60*(1 + theta_local + xi_local),
]
check('8.1--8.7 response and parent identities on full source coordinates',
      all(sp.cancel(identity) == 0 for identity in identities))
xa = sp.Matrix([-B*n - C, A*n, A])
qa = sp.Poly((xa.T*Q*xa)[0], n)
ta = sp.Poly((xa.T*T*xa)[0], n)
q2, q1, q0 = qa.all_coeffs()
t2, t1, t0 = ta.all_coeffs()
N, D = sp.expand(t2*q0 - q2*t0), sp.expand(q2*t1 - t2*q1)
nhat, dhat = [sp.cancel(-203*z/A) for z in (N, D)]
check('8.17 polynomial factors and degrees',
      sp.denom(nhat) == 1 and sp.denom(dhat) == 1 and
      (sp.degree(nhat, u), sp.degree(dhat, u), sp.degree(nhat, b), sp.degree(dhat, b)) == (10, 9, 3, 3))
chart_resultant = sp.expand((q2*t0 - q0*t2)**2 - (q2*t1 - q1*t2)*(q1*t0 - q0*t1))
check('8.17 chart resultant equals A^4 R', sp.Poly(chart_resultant - A**4*R.as_expr(), u, b).is_zero)
check('8.17 reconstructed n solves both quadratics at every R root', all(
    sp.div(sp.Poly(lead*nhat**2 + middle*nhat*dhat + const*dhat**2, u, b), R)[1].is_zero
    for lead, middle, const in ((q2, q1, q0), (t2, t1, t0))))
check('8.16 degree 16', R.degree(u) == 16)

b_ends=intervals(display(section,'9.1'))[0]
b_box=Box(*b_ends)
check('literal Bprime endpoints', b_ends==(F(16423371109,10**10),F(1642337111,10**9)))
base=data.decode().split('## 8. ',1)[1]
broad=re.findall(r'\[([^,\]]+),([^\])]+)[\])]',display(base,'8.16'))
boundaries=[rational(x[0]) for x in broad]
check('complete positive half-line partition',len(broad)==13 and boundaries[0]==0 and all(rational(broad[i][1])==boundaries[i+1] for i in range(12)) and broad[-1][1]==r'\infty')
expected=['+'*17,'+'*3+'-'*14,'-'*17,'-'*3+'+'*14,'+'*17,'+'*10+'-'*7,'-'*17,'-'*7+'+'*10,'+'*17,'+'*7+'-'*10,'-'*17,'-'*7+'+'*10,'+'*17]
shown=[''.join(s*int(k or 1) for s,k in re.findall(r'([+-])(?:\^\{?(\d+)\}?)?',word)) for word in display(section,'9.4').replace(r'\\','').split(',')]
check('literal 9.4 coefficient sign strings',shown==expected)
coeff=[[F(R.coeff_monomial(u**k*b**j)) for k in range(17)] for j in range(R.degree(b)+1)]
for row,lo in enumerate(boundaries):
    signs=[]
    for power in range(17):
        by=[]
        for j,cs in enumerate(coeff):
            val=F(0)
            for k,c in enumerate(cs):
                if row<12:
                    hi=boundaries[row+1]
                    val+=c*sum(F(comb(k,e)*comb(16-k,power-e))*lo**(k-e)*hi**e for e in range(max(0,power-16+k),min(k,power)+1))
                elif k>=16-power: val+=c*comb(k,16-power)*lo**(k-16+power)
            by.append(val)
        val=Box(0)
        for cc in reversed(by): val=val*b_box+cc
        signs.append(val.sign())
    check(f'uniform Descartes row {row+1}', ''.join(signs)==expected[row])
    endsign=[centered_poly(R.as_expr().subs(u,sp.Rational(lo)),[b],{b:b_box}).sign()]
    endsign.append(centered_poly(R.as_expr().subs(u,sp.Rational(boundaries[row+1])) if row<12 else sp.Poly(R.as_expr(),u).LC(),[b],{b:b_box}).sign())
    check(f'uniform broad endpoints {row+1}',endsign==[expected[row][0],expected[row][-1]])
rootrows=re.findall(r'([1-6])&\[([0-9.]+),([0-9.]+)\]&([+-]{4})', display(section,'9.5'))
check('literal six source boxes',len(rootrows)==6)
terminal=intervals(display(section,'9.7'))
check('literal MJ intervals ordered',len(terminal)==2 and all(lo<hi for lo,hi in terminal))
roots=[]
for ix,(_,lo,hi,source_signs) in enumerate(rootrows):
    lo,hi=F(lo),F(hi); ub=Box(lo,hi); ranges={u:ub,b:b_box}
    endsign=[centered_poly(R.as_expr().subs(u,sp.Rational(ep)),[b],{b:b_box}).sign() for ep in (lo,hi)]
    row=2*ix+1
    check(f'root {ix+1} unique within broad cell',boundaries[row]<=lo<hi<=boundaries[row+1] and endsign==[expected[row][0],expected[row][-1]])
    av=centered_poly(A,[u,b],ranges); dv=centered_poly(dhat,[u,b],ranges)
    check(f'root {ix+1} chart and linear subresultant legal',av.sign()!='0' and dv.sign()!='0')
    nv=centered_poly(nhat,[u,b],ranges)/dv
    aa=-(centered_poly(B,[u,b],ranges)*nv+centered_poly(C,[u,b],ranges))/av
    leads=[centered_poly(z,[u,b],ranges).sign() for z in (q2,t2)]
    check(f'root {ix+1} quadratic leading coefficients legal',all(s!='0' for s in leads))
    pvals=[nv+aa-3*ub-1,nv+2*aa-4*ub-3,1+2*ub-aa,2+3*ub-nv-aa]
    check(f'root {ix+1} all source filters', ''.join(z.sign() for z in pvals)==source_signs)
    ph=(3*ub*(ub+1)-10)*aa+(8*ub**2+18*ub-30)*nv-23*ub*(ub+1)
    check(f'root {ix+1} Phi and n nonzero',ph.sign()!='0' and nv.sign()!='0')
    result=dict(root=ix+1,endpoint_signs=endsign,source_signs=source_signs,n=nv.outward(12),a=aa.outward(12),Phi=ph.outward(12))
    if ix in (1,3):
        check(f'root {ix+1} positive legal traces', nv.sign()==aa.sign()==ub.sign()=='+' and ph.sign()=='-')
        drift=120*ub*(ub+1)*XI/ph
        gain7=210/nv+drift*((20*ub**2+34*ub-70)/(ub*(ub+1))-35/nv)
        mj=drift+(4*nv-ub+1)*gain7
        shown=terminal[0 if ix==1 else 1]
        check(f'root {ix+1} literal displayed MJ enclosure',shown[0]<mj.endpoints()[0]<=mj.endpoints()[1]<shown[1],displayed=list(map(str,shown)),computed=mj.outward(12))
        result['actual_MJ']=mj.outward(12)
    roots.append(result)
check('literal separated MJ images',terminal[1][1]<terminal[0][0])
# The known defect is rejected even if the disjointness-only check passes.
bad=section.replace(r'2\cdot10^{11}',r'2\cdot10^{13}',1)
badints=intervals(display(bad,'9.7'))
check('regression rejects reversed original display',badints[0][0]>badints[0][1])
measurements['roots']=roots
measurements['displayed_MJ']=list(map(lambda x:list(map(str,x)),terminal))
report=dict(passed=all(c['passed'] for c in checks),checks=checks,measurements=measurements,document_sha256=hashlib.sha256(data).hexdigest(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),arithmetic='symbolic rational identities and outward integer intervals at scale 10^80; full Bprime, no samples',scope='fixed father (3/7,-35/29), all complete worlds in every M-child fibre b in Bprime',not_established=['independent review','Lean verification','other responses','all-parent or ancestor compatibility','root8/9 or nine-read upper bound'])
Path(sys.argv[2]).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(passed=report['passed'],checks=len(checks),failures=[c for c in checks if not c['passed']],measurements=measurements),indent=2),flush=True)
sys.exit(0 if report['passed'] else 1)
