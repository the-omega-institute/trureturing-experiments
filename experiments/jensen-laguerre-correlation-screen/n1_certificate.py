"""Directed screen of one proposed sufficient condition, not an RH test."""
import ast
import json
import hashlib
import flint
import sys
from pathlib import Path
from flint import acb, arb, ctx, fmpq

ctx.prec = 160
BASE = Path(__file__).resolve().parents[2] / 'docs/reports/theta-mixed-matrix'
supplier_bytes = (BASE / 'derivative_bandwidth.py').read_bytes()
supplier_sha256 = hashlib.sha256(supplier_bytes).hexdigest()
assert supplier_sha256 == 'a74b0d95fdc546c240b605df15ed69fe7c50f4963deeff2a4ad6498fb687bb40'
ns = {'fmpq': fmpq}
nodes = [node for node in ast.parse(supplier_bytes).body
         if isinstance(node, ast.FunctionDef) and node.name == 'polys']
assert len(nodes) == 1
exec(compile(ast.Module(body=nodes, type_ignores=[]), 'reused_theta_derivative_polys', 'exec'), ns)
polys = ns['polys']
pi = arb.pi()
terms = 8
coeff = {(j, a): polys(a, j) for j in range(3)
         for a in (fmpq(9, 2), fmpq(5, 2))}

def tail_series(power, start=terms + 1):
    ratio = arb(fmpq(start + 1, start))**power * (-pi * (2 * start + 1)).exp()
    assert ratio < 1
    return arb(start)**power * (-pi * (start * start - 1)).exp() / (1 - ratio)

def constant(j, omitted):
    total = arb(0)
    for a, pref, degree in ((fmpq(9, 2), 4*pi*pi, 4),
                           (fmpq(5, 2), 6*pi, 2)):
        for k, v in enumerate(coeff[j, a]):
            power = degree + 2*k
            series = tail_series(power)
            if not omitted:
                series += sum(arb(m)**power * (-pi * (m*m - 1)).exp()
                              for m in range(1, terms + 1))
            total += pref * abs(arb(v)) * pi**k * series
    return total.upper()

ct = [constant(j, True) for j in range(3)]
cs = [constant(j, False) for j in range(3)]
assert cs[0] < 64 and cs[1] < 800

def jet(x, order):
    u = (2*x).exp()
    value = arb(0)
    for m in range(1, terms + 1):
        z = pi*m*m*u
        ps = []
        for a in (fmpq(9, 2), fmpq(5, 2)):
            q = arb(0)
            for c in reversed(coeff[order, a]):
                q = q*z + arb(c)
            ps.append(q)
        value += (4*pi*pi*m**4 * (arb(fmpq(9,2))*x).exp()*ps[0]
                  - 6*pi*m*m * (arb(fmpq(5,2))*x).exp()*ps[1]) * (-z).exp()
    err = ct[order] * (arb(fmpq(9,2)+2*order)*x).exp() * (-pi*u).exp()
    return value + arb(0, err.upper())

tail = arb(2)**-1000
# For x>=3 and 20<alpha<22, all six raw moments have this common tail.
# R_3/2(y)<=exp(y) when y>=1, and y<=11*x^2.
# 22*x^2<=e^(2x)/2 and (25/2)*x<=e^(2x)/8 for x>=3;
# both inequalities propagate because e^(2x)/x^2 and e^(2x)/x increase there.
# x^6<=exp(6x), so every integrand is <=121*864*exp(-2*e^(2x)).
# Integrate using e^(2x)>=e^6*(1+2*(x-3)).
assert arb(198) < arb(6).exp()/2
assert arb(fmpq(75,2)) < arb(6).exp()/8
assert 121*864*(-2*arb(6).exp()).exp()/(4*arb(6).exp()) < tail

def truncated_phi(z):
    u = (2*z).exp()
    return sum((4*pi*pi*m**4*(arb(fmpq(9,2))*z).exp()
                -6*pi*m*m*(arb(fmpq(5,2))*z).exp()) * (-pi*m*m*u).exp()
               for m in range(1, terms+1))

mom = []
for j in (0,1,2):
    core = acb(0)
    for k in range(12):
        def fun(z, analytic, j=j):
            return truncated_phi(z)*z**(2*j)
        core += acb.integral(fun, acb(fmpq(k,4)), acb(fmpq(k+1,4)),
                             abs_tol=arb(2)**-120, rel_tol=arb(2)**-120)
    # Complete m>=9 source tail on [0,3], and all spatial x>3.
    err = 3*3**(2*j)*ct[0]*arb(fmpq(27,2)).exp() + tail
    mom.append(core.real + arb(err/2, (err/2).upper()))
assert all(v > 0 and v.is_finite() for v in mom)
alpha = mom[0]/mom[1]
assert alpha > 20 and alpha < 22
node = 2*alpha
print(json.dumps({'phase':'source_moments','alpha':str(alpha),'node':str(node)}), flush=True)

def hull(lo, hi):
    return arb((lo+hi)/2, ((hi-lo)/2).upper())

def Rpoint(y):
    r = y.sqrt()
    return pi.sqrt()/2*r.erf()*y.exp()/r-1

def Rrange(y):
    # R(y)=y int_0^1 v^(1/2)e^(y(1-v))dv is increasing.
    return hull(Rpoint(y.lower()).lower(), Rpoint(y.upper()).upper())

grid = 131072
eps = fmpq(1,4096)
epsa = arb(eps)
length = fmpq(3)-eps
width = arb(length/grid)
tot = [arb(0) for _ in range(6)]
for k in range(grid):
    lo = eps + length*fmpq(k,grid)
    hi = eps + length*fmpq(k+1,grid)
    x = arb((lo+hi)/2,arb((hi-lo)/2).upper())
    H = jet(x,0)+jet(x,1)/(alpha*x)
    density = abs(H)*x**2
    Y = alpha*x**2/2
    R = Rrange(Y)
    for j, weight in enumerate((1,Y,Y**2,R,R**2,Y*R)):
        tot[j] += width*density*weight
    if (k+1)%32768 == 0:
        print(json.dumps({'phase':'directed_rectangles','completed':k+1,'grid':grid}),flush=True)

# Complete the small interval [0,eps] by Phi'(0)=0 and the mean value bound.
box = arb(epsa/2,(epsa/2).upper())
hbound = abs(jet(box,0)).upper()+abs(jet(box,2)).upper()/alpha.lower()
mass0 = hbound*epsa**3/3
ymax = alpha.upper()*epsa**2/2
rmax = ymax*ymax.exp()/arb(fmpq(3,2))
for j, weight in enumerate((1,ymax,ymax*ymax,rmax,rmax*rmax,ymax*rmax)):
    err = (mass0*weight+tail).upper()
    tot[j] += arb(err/2,(err/2).upper())

M = tot[0]
EY,EY2,ER,ER2,EYR = [v/M for v in tot[1:]]
varY,varR,cov = EY2-EY**2,ER2-ER**2,EYR-EY*ER
assert M > 0 and varY > 0 and varR > 0 and cov > 0, str((M,varY,varR,cov,tot,mass0))
# At alpha=M0/M2 the signed raw mass is exactly zero; integration by parts:
# int Y*H*x^2 dx = alpha*M4/2-3*M2/2.
ESY = (alpha*mom[2]/2-arb(fmpq(3,2))*mom[1])/M
c2 = ESY**2/varY
r2 = cov**2/(varY*varR)
score = c2+r2-1
result = {'scope':'One exact n=1 actual theta Jensen node; stronger correlation condition only',
          'precision_bits':ctx.prec,'theta_terms_retained':terms,'grid':grid,
          'runtime':{'python':sys.version.split()[0],'python_flint':flint.__version__},
          'supplier':{'path':'docs/reports/theta-mixed-matrix/derivative_bandwidth.py',
                      'sha256':supplier_sha256},
          'source_moments_half':[str(v) for v in mom],
          'alpha':str(alpha),'node':str(node),'source_tail_constants':[str(v) for v in ct],
          'source_envelopes':[str(v) for v in cs],
          'spatial_raw_moment_tail':str(tail),'small_core_mass_upper':str(mass0.upper()),
          'raw_moments':[str(v) for v in tot],
          'c_squared':str(c2),'r_squared':str(r2),'criterion_score':str(score),
          'criterion_negative_certified':bool(score < 0),
          'score_endpoints_dyadic':{
              'lower':[str(v) for v in score.lower().man_exp()],
              'upper':[str(v) for v in score.upper().man_exp()]},
          'no_RH_or_response_sign_conclusion':True,'Lean':'not formalized'}
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
assert score < 0, 'Directed enclosure did not settle the proposed condition'
