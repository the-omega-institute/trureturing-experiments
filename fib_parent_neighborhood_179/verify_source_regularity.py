#!/usr/bin/env python3
"""Exact source-family checks for the FIB parent-neighborhood increment."""
from __future__ import annotations
import hashlib, json, sys
from fractions import Fraction
from pathlib import Path
assert sys.version_info >= (3, 12), sys.version
CERT = Path(__file__).resolve().with_name('source_regularity_macaulay.certificate.json')
OUT = Path(__file__).resolve().with_name('source_checks.json')
Z = Fraction(0)
class Poly:
    def __init__(self,d=None): self.d={tuple(k):Fraction(v) for k,v in (d or {}).items() if Fraction(v)}
    @staticmethod
    def const(c): return Poly({(0,0,0):Fraction(c)})
    @staticmethod
    def var(i):
        e=[0,0,0]; e[i]=1; return Poly({tuple(e):Fraction(1)})
    def __add__(self,o):
        o=aspoly(o); d=dict(self.d)
        for e,c in o.d.items(): d[e]=d.get(e,Z)+c
        return Poly(d)
    __radd__=__add__
    def __neg__(self): return Poly({e:-c for e,c in self.d.items()})
    def __sub__(self,o): return self+(-aspoly(o))
    def __rsub__(self,o): return aspoly(o)-self
    def __mul__(self,o):
        o=aspoly(o); d={}
        for e,c in self.d.items():
            for f,g in o.d.items():
                h=tuple(e[i]+f[i] for i in range(3)); d[h]=d.get(h,Z)+c*g
        return Poly(d)
    __rmul__=__mul__
    def __pow__(self,k):
        out=Poly.const(1); x=self
        while k:
            if k&1: out=out*x
            x=x*x; k//=2
        return out
    def deriv(self,i):
        d={}
        for e,c in self.d.items():
            if e[i]:
                h=list(e); h[i]-=1; h=tuple(h); d[h]=d.get(h,Z)+c*e[i]
        return Poly(d)
    def eq(self,o): return (self-aspoly(o)).d=={}
def aspoly(x): return x if isinstance(x,Poly) else Poly.const(x)
def from_terms(terms):
    d={}
    for exp,c in terms:
        e=tuple(exp); d[e]=d.get(e,Z)+Fraction(c)
    return Poly(d)
def check(name,ok):
    if not ok: raise AssertionError(name)
    CHECKS.append(name)
CHECKS=[]
data=json.loads(CERT.read_text())
check('certificate_sha256',hashlib.sha256(CERT.read_bytes()).hexdigest()=='f4a6f91da3e833ede44549a532a9b32d7695e3e24d744b8e6bac8e5e46f654ab')
check('variables',tuple(data['variables'])==('a','n','u'))
gens=[from_terms(g) for g in data['generators']]; mult=[from_terms(m) for m in data['multipliers']]
a,n,u=(Poly.var(i) for i in range(3))
h=u*(u+1); ell=u**2+3*u-4
p=n+a-3*u-1; q=n+2*a-4*u-3; rr=1+2*u-a
Q=p*(1-p)-q*rr
Phi=(3*h-10)*a+(8*u**2+18*u-30)*n-23*h
Delta=n*h+a*n*ell-4*a*h
Tnum=(87*n-203*a)*Phi+1225*Delta
JQ=[Q.deriv(i) for i in range(3)]; JT=[Tnum.deriv(i) for i in range(3)]
mins=[JQ[0]*JT[1]-JQ[1]*JT[0],JQ[0]*JT[2]-JQ[2]*JT[0],JQ[1]*JT[2]-JQ[2]*JT[1]]
for i,(got,want) in enumerate(zip(gens,[Q,Tnum,*mins])): check(f'generator_{i}_original_formula',got.eq(want))
identity=Poly.const(0)
for m,g in zip(mult,gens): identity += m*g
check('macaulay_sum_equals_one',identity.eq(1))
check('max_product_degree_8',max(sum(e) for m,g in zip(mult,gens) for e in (m*g).d)==8)
check('five_generators',len(gens)==len(mult)==5)
th=Fraction(3,7); xi=Fraction(-35,29); T=(th*n-a)*Phi-5*xi*Delta
check('anchor_scale', (Tnum-203*T).eq(0))
def leading_on_tail(polynomial, ca, cn):
    out = Poly.const(0)
    for e, coefficient in polynomial.d.items():
        out += coefficient * (ca*u)**e[0] * (cn*u)**e[1] * u**e[2]
    return out.d.get((0, 0, 4), Z)

for name, ca, cn, expected in (
    ('q_tail', 2, 1, (14, -28, -10)),
    ('r_tail', 1, 2, (38, -19, -10)),
):
    got = tuple(leading_on_tail(f, ca, cn) for f in (n*Phi, -a*Phi, -5*Delta))
    check(name+'_general_leading_coefficients', got == expected)
    anchor = got[0]*th + got[1] + got[2]*xi
    perturbation = (abs(got[0]) + abs(got[2])) / Fraction(100)
    check(name+'_uniform_nonzero_on_closed_box', abs(anchor) > perturbation)

corner = sum(c for e,c in T.d.items() if e[2] == 0)
check('corner_u0', corner == Fraction(-260,203))
check('corner_uniform_nonzero_on_closed_box', abs(corner) > Fraction(60,100))
check('xi_uniform_nonzero_on_closed_box', abs(xi) > Fraction(1,100))
result = {
    'check_count': len(CHECKS),
    'checks': CHECKS,
    'certificate_sha256': hashlib.sha256(CERT.read_bytes()).hexdigest(),
    'anchor': {'theta':'3/7', 'xi':'-35/29', 'scale':'Tnum=203*T'},
    'macaulay_identity': 'sum_i multiplier_i*generator_i=1',
    'max_product_degree': 8,
    'tail_leading_coefficients': ['14*theta-28-10*xi', '38*theta-19-10*xi'],
    'u0_corner': '-260/203',
    'uniform_box_radius': '1/100',
    'scope': 'exact arithmetic; source regularity and compact-envelope coefficients; not response immersion, neighborhood proof, or Lean verification'
}
OUT.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, sort_keys=True))
