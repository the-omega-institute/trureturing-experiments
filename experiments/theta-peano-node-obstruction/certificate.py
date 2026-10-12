#!/usr/bin/env python3
"""Enclose a degree-seven theta node and the local Peano obstruction.

The two analytic truncation bounds are derived in README.md. This is an
Arb interval certificate, not a Lean proof or an RH/Robin counterexample.
Requires python-flint==0.8.0; output is an explicit command-line argument.
"""
import argparse
import json
import math
from pathlib import Path
from flint import arb, acb, acb_series, ctx
parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--bits',type=int,default=200)
args=parser.parse_args()
if args.bits<180: parser.error('at least 180 bits required')
ctx.prec=args.bits
ctx.cap=13
pi=arb.pi()
terms=16
I=[]
def K(z):
 e=(4*z).exp()
 return sum((pi*m*m*(2*pi*m*m*e-3)*(5*z-pi*m*m*e).exp() for m in range(1,terms+1)),acb(0))
for k in range(8):
 value=acb(0)
 for lo,hi in [('0','0.25'),('0.25','0.5'),('0.5','1'),('1','2')]:
  value+=acb.integral(lambda z,analytic: 4*z**(2*k)*K(z),acb(lo),acb(hi),rel_tol=arb(2)**(-args.bits+30),abs_tol=arb(2)**(-args.bits+30),eval_limit=100000)
 if not value.is_finite() or not value.imag.contains(0): raise ArithmeticError('invalid integral enclosure')
 moment=arb(value.real.mid(),value.real.rad()+arb(2)**(-500))
 I.append(moment)
 print('I',k,str(moment),flush=True)
u=acb_series([0,1],13)
e=(4*u).exp()
ks=sum((pi*m*m*(2*pi*m*m*e-3)*(5*u-pi*m*m*e).exp() for m in range(1,terms+1)),acb_series([0],13))
Sj=[]
for j in range(7):
 val=2*math.factorial(j)*ks[2*j]
 if not val.is_finite() or not val.imag.contains(0): raise ArithmeticError('invalid series enclosure')
 Sj.append(arb(val.real.mid(),val.real.rad()+arb(2)**(-200)))
 print('S derivative',j,str(Sj[-1]),flush=True)
def ff(x,k):
 value=arb(1)
 for i in range(k): value*=x-i
 return value
def binomial(x,k): return ff(x,k)/math.factorial(k)
n=7
coeffs=[(-1)**(n+j)*math.factorial(n)*binomial(arb('6.5'),n-j)/math.factorial(j)*I[j] for j in range(n+1)]
def polynomial(v):
 q=arb(0)
 for c in coeffs: q=q*v+c
 return q
vlo=arb(3115806)/100000000
vhi=arb(3115808)/100000000
qlo=polynomial(vlo); qhi=polynomial(vhi)
v=arb((vlo.mid()+vhi.mid())/2,(vhi.mid()-vlo.mid())/2+vlo.rad()+vhi.rad())
D=sum((math.comb(6,j)*v**(-(6-j))*(-1)**j*ff(arb('5.5'),j)*I[6-j] for j in range(7)),arb(0))
Z=v**6*v.sqrt()*arb('6.5').gamma()
H0=sum((math.comb(6,j)*v**(-(6-j))*Sj[j] for j in range(7)),arb(0))
endpoint=H0-D/Z
checks={'q_at_lower_negative':bool(qlo<0), 'q_at_upper_positive':bool(qhi>0), 'D_less_than_minus_one_eighth':bool(D < -arb(1)/8), 'endpoint_less_than_minus_sixty_million':bool(endpoint < -60000000), 'product_greater_than_nine_million':bool(D*endpoint > 9000000)}
result={'kind':'computer-assisted interval certificate, not Lean kernel proof','normalization':{'S':'2K(sqrt(T))','T':'x^2/4','node':'v=7*lambda, sourceQ_7(lambda)=0','weights':'unnormalized W=T^(11/2)*exp(-T/v), H=(exp(T/v)*S)^(6)'},'python_flint_version': __import__('flint').__version__, 'bits':args.bits,'theta_terms':terms,'source_tail_bounds':{'moments_abs':'2^-500 (orders 0..7)','S_derivatives_abs':'2^-200 (orders 0..6)'},'node_bracket':['3115806/100000000','3115808/100000000'],'moment_integral_cutoff':'u=2', 'I':[str(z) for z in I],'S_derivatives_0':[str(z) for z in Sj],'qlo':str(qlo),'qhi':str(qhi),'D':str(D),'Z':str(Z),'H0':str(H0),'H0_minus_D_over_Z':str(endpoint),'D_times_endpoint':str(D*endpoint),'checks':checks}
print(json.dumps(result,indent=2),flush=True)
if not all(checks.values()): raise SystemExit('certificate failed')
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(result,indent=2)+'\n')
