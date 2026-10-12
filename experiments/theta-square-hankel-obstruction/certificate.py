#!/usr/bin/env python3
"""Bound the actual square-variable theta derivative determinant.

Computer-assisted sign certificate, not a Lean proof or RH counterexample.
The omitted-series Cauchy bounds are given in README.md.
"""
import argparse
import json
import math
from pathlib import Path
from flint import arb, acb_series, arb_mat, ctx

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--bits',type=int,default=256)
args=parser.parse_args()
if args.bits<200: parser.error('at least 200 bits required')
ctx.prec=args.bits
pi=arb.pi()
terms=16
error=arb(2)**(-180)

def K_series(u,cap):
    e=(4*u).exp()
    return sum((pi*m*m*(2*pi*m*m*e-3)*(5*u-pi*m*m*e).exp()
                for m in range(1,terms+1)),acb_series([0],cap))

def enclose_real(value):
    if not value.is_finite() or not value.imag.contains(0):
        raise ArithmeticError('nonfinite or nonreal series coefficient')
    return arb(value.real.mid(),value.real.rad()+error)

ctx.cap=17
u=acb_series([0,1],17)
ku=K_series(u,17)
endpoint=[enclose_real(2*math.factorial(j)*ku[2*j]) for j in range(9)]
ctx.cap=9
t0=arb(1)/10000
t=acb_series([t0,1],9)
u=t.sqrt()
# Evenize the finite sum: the full original kernel is even, not its truncation.
s=K_series(u,9)+K_series(-u,9)
positive_point=[enclose_real(math.factorial(j)*s[j]) for j in range(9)]
rows=[]
for point,ds in [('0',endpoint),('1/10000',positive_point)]:
    if not ds[0]>0: raise ArithmeticError('positivity of S is not enclosed')
    matrix=arb_mat([[ds[i+j]/ds[0] for j in range(5)] for i in range(5)])
    determinant=matrix.det()
    tail_bound=arb(2)**(-279)*math.factorial(8)/(arb(1)/1024-t0)**8
    row={'T':point,'derivatives':[str(z) for z in ds],
         'normalized_delta5':str(determinant),
         'checks':{'S_positive':bool(ds[0]>0),
                   'delta5_over_S5_less_than_minus_10_to_25':bool(determinant < -arb(10)**25),
                   'cauchy_derivative_bound_less_than_added_radius':bool(tail_bound<error)}}
    rows.append(row)
    print(json.dumps(row),flush=True)
result={'kind':'computer-assisted interval certificate, not Lean kernel proof',
        'normalization':{'S':'2K(sqrt(T))','T':'x^2/4',
                         'delta5':'det[S^(i+j)(T)]_(i,j=0..4)'},
        'python_flint_version':__import__('flint').__version__,'bits':args.bits,
        'theta_terms':terms,'derivative_orders':[0,8],
        'derivative_error_radius':'2^-180',
        'analytic_tail_bound':{'u_circle_radius':'1/32','K_tail_bound':'2^-280',
                              'T_circle_radius':'1/1024','evenized_tail_bound':'2^-279'},
        'points':rows}
if not all(check for row in rows for check in row['checks'].values()):
    raise SystemExit('certificate failed')
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,indent=2)+'\n')
