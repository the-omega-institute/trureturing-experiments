#!/usr/bin/env python3
"""Check the retained discriminant factorization without recomputing a resultant."""
from __future__ import annotations
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

assert sys.version_info >= (3, 12), sys.version
ROOT = Path(__file__).resolve().parent
manifest = json.loads((ROOT/'coefficient_manifest.json').read_text())
assert manifest['input_commit'] == '7c45a21dbea5cdc0974f878a96dc64281cc5db07'
checks = []

def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)

def read(name):
    raw = (ROOT/'coefficients'/name).read_bytes()
    check(name+'_sha256', hashlib.sha256(raw).hexdigest() == manifest['sha256'][name])
    data = json.loads(raw)
    check(name+'_variables', data['variables'] == ['u','b','e'])
    return data['terms']

def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def univariate(terms):
    check('only_b_terms', all(e[0] == e[2] == 0 for e,c in terms))
    p = [Fraction(0)] * (max(e[1] for e,c in terms)+1)
    for e,c in terms:
        p[e[1]] += Fraction(c)
    return trim(p)

def mul(p,q):
    out = [Fraction(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):
            out[i+j] += a*b
    return trim(out)

def rem_mod(p,q,prime):
    p = trim(p)
    while p != [0] and len(p) >= len(q):
        k = len(p)-len(q)
        c = p[-1]*pow(q[-1],-1,prime) % prime
        for j,b in enumerate(q):
            p[k+j] = (p[k+j]-c*b) % prime
        p = trim(p)
    return p

modular_witnesses = {}
def certify_coprime(name,p,q):
    # A preserved-degree resultant nonzero in one good reduction is nonzero
    # over Q. This is an exact nonvanishing certificate, not probabilistic testing.
    for prime in (101,103,107,109,113,127,131,137,139,149):
        if any(c.denominator % prime == 0 for c in (*p,*q)):
            continue
        a = [c.numerator*pow(c.denominator,-1,prime)%prime for c in p]
        b = [c.numerator*pow(c.denominator,-1,prime)%prime for c in q]
        if a[-1] == 0 or b[-1] == 0:
            continue
        while b != [0]:
            a,b = b,rem_mod(a,b,prime)
        if len(a) == 1 and a[0] != 0:
            modular_witnesses[name] = prime
            check(name,True)
            return
    raise AssertionError('no retained good prime: '+name)

def deriv(p):
    return [i*p[i] for i in range(1,len(p))] or [Fraction(0)]

D = univariate(read('H_main_discriminant.json'))
factors = {i:univariate(read(f'support_factor_{i}.json')) for i in (0,1,2,3,7,8)}
expected_degrees = {0:1,1:2,2:2,3:2,7:44,8:91}
for i, expected in {0:[-29,23],1:[-75,-135,98],2:[297,-1436,784],3:[-342,-1910,1225]}.items():
    F = factors[i]
    ratio = F[-1]/expected[-1]
    check(f'F{i}_matches_volume_definition', F == [ratio*c for c in expected])
for i,F in factors.items():
    check(f'F{i}_degree', len(F)-1 == expected_degrees[i])
    certify_coprime(f'F{i}_squarefree', F,deriv(F))
for i,F in factors.items():
    for j,G in factors.items():
        if i < j:
            certify_coprime(f'F{i}_F{j}_coprime', F,G)
product = [Fraction(1)]
for i,power in ((0,2),(1,2),(2,2),(3,2),(7,1),(8,2)):
    for _ in range(power):
        product = mul(product,factors[i])
check('discriminant_degree_240', len(D)-1 == len(product)-1 == 240)
scale = D[-1]/product[-1]
check('discriminant_nonzero_rational_scale', scale != 0)
check('exact_discriminant_product', D == [scale*c for c in product])
H = read('H_factor_1.json')
check('H_e_degree_16', max(e[2] for e,c in H) == 16)
leading = [(e,c) for e,c in H if e[2] == 16]
check('H_e_leading_constant', len(leading) == 1 and leading[0][0] == [0,0,16]
      and Fraction(leading[0][1]) == 4056724031542677933477076950962364655200000000)
result = {
    'check_count':len(checks), 'checks':checks,
    'factorization':'D_H=c*F0^2*F1^2*F2^2*F3^2*F44*F91^2',
    'degree':240, 'F44_valuation':1, 'modular_nonvanishing_witnesses':modular_witnesses,
    'scale_numerator_digits':len(str(abs(scale.numerator))),
    'scale_denominator_digits':len(str(scale.denominator)),
    'scope':'exact coefficient checks; reuses the retained resultant; analytic local inverse and neighborhood theorem are proved in volume X section 13, without Lean verification'
}
(ROOT/'discriminant_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,sort_keys=True))
