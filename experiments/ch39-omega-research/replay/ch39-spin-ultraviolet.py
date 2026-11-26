"""Weak running above the matched spin-balanced crossing and a lower bound for a simple group."""
from pathlib import Path
import runpy
from flint import arb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/'ch39-spin-balanced-rge.py'))
q=m['q']; pi=arb.pi()
r=q('8.37369').union(q('8.53494'))
geo=m['m']['data']['m']['model']
s=geo['dirac_s'](r)
logs=geo['geometry'](r)
def fermion_budget(p):
    total=q(0)
    for j in range(1,101):
        # The function (1+j*s^2)/(1+s^2) is increasing in s for j>=1.
        mx=(-logs[0].lower()).exp()*q(j).sqrt()*(1+j*s.upper()**2)/(1+s.upper()**2)
        denominator=max(q(1),mx.upper())
        if denominator<p:total+=(p/denominator).log()
    return total
def weak(p):
    extra=18*q('10/3')*p.log()+q('8/3')*fermion_budget(p)
    return q('30.07107')+((q('19/6')+21/(8*pi))*p.log()-extra)/(2*pi)
def simple(p):
    positive=(q('25/2')+180)*p.log()+8*fermion_budget(p)
    return (positive-6*pi*q('30.07107'))/(11*p.log())
for ratio in (10,14,16,18,20,30,50,100):
    p=q(ratio)
    print('scale ratio',ratio,'weak inverse upper',weak(p).upper().str(10),
          'simple I*C_A necessary above',simple(p).lower().str(10))
assert weak(q(14))<1 and simple(q(20))>q('5.49734') and simple(q(30))>q('8.69866')
