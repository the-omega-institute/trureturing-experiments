#!/usr/bin/env python3
"""Exact reduced-incidence calculations; no repository writes or numerical claims.

Run from any cwd with Python >=3.12 and python-flint 0.9.
Each command runs one finite stage. A parent process must impose its deadline.
The fixed-parent calculation is a trial, not a generic certificate.
"""
import argparse
import json
import resource
import sys
import time
from pathlib import Path

import flint

ROOT = Path(__file__).resolve().parent


def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, sort_keys=True, indent=2) + "\n")


def encode(p):
    return {"variables": list(p.context().names()),
            "terms": [[[int(x) for x in m], str(c)] for m, c in p.terms()]}


def info(p):
    return {"degrees": [int(x) for x in p.degrees()], "total_degree": int(p.total_degree()),
            "terms": len(p.to_dict())}


def model(ctx, theta=None, xi=None):
    u, b, e, j = ctx.gens()[:4]
    h = u*(u+1)
    ell = u*u+3*u-4
    al = 3*h-10
    be = 8*u*u+18*u-30
    A = (10-3*b)*u*u+(8-3*b)*u+10*b-32
    B = (10-8*b)*u*u+(19-18*b)*u+30*b-36
    C = (23*b-29)*h
    if theta is None: theta = ctx.constant(3)/7
    if xi is None: xi = -ctx.constant(35)/29
    # Substitute a=(-Bn-C)/A and clear A^2 in the original quadrics.
    q2 = B*B+A*B-A*A
    q1 = 2*B*C+A*C+2*(u+1)*A*B+2*(2*u+1)*A*A
    q0 = C*C+2*(u+1)*A*C+(-u*u+u+1)*A*A
    t2 = -al*B*B-(theta*al-be-5*xi*ell)*A*B+theta*be*A*A
    t1 = -2*al*B*C-(theta*al-be-5*xi*ell)*A*C-(23+20*xi)*h*A*B-(23*theta+5*xi)*h*A*A
    t0 = -al*C*C-(23+20*xi)*h*A*C
    # Exact divisions are polynomial operations, not specialization assumptions.
    R = ((q2*t0-q0*t2)**2-(q2*t1-q1*t2)*(q1*t0-q0*t1))/(A**4)
    Nd = (q2*t0-t2*q0)/A
    Dd = (t2*q1-q2*t1)/A
    # n=Nd/Dd; a=Ad/Dd; A divides this numerator identically.
    Ad = -(B*Nd+C*Dd)/A
    phi = al*Ad+be*Nd-23*h*Dd
    eta = (4*(2*Nd+3*Ad)*(20*u*u+34*u-70)-105*Ad*ell-381*h*Dd)
    E = e*phi-eta
    jd = Nd*phi
    jn = 120*h*xi*Nd*Dd+(4*Nd+(1-u)*Dd)*(210*phi+120*xi*(Nd*(20*u*u+34*u-70)-35*h*Dd))
    J = j*jd-jn
    return locals()


def decode(ctx, data):
    coef = int if isinstance(ctx,flint.nmod_mpoly_ctx) else flint.fmpq
    return ctx.from_dict({tuple(m): coef(c) for m,c in data['terms']})


def stage_singular(modulus):
    ctx = (flint.nmod_mpoly_ctx.get(('u','b','e','j'), modulus) if modulus
           else flint.fmpq_mpoly_ctx.get(('u','b','e','j')))
    supports = []
    results = {}
    for label,ind in [('E',2),('J',3)]:
        data = json.loads((ROOT/f'fixed_{modulus}_{label}.json').read_text())
        images = [d for d in data['factors'] if d['info']['degrees'][ind]>0]
        assert len(images)==1 and images[0]['multiplicity']==1
        H = decode(ctx,images[0]['polynomial'])
        t = time.monotonic()
        Dy = H.resultant(H.derivative(ind),ind)
        Db = H.resultant(H.derivative(1),ind)
        G = Dy.gcd(Db)
        co,fs=G.factor()
        support=ctx.constant(1)
        for f,m in fs: support*=f
        results[label]={'seconds':time.monotonic()-t,'image':info(H),
            'discriminant':info(Dy),'horizontal_resultant':info(Db),
            'support':encode(support),'factors':[{'info':info(f),'multiplicity':m,'polynomial':encode(f)} for f,m in fs]}
        supports.append(support)
        print(json.dumps({'stage':'singular','label':label,'seconds':results[label]['seconds'],
             'degrees':{k:info(v) for k,v in [('Dy',Dy),('Db',Db),('G',G)]},
             'factor_degrees':[(info(f),m) for f,m in fs]}),flush=True)
    common=supports[0].gcd(supports[1])
    co,fs=common.factor()
    results['common']={'polynomial':encode(common),'factors':[{'polynomial':encode(f),'info':info(f),'multiplicity':m} for f,m in fs]}
    dump(f'fixed_{modulus}_singular.json',results)
    print(json.dumps({'stage':'common_singular_support','info':info(common),'factors':[(str(f),m) for f,m in fs]}),flush=True)


def stage_fixed(modulus):
    if modulus:
        ctx = flint.nmod_mpoly_ctx.get(('u','b','e','j'), modulus)
    else:
        ctx = flint.fmpq_mpoly_ctx.get(('u','b','e','j'))
    t = time.monotonic()
    m = model(ctx)
    measurements = {x: info(m[x]) for x in ('A','Nd','Dd','Ad','R','phi','eta','E','J')}
    measurements['setup_seconds'] = time.monotonic()-t
    print(json.dumps({'stage':'setup','modulus':modulus,'measurements':measurements}),flush=True)
    dump(f'fixed_{modulus}_input.json',{'modulus':modulus,'measurements':measurements,
         'polynomials':{x:encode(m[x]) for x in ('R','E','J','A','Dd','phi')}})
    for label in ('E','J'):
        t = time.monotonic()
        res = m['R'].resultant(m[label], 0)
        elapsed = time.monotonic()-t
        print(json.dumps({'stage':'resultant','label':label,'seconds':elapsed,'info':info(res)}),flush=True)
        t = time.monotonic()
        coeff, factors = res.factor()
        factor_data = [{'multiplicity':v,'info':info(p),'polynomial':encode(p)} for p,v in factors]
        dump(f'fixed_{modulus}_{label}.json', {'resultant_seconds':elapsed,'factor_seconds':time.monotonic()-t,
             'resultant_info':info(res),'constant':str(coeff),'factors':factor_data})
        print(json.dumps({'stage':'factors','label':label,'factor_seconds':time.monotonic()-t,
             'factors':[{k:v for k,v in d.items() if k!='polynomial'} for d in factor_data]}),flush=True)


def stage_universal(modulus):
    names=('u','b','e','j','theta','xi')
    ctx = (flint.nmod_mpoly_ctx.get(names,modulus) if modulus
           else flint.fmpq_mpoly_ctx.get(names))
    t=time.monotonic()
    m=model(ctx,ctx.gens()[4],ctx.gens()[5])
    measures={x:info(m[x]) for x in ('R','E','J','Nd','Dd','Ad','phi')}
    dump(f'universal_{modulus}_input.json',{'modulus':modulus,'measurements':measures,
        'polynomials':{x:encode(m[x]) for x in ('R','E','J','Nd','Dd','Ad','phi')}})
    print(json.dumps({'stage':'universal_setup','seconds':time.monotonic()-t,'inputs':measures}),flush=True)
    t=time.monotonic()
    res=m['R'].resultant(m['E'],0)
    elapsed=time.monotonic()-t
    print(json.dumps({'stage':'universal_resultant','seconds':elapsed,'info':info(res)}),flush=True)
    dump(f'universal_{modulus}_E_raw.json',{'modulus':modulus,'seconds':elapsed,'polynomial':encode(res),'info':info(res)})
    t=time.monotonic()
    co,fs=res.factor()
    factors=[{'multiplicity':v,'info':info(p),'polynomial':encode(p)} for p,v in fs]
    dump(f'universal_{modulus}_E.json',{'constant':str(co),'factor_seconds':time.monotonic()-t,'factors':factors})
    print(json.dumps({'stage':'universal_factors','seconds':time.monotonic()-t,'factors':[{k:v for k,v in d.items() if k!='polynomial'} for d in factors]}),flush=True)


if __name__ == '__main__':
    assert sys.version_info >= (3,12)
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=('fixed','singular','universal'))
    parser.add_argument('--modulus',type=int,default=0)
    args = parser.parse_args()
    {'fixed':stage_fixed,'singular':stage_singular,'universal':stage_universal}[args.stage](args.modulus)
