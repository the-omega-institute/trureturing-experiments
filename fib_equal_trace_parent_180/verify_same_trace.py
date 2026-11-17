#!/usr/bin/env python3
"""Independent SymPy construction and exact check of the FLINT certificate.

Use an owned parent process with a 60-second timeout. No repository access.
The arithmetic check is independent implementation, not independent peer review.
"""
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import sympy as s

ROOT=Path(__file__).resolve().parent


def read_poly(data, syms):
    return s.Poly.from_dict({tuple(m):s.Rational(co) for m,co in data['terms']},syms,domain=s.QQ)


def main():
    t=time.monotonic()
    u,b,e,j,theta,xi,n,a=s.symbols('u b e j theta xi n a')
    allsyms=(u,b,e,j,theta,xi)
    h=u*(u+1);ell=u*u+3*u-4
    al=3*h-10;be=8*u*u+18*u-30
    phi=al*a+be*n-23*h
    delta=n*h+a*n*ell-4*a*h
    psi6=(n+2*a)*(10*u*u+19*u-36)-10*a*ell-29*h
    psi7=4*(2*n+3*a)*(20*u*u+34*u-70)-105*a*ell-381*h
    M=s.Matrix([[s.diff(f,a),s.diff(f,n),f.subs({a:0,n:0})] for f in (phi,psi6,psi7)])
    K=205*u**4+313*u**3-1211*u**2-894*u+1940
    assert s.expand(M.det()+20*h*K)==0
    assert s.gcd(s.Poly(K,u),s.Poly(al,u)).degree()==0
    L=s.Poly(psi6-b*phi,a,n)
    A=L.coeff_monomial(a);B=L.coeff_monomial(n);C=L.coeff_monomial(1)
    Q=a*a-a*n-2*(u+1)*a-n*n+2*(2*u+1)*n-u*u+u+1
    T=(theta*n-a)*phi-5*xi*delta
    # Homogenize in a,n to substitute the line without rational functions.
    Qh=s.Poly(Q,a,n)
    Th=s.Poly(T,a,n)
    def restrict(p):
        return s.Poly(s.expand(sum(co*(-B*n-C)**i*(A*n)**k*A**(2-i-k)
                  for (i,k),co in p.terms())),n)
    qa=restrict(Qh);ta=restrict(Th)
    q2,q1,q0=[qa.nth(i) for i in (2,1,0)]
    t2,t1,t0=[ta.nth(i) for i in (2,1,0)]
    syms=(u,b,theta,xi)
    def divexact(p):
        q,r=s.div(s.Poly(p,*syms),s.Poly(A,*syms))
        assert r.is_zero
        return q
    N=divexact(q2*t0-t2*q0)
    D=divexact(t2*q1-q2*t1)
    data=json.loads((ROOT/'same_trace_elimination.json').read_text())
    inputs=json.loads((ROOT/'universal_0_input.json').read_text())
    for name,p in [('Nd',N),('Dd',D)]:
        expected=read_poly(inputs['polynomials'][name],allsyms).as_expr()
        assert s.Poly(expected,*syms)==p
    print(json.dumps({'stage':'independent_source_and_chart_identities','seconds':time.monotonic()-t}),flush=True)
    t1=time.monotonic()
    # An independent prime-field determinant witnesses nonvanishing of the
    # universal eliminant. Neither input drops its b degree at this point.
    Ns=s.Poly(N.as_expr().subs({theta:0,xi:1}),b,u,modulus=103)
    Ds=s.Poly(D.as_expr().subs({theta:0,xi:1}),b,u,modulus=103)
    assert Ns.degree(b)==3 and Ds.degree(b)==3
    Ss=Ns.resultant(Ds)
    Ks=s.Poly(K,u,modulus=103)
    residue=int(Ks.resultant(Ss))%103
    assert residue==87
    assert Ss.degree()==55
    z=read_poly(data['parent_polynomial'],allsyms)
    zz=z.as_expr().subs({theta:0,xi:1})
    assert int(zz)%103==residue
    zdiv,zrem=s.div(z,s.Poly(xi**4,*allsyms))
    assert zrem.is_zero and zdiv.total_degree()==20
    assert zdiv.degree(theta)==20 and zdiv.degree(xi)==20
    print(json.dumps({'stage':'independent_nonvanishing','residue_mod_103':residue,
          'specialized_S_u_degree':Ss.degree(),'seconds':time.monotonic()-t1}),flush=True)
    t2=time.monotonic()
    # Full exact reconstruction of the two eliminants, not sample equality.
    S=s.Poly(N.as_expr(),b,u,theta,xi).resultant(s.Poly(D.as_expr(),b,u,theta,xi))
    expected=read_poly(data['b_resultant'],allsyms).as_expr()
    assert S==s.Poly(expected,u,theta,xi)
    print(json.dumps({'stage':'independent_full_b_resultant','seconds':time.monotonic()-t2}),flush=True)
    t3=time.monotonic()
    Z=s.Poly(K,u,theta,xi).resultant(S)
    assert Z==s.Poly(z.as_expr(),theta,xi)
    print(json.dumps({'stage':'independent_full_parent_resultant','seconds':time.monotonic()-t3}),flush=True)
    q,r=s.div(Z,s.Poly(xi**4,theta,xi))
    assert r.is_zero
    # Primitive integer normalization of the delivered parent constraint.
    _,f=q.clear_denoms(convert=True)
    _,f=f.primitive()
    if f.LC()<0: f=-f
    factors=[]
    for item in data['parent_factors']:
        ff=read_poly(item['polynomial'],allsyms)
        if ff.as_expr()==xi: continue
        assert item['multiplicity']==1
        factors.append(s.Poly(ff.as_expr(),theta,xi))
    degrees=sorted(ff.total_degree() for ff in factors)
    assert degrees==[4,8,8]
    product=s.Poly(1,theta,xi)
    for ff in factors: product*=ff
    _,product=product.primitive()
    if product.LC()<0: product=-product
    assert product==f
    out={"variables":["theta","xi"],"terms":[[list(m),str(co)] for m,co in f.terms()]}
    encoded=(json.dumps(out,sort_keys=True,indent=2)+'\n').encode()
    coefficient_sha256=hashlib.sha256(encoded).hexdigest()
    assert coefficient_sha256=='0fae329c1795f14dc7ce1374d5c67c4313a35229af9de7f2fc7ea0afcf7c4c51'
    assert (ROOT/'same_trace_parent_F.json').read_bytes()==encoded
    (ROOT/'same_trace_parent_F.json').write_bytes(encoded)
    usage=resource.getrusage(resource.RUSAGE_SELF)
    result={'status':'passed','kind':'independent exact arithmetic implementation',
       'not_peer_review':True,'source_revision':'daec012452c47287e1a48afc68988b73e8ba9642',
       'quartic_identity':True,'chart_identities':True,'both_full_resultants_equal':True,
       'nonzero_residue_mod_103':residue,'nonzero_residue_parent':[0,1],
       'specialized_N_b_degree':Ns.degree(b),'specialized_D_b_degree':Ds.degree(b),
       'specialized_S_u_degree':Ss.degree(),
       'K_alpha_resultant':int(s.resultant(K,al,u)),
       'F_total_degree':f.total_degree(),'F_degree_theta':f.degree(theta),'F_degree_xi':f.degree(xi),
       'F_terms':len(f.terms()),'F_factors_degrees':degrees,
       'wall_seconds':time.monotonic()-t,'user_seconds':usage.ru_utime,
       'system_seconds':usage.ru_stime,
       'max_rss_bytes':usage.ru_maxrss if sys.platform=='darwin' else usage.ru_maxrss*1024,
       'F_coefficient_sha256':coefficient_sha256,
       'sympy_version':s.__version__,'python':sys.version.split()[0]}
    (ROOT/'same_trace_verification.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
