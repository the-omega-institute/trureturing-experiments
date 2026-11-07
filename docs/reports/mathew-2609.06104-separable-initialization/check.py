#!/usr/bin/env python3
"""Exact Pauli-algebra checks for the literal complete-graph/Jordan-Wigner family.

Coefficients are SymPy exact expressions. A word (x,z) means X^x Z^z,
so multiplication contributes (-1)^(z dot x'). No floating-point arithmetic.
Finite checks support a preregistration; they do not establish the all-size claim.
"""
from itertools import combinations, product
import sympy as s

I = s.I
ONE = {(0, 0): s.Integer(1)}


def clean(a):
    return {k: s.expand(v) for k, v in a.items() if s.expand(v) != 0}


def add(*aa):
    out = {}
    for a in aa:
        for k, v in a.items():
            out[k] = out.get(k, 0) + v
    return clean(out)


def scale(c, a):
    return clean({k: c * v for k, v in a.items()})


def mul(a, b):
    out = {}
    for (x, z), u in a.items():
        for (xx, zz), v in b.items():
            key = (x ^ xx, z ^ zz)
            sign = -1 if (z & xx).bit_count() % 2 else 1
            out[key] = out.get(key, 0) + sign * u * v
    return clean(out)


def adj(a):
    return clean({(x,z): s.conjugate(v) * (-1 if (x & z).bit_count() % 2 else 1)
                  for (x,z),v in a.items()})


def comm(a, b):
    return add(mul(a,b), scale(-1,mul(b,a)))


def anti(a, b):
    return add(mul(a,b), mul(b,a))


def equal(a, b):
    assert not add(a, scale(-1,b))


def tr(a, n):
    return s.expand(2**n * a.get((0,0),0))


def tensor(a, b, n):
    out = {}
    for (x,z),u in a.items():
        for (xx,zz),v in b.items():
            out[(x | (xx << n), z | (zz << n))] = u*v
    return clean(out)


def word(labels):
    x = z = y = 0
    for q, label in enumerate(labels):
        if label in 'XY': x |= 1 << q
        if label in 'ZY': z |= 1 << q
        if label == 'Y': y += 1
    return {(x,z): I**y}


def model(m):
    edges = list(combinations(range(2*m),2))
    nb = len(edges)
    nq = (m+2)//2
    pp = [word(['X' if a==v else 'Z' if b==v else 'I' for a,b in edges])
          for v in range(2*m)]
    yy = [word(['Y' if e==(2*r,2*r+1) else 'I' for e in edges]) for r in range(m)]
    cc = [word(['Z' if q<j//2 else ('X' if j%2==0 else 'Y') if q==j//2 else 'I'
                for q in range(nq)]) for j in range(m+1)]
    aa = [scale(I,mul(cc[0],cc[r+1])) for r in range(m)]
    hh = [scale(1/s.sqrt(2),add(pp[2*r],pp[2*r+1])) for r in range(m)]
    kk = [scale(I,mul(yy[r],hh[r])) for r in range(m)]
    hb = add(*pp)
    c = add(*(tensor(yy[r],aa[r],nb) for r in range(m)))
    d = add(*(tensor(kk[r],aa[r],nb) for r in range(m)))
    return nb,nq,pp,yy,cc,aa,hh,kk,hb,c,d


def check(m):
    nb,nq,pp,yy,cc,aa,hh,kk,hb,c,d = model(m)
    # Fidelity: the source states (HB^NC)^2=M I and HC^2=m I.
    equal(mul(hb,hb),scale(2*m,ONE))
    equal(mul(c,c),scale(m,ONE))
    for family in (pp,cc,aa,hh,kk):
        for r,a in enumerate(family):
            equal(adj(a),a)
            equal(mul(a,a),ONE)
            for b in family[r+1:]: equal(anti(a,b),{})
    power_op = scale(I,comm(c,tensor(hb,ONE,nb)))
    equal(power_op,scale(2*s.sqrt(2),d))
    equal(anti(c,d),{})
    equal(comm(c,mul(d,d)),{})
    # Literal D^2 decomposition consumed by the covariance argument.
    expansion = scale(m,ONE)
    for r in range(m):
        for q in range(m):
            if r != q:
                expansion=add(expansion,tensor(mul(kk[r],kk[q]),mul(aa[r],aa[q]),nb))
    equal(mul(d,d),expansion)
    # Ground-supported product attainment from commuting paired involutions.
    ground = scale(s.Rational(1,2),add(ONE,scale(-1/s.sqrt(2*m),hb)))
    equal(mul(ground,ground),ground)
    assert tr(ground,nb)==2**(nb-1)
    jj=[scale(I,mul(kk[2*j],kk[2*j+1])) for j in range(m//2)]
    tt=[scale(I,mul(aa[2*j],aa[2*j+1])) for j in range(m//2)]
    for j in jj:
        equal(adj(j),j);equal(mul(j,j),ONE);equal(comm(j,ground),{})
    for a,b in combinations(jj,2):equal(comm(a,b),{})
    for signs in product((-1,1), repeat=m//2):
        proj=ground
        for sign,j in zip(signs,jj):proj=mul(proj,scale(s.Rational(1,2),add(ONE,scale(sign,j))))
        weight=tr(proj,nb)
        if weight != 0:break
    assert weight > 0
    equal(mul(proj,proj),proj)
    tau=scale(1/weight,proj)
    aux=ONE
    for sign,t in zip(signs,tt):aux=mul(aux,scale(s.Rational(1,2),add(ONE,scale(-sign,t))))
    equal(mul(aux,aux),aux)
    assert tr(aux,nq)>0
    sigma=scale(1/tr(aux,nq),aux)
    equal(adj(tau),tau);equal(adj(sigma),sigma)
    assert tr(tau,nb)==tr(sigma,nq)==1
    covariance_inner=s.Integer(0)
    for r in range(m):
        for q in range(m):
            if r != q:
                gb=s.expand(I*tr(mul(tau,mul(kk[r],kk[q])),nb))
                ga=s.expand(I*tr(mul(sigma,mul(aa[r],aa[q])),nq))
                covariance_inner += gb*ga
    ell=2*(m//2)
    assert s.simplify(covariance_inner+ell)==0
    rho=tensor(tau,sigma,nb)
    second=tr(mul(rho,mul(d,d)),nb+nq)
    assert s.simplify(second-(m+ell))==0
    # U at the advertised peak time pi/(4 sqrt(m)).
    u=scale(1/s.sqrt(2),add(ONE,scale(-I/s.sqrt(m),c)))
    equal(mul(adj(u),u),ONE)
    evolved=mul(mul(u,rho),adj(u))
    peak=s.simplify(tr(mul(evolved,power_op),nb+nq))
    target=2*s.sqrt(2)*(1+s.Rational(ell,m))
    assert s.simplify(peak-target)==0
    assert s.simplify(tr(mul(evolved,mul(d,d)),nb+nq)-second)==0
    assert m+ell <= 2*m
    print(f'm={m}; N={nb}; auxiliary_qubits={nq}; source_squares=OK; '
          f'Clifford_and_conserved_square=OK; covariance_inner={covariance_inner}; '
          f'ground_product_peak={peak}; separable_bound_squared={8*(m+ell)}; benchmark_squared={16*m}')


if __name__=='__main__':
    for size in range(1,5):check(size)
    print('ALL_CHECKS_OK')
