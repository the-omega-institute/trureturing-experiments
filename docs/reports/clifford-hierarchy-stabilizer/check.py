#!/usr/bin/env python3
"""Exact finite checks for Fechisin et al., arXiv:2610.00464v1, Appendix B.

No floating-point arithmetic is used. Phases are integer residues modulo Q;
Hermite normal forms compare the entire finite phase groups, not samples.
"""
from itertools import product
from math import comb, prod
from sympy import Matrix, eye
from sympy.matrices.normalforms import hermite_normal_form


def points(p, n):
    return list(product(range(p), repeat=n))


def add(x, y, p):
    return tuple((a + b) % p for a, b in zip(x, y))


def delta(f, h, xs, index, p, q):
    return tuple((f[index[add(x, h, p)]] - f[i]) % q
                 for i, x in enumerate(xs))


def source_fidelity():
    # B.3: the unextended CCZ state has only the identity onsite stabilizer.
    xs = points(2, 3)
    phase = {x: 2 * prod(x) for x in xs}  # fourth roots, CCZ = i^(2xyz)
    stabilizers = []
    for v in xs:
        for a, b, c, d in product(range(4), repeat=4):
            if all((phase[add(x, v, 2)] - phase[x]
                    - a*x[0] - b*x[1] - c*x[2] - d) % 4 == 0 for x in xs):
                stabilizers.append((v, a, b, c, d))
    assert stabilizers == [((0, 0, 0), 0, 0, 0, 0)]
    # B.4: literal six-qubit generators, encoded as monomial matrices.
    ys = points(2, 6)
    shifts = [(1,0,0,1,1,0), (0,1,0,1,0,1), (0,0,1,0,1,1)]
    def exponent(j, y):
        return [(y[1]+y[2]-y[5]), (y[0]+y[2]-y[4]),
                (y[0]+y[1]-y[3])][j] % 4
    def action(j, y):
        return add(y, shifts[j], 2), exponent(j, y)
    for j in range(3):
        square_nonidentity = False
        for y in ys:
            z, a = action(j, y); w, b = action(j, z)
            square_nonidentity |= (w != y or (a+b) % 4 != 0)
            u, c = action(j, w); t, d = action(j, u)
            assert t == y and (a+b+c+d) % 4 == 0
        assert square_nonidentity
    noncommuting_pairs = 0
    for j in range(3):
        for k in range(j):
            noncommuting_pairs += any(
                (exponent(k,y)+exponent(j,add(y,shifts[k],2))
                 -exponent(j,y)-exponent(k,add(y,shifts[j],2))) % 4
                for y in ys)
    assert noncommuting_pairs == 3
    support = {(*x, (x[0]+x[1])%2, (x[0]+x[2])%2,
                (x[1]+x[2])%2): 2*prod(x) for x in xs}
    for j in range(3):
        for y, a in support.items():
            z, b = action(j, y)
            assert z in support and (a+b-support[z]) % 4 == 0
    # Generator squares force precisely this eight-point support. Translation
    # edges connect it; invariant amplitudes therefore have one free scalar.
    surviving = []
    for y in ys:
        if all((exponent(j,y)+exponent(j,add(y,shifts[j],2)))%4 == 0
               for j in range(3)):
            surviving.append(y)
    assert set(surviving) == set(support)
    seen = {surviving[0]}
    pending = list(seen)
    while pending:
        y = pending.pop()
        for v in shifts:
            z = add(y,v,2)
            if z not in seen: seen.add(z); pending.append(z)
    assert seen == set(support)
    for x,y in product(range(2),repeat=2):
        assert (2*x*y - x-y+(x+y)%2) % 4 == 0  # B15
    print('FIDELITY: B.3 identity-only; B.4 order four, 3 noncommuting pairs,')
    print('          eight-point support, one invariant ray; B15 exact for all bits')


def compositions(total, n):
    if n == 1:
        yield (total,)
    else:
        for a in range(total+1):
            for tail in compositions(total-a,n-1): yield (a,*tail)


def difference_degree(f, p, n, r, q):
    xs = points(p,n); idx = {x:i for i,x in enumerate(xs)}
    # All (r+1)-fold coordinate differences, at every base point.
    for beta in compositions(r+1,n):
        g = tuple(f)
        for i, multiplicity in enumerate(beta):
            h = tuple(int(j==i) for j in range(n))
            for _ in range(multiplicity): g = delta(g,h,xs,idx,p,q)
        assert not any(g), (p,n,r,beta)


def phase_group_check(p,n,r):
    xs=points(p,n); q=p**(1+(r-1)//(p-1))
    canonical=[]
    for alpha in product(range(p),repeat=n):
        weight=sum(alpha)
        if not weight: continue
        for m in range(1,1+(r-weight)//(p-1)+1):
            f=[(q//p**m * prod(xi**ai for xi,ai in zip(x,alpha)))%q
               for x in xs]
            difference_degree(f,p,n,r,q)
            canonical.append(f)
    ridges=[]
    for b in xs[1:]:
        for a in range(1,p):
            for m in range(1,1+(r-a)//(p-1)+1):
                f=[q//p**m * (sum(bi*xi for bi,xi in zip(b,x))%p)**a%q
                   for x in xs]
                difference_degree(f,p,n,r,q)
                ridges.append(f)
    constants=Matrix.ones(len(xs),1)
    modulus=q*eye(len(xs))
    canon=Matrix.hstack(modulus,constants,Matrix(canonical).T)
    ridge=Matrix.hstack(modulus,constants,Matrix(ridges).T)
    assert hermite_normal_form(canon) == hermite_normal_form(ridge), (p,n,r)
    # Independently check the proposed tensor rank over F_p, avoiding r!.
    betas=list(compositions(r,n))
    rows=[[prod(bi**be for bi,be in zip(b,beta))%p for b in xs]
          for beta in betas]
    rank=0
    for col in range(len(xs)):
        pivot=next((j for j in range(rank,len(rows)) if rows[j][col]),None)
        if pivot is None: continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        inv=pow(rows[rank][col],-1,p)
        rows[rank]=[v*inv%p for v in rows[rank]]
        for j in range(len(rows)):
            if j!=rank:
                v=rows[j][col]
                rows[j]=[(x-v*y)%p for x,y in zip(rows[j],rows[rank])]
        rank+=1
    expected=sum(0<sum(a)<=r and (sum(a)-r)%(p-1)==0
                 for a in product(range(p),repeat=n))
    assert rank==expected,(p,n,r,rank,expected)
    a=1+(r-1)%(p-1); m=1+(r-a)//(p-1); denominator=p**m
    f=tuple(t**a for t in range(p)); one_xs=[(t,) for t in range(p)]
    index={x:i for i,x in enumerate(one_xs)}
    for _ in range(r):f=delta(f,(1,),one_xs,index,p,denominator)
    assert len(set(f))==1 and f[0]%(denominator//p)==0
    assert f[0]//(denominator//p)%p != 0
    print(f'SPAN p={p} ell={n} r={r}: full phase-group HNF equal; tensor rank={rank}; top derivative nonzero')


def main():
    source_fidelity()
    cases=0
    for p,n,limit in [(2,2,7),(2,3,6),(3,2,6),(3,3,4),(5,2,6)]:
        for r in range(1,limit+1):phase_group_check(p,n,r);cases+=1
    print(f'EXACT_CASES={cases}; every phase in the canonical finite groups is covered, modulo constants')
    print('BOUNDARY: finite checks support the route; they do not prove B.5 for all p, ell, k')
    print('ALL_CHECKS_OK')

if __name__=='__main__': main()
