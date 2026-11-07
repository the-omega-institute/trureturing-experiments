#!/usr/bin/env python3
"""Exact checks for Verbeken's odd case and the proposed characteristic-two family.

Only Python integer arithmetic is used. Finite checks do not certify the
unbounded family: its algebraic argument belongs to the preregistration.
"""
from itertools import product


def span(vectors):
    pivots = {}
    for v in vectors:
        while v:
            k = v.bit_length() - 1
            if k in pivots:
                v ^= pivots[k]
            else:
                pivots[k] = v
                break
    ans = {0}
    for v in pivots.values():
        ans |= {x ^ v for x in ans}
    return ans


class BinaryAlgebra:
    """J in F2[x,y]/(xy,x^(m+1)+y^(m+1)), in its normal-form basis."""
    def __init__(self, m):
        self.m = m
        self.dim = 2 * m + 1
        self.z = 1 << (2 * m)
        self.basis = [1 << i for i in range(self.dim)]
        self.x, self.y = self.basis[0], self.basis[m]
        self.xm, self.ym = self.basis[m-1], self.basis[2*m-1]
        self.bb = [[self.basis_product(i,j) for j in range(self.dim)]
                   for i in range(self.dim)]

    def basis_product(self, i, j):
        m = self.m
        if i == 2*m or j == 2*m or (i < m) != (j < m):
            return 0
        degree = i % m + j % m + 2
        if degree > m+1:
            return 0
        if degree == m+1:
            return self.z
        return 1 << ((0 if i < m else m) + degree - 1)

    def mul(self, a, b):
        out = 0
        while a:
            i = (a & -a).bit_length()-1
            a &= a-1
            bb = b
            while bb:
                j = (bb & -bb).bit_length()-1
                bb &= bb-1
                out ^= self.bb[i][j]
        return out

    def unit(self, a, b):
        # a denotes the principal unit 1+a.
        return a ^ b ^ self.mul(a,b)

    def upow(self, a, n):
        out = 0
        while n:
            if n & 1:
                out = self.unit(out,a)
            a = self.unit(a,a)
            n >>= 1
        return out

    def D(self, a):
        return (self.ym if a & self.x else 0) ^ (self.xm if a & self.y else 0)

    def phi(self, a):
        return a ^ self.D(a)

    def epsilon(self, a):
        # Triangular factorization into (1+x^i), (1+y^i), then (1+z).
        # Remove a factor by its finite geometric-series inverse.
        for b in self.basis[:-1]:
            if a & b:
                inv = 0
                power = b
                while power:
                    inv ^= power
                    power = self.mul(power,b)
                a = self.unit(a,inv)
        assert a in (0,self.z)
        return int(a == self.z)

    def star(self, a, b):
        return self.mul(a,b) ^ ((self.D(b) ^ self.mul(a,self.D(b)))
                                if self.epsilon(a) else 0)


def check_ring(A):
    B = A.basis
    # Associativity and Leibniz are multilinear identities, so check a basis.
    for a,b,c in product(B,repeat=3):
        assert A.mul(A.mul(a,b),c) == A.mul(a,A.mul(b,c))
    for a,b in product(B,repeat=2):
        assert A.D(A.mul(a,b)) == A.mul(A.D(a),b) ^ A.mul(a,A.D(b))
        assert A.mul(A.D(a),A.D(b)) == 0
    assert A.D(A.x) == A.ym and A.D(A.y) == A.xm
    assert A.D(A.x) != 0
    assert all(A.D(A.D(a)) == 0 for a in B)
    # Adjoining 1 adds D(1)=0; all extra unit Leibniz/associativity cases follow.
    current = set(B)
    for _ in range(A.m+1):
        current = {A.mul(a,b) for a in current for b in B} - {0}
    assert current == set()


def odd_source_fidelity():
    # Source section 3 and Proposition 4.1, p=3: dim J=7,
    # J^4=<z>, J^5=0 and U^3=1+<x^3,y^3>, of order 9.
    p,m = 3,3
    n = 2*m+1
    zero = (0,)*n
    def add(a,b):
        return tuple((x+y)%p for x,y in zip(a,b))
    def mul(a,b):
        out = [0]*n
        for i,x in enumerate(a):
            if not x: continue
            for j,y in enumerate(b):
                if not y or i == 2*m or j == 2*m or (i < m) != (j < m): continue
                d = i % m + j % m + 2
                if d <= m+1:
                    k = 2*m if d == m+1 else (0 if i < m else m)+d-1
                    out[k] = (out[k]+x*y)%p
        return tuple(out)
    basis = [tuple(int(i==j) for j in range(n)) for i in range(n)]
    elems = list(product(range(p),repeat=n))
    cubes = {mul(mul(a,a),a) for a in elems}
    W = {tuple(a if i==m-1 else b if i==2*m-1 else 0 for i in range(n))
         for a,b in product(range(p),repeat=2)}
    assert cubes == W and len(cubes)==9 and len(elems)==2187
    powers = set(basis)
    for k in range(2,6):
        powers = {mul(a,b) for a in powers for b in basis} - {zero}
        if k==4: assert powers == {basis[-1]}
        if k==5: assert powers == set()
    print('FIDELITY p=3: dim(B)=8; |J|=2187; J^4=<z>; J^5=0; |U^3|=9')


def principal():
    A = BinaryAlgebra(4)
    check_ring(A)
    N = 1 << A.dim
    E = list(range(N))
    eps = [A.epsilon(a) for a in E]
    def formula(a):
        aa = [(a>>i)&1 for i in range(4)]
        bb = [(a>>(4+i))&1 for i in range(4)]
        return ((a>>8) ^ aa[0]*aa[3] ^ (aa[0]^aa[1])*aa[2]
                ^ bb[0]*bb[3] ^ (bb[0]^bb[1])*bb[2]) & 1
    assert eps == [formula(a) for a in E]
    M = [[A.mul(a,b) for b in E] for a in E]
    U = [[a ^ b ^ M[a][b] for b in E] for a in E]
    ph = [A.phi(a) for a in E]
    circle = [[U[a][ph[b] if eps[a] else b] for b in E] for a in E]
    lam = [[circle[a][b]^a for b in E] for a in E]
    inv = []
    for a in E:
        assert circle[a][0] == a and circle[0][a] == a
        candidates = [b for b in E if circle[a][b]==0 and circle[b][a]==0]
        assert len(candidates)==1
        inv.append(candidates[0])
        assert A.D(A.D(a))==0 and M[a][A.D(a)]==0
        assert eps[ph[a]]==eps[a]
        # Every lambda is checked against the linear map specified by its basis.
        for b in E:
            assert lam[a][b] == (0 if b==0 else
                lam[a][b & (b-1)] ^ lam[a][b & -b])
    comms,fourths = set(),set()
    for a in E:
        sq = circle[a][a]
        fourths.add(circle[sq][sq])
        assert circle[sq][sq] == A.upow(a,4)
        for b in E:
            assert eps[U[a][b]] == eps[a]^eps[b]
            assert ph[U[a][b]] == U[ph[a]][ph[b]]
            c = circle[a][b]
            assert eps[c] == eps[a]^eps[b]
            # Equality of linear parts on a basis certifies affine composition.
            for v in A.basis:
                assert lam[c][v] == lam[a][lam[b][v]]
            comms.add(circle[circle[c][inv[a]]][inv[b]])
    V = span([A.xm,A.ym])
    assert comms == fourths == V and len(V)==4
    assert all(circle[a][b] in V for a,b in product(V,repeat=2))
    assert all(circle[v][a] == circle[a][v] for v in V for a in E)
    assert all(circle[circle[a][a]][circle[a][a]] in V for a in E)
    assert all(circle[A.upow(a,4)][A.upow(a,4)]==0 for a in E)
    assert any(circle[a][a] != 0 and circle[circle[a][a]][circle[a][a]] != 0 for a in E)
    T = span([A.xm,A.ym,A.z])
    star = [[circle[a][b]^a^b for b in E] for a in E]
    assert span(star[t][b] for t in T for b in E)==T
    assert star[A.xm][A.x]==A.z
    assert star[A.z][A.x]==A.ym and star[A.z][A.y]==A.xm
    R = set(E)
    dimensions = [A.dim]
    while True:
        new = span(star[a][b] for a in R for b in E)
        dimensions.append(len(new).bit_length()-1)
        assert T <= new
        if new==R: break
        R = new
    assert R==T
    print('PRINCIPAL m=4: |J|=512; algebra and nonzero derivation hypotheses OK')
    print('PRINCIPAL epsilon: explicit coefficient formula agrees on all 512 elements')
    print('PRINCIPAL all 262144 pairs: character, phi, graph closure and affine composition OK')
    print('PRINCIPAL group: |G\'|=|G^4|=4; G\'=G^4; class=2; exponent=8')
    print('PRINCIPAL right-series dimensions:', ','.join(map(str,dimensions)))
    print('PRINCIPAL inductive certificate: dim(T)=3; T star J=T; T is nonzero')


def family_checks():
    for r in (2,3):
        A = BinaryAlgebra(4*r)
        check_ring(A)
        assert A.upow(A.basis[r-1],4)==A.xm
        assert A.upow(A.basis[A.m+r-1],4)==A.ym
        state = 43
        mask = (1<<A.dim)-1
        for _ in range(4096):
            state = (1103515245*state+12345)&0x7fffffff
            a = state & mask
            state = (1103515245*state+12345)&0x7fffffff
            b = state & mask
            assert A.epsilon(A.unit(a,b))==A.epsilon(a)^A.epsilon(b)
            assert A.epsilon(A.phi(a))==A.epsilon(a)
            assert A.mul(a,A.D(a))==0
            assert A.D(A.D(a))==0
            assert A.phi(A.unit(a,b))==A.unit(A.phi(a),A.phi(b))
            def circ(x,y): return A.unit(x,A.phi(y) if A.epsilon(x) else y)
            sq = circ(a,a)
            assert circ(sq,sq)==A.upow(a,4)
        T = span([A.xm,A.ym,A.z])
        assert span(A.star(t,b) for t in T for b in A.basis)==T
        print(f'FAMILY r={r}, m={4*r}, |J|=2^{A.dim}: basis identities and 4096 exact pairs OK')


if __name__ == '__main__':
    odd_source_fidelity()
    principal()
    family_checks()
    print('ALL_CHECKS_OK')
