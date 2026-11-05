#!/usr/bin/env python3
"""Exact finite checks for Krasula, arXiv:2609.15423v1, Section 5.1.

All arithmetic is in integers modulo a prime. Every vector subspace is
visited via its unique reduced row echelon basis. No random samples,
floating point operations, external packages, or Lean proofs are used.
"""
from itertools import combinations, product


def require(ok, message):
    if not ok:
        raise AssertionError(message)


class Algebra:
    def __init__(self, name, p, d, multiply, one, functional):
        self.name, self.p, self.d = name, p, d
        self.vectors = list(product(range(p), repeat=d))
        self.index = {v: i for i, v in enumerate(self.vectors)}
        self.elements = frozenset(range(len(self.vectors)))
        self.zero = 0
        self.one = self.index[one]
        self.basis = tuple(self.index[tuple(int(i == j) for i in range(d))]
                           for j in range(d))
        self.add = [[self.index[tuple((a + b) % p for a, b in zip(x, y))]
                     for y in self.vectors] for x in self.vectors]
        self.mul = [[self.index[multiply(x, y)] for y in self.vectors]
                    for x in self.vectors]
        self.functional = functional

    def span(self, basis):
        result = {0}
        for v in basis:
            new = set(result)
            w = 0
            for _ in range(1, self.p):
                w = self.add[w][v]
                new.update(self.add[x][w] for x in result)
            result = new
        return frozenset(result)

    def subspaces(self):
        """Exactly one RREF basis per vector subspace, including 0 and R."""
        for k in range(self.d + 1):
            for pivots in combinations(range(self.d), k):
                free = [(i, j) for i, pivot in enumerate(pivots)
                        for j in range(pivot + 1, self.d) if j not in pivots]
                for entries in product(range(self.p), repeat=len(free)):
                    rows = [[int(j == pivot) for j in range(self.d)]
                            for pivot in pivots]
                    for (i, j), a in zip(free, entries):
                        rows[i][j] = a
                    basis = tuple(self.index[tuple(row)] for row in rows)
                    yield self.span(basis), basis

    def left_ann(self, subset):
        return frozenset(a for a in self.elements
                         if all(self.mul[a][x] == 0 for x in subset))

    def right_ann(self, subset):
        return frozenset(a for a in self.elements
                         if all(self.mul[x][a] == 0 for x in subset))

    def validate(self):
        E, add, mul = self.elements, self.add, self.mul
        for x in E:
            require(mul[self.one][x] == x == mul[x][self.one], self.name + ': unit')
            for y in E:
                for z in E:
                    require(mul[mul[x][y]][z] == mul[x][mul[y][z]],
                            self.name + ': associativity')
                    require(mul[x][add[y][z]] == add[mul[x][y]][mul[x][z]],
                            self.name + ': left distributivity')
                    require(mul[add[x][y]][z] == add[mul[x][z]][mul[y][z]],
                            self.name + ': right distributivity')
        # Nondegenerate associative pairing lambda(xy) is an exact Frobenius
        # certificate for these finite-dimensional algebras over F_p.
        gram = [[self.functional(self.vectors[mul[x][y]]) % self.p
                 for y in self.basis] for x in self.basis]
        rank = 0
        for col in range(self.d):
            pivot = next((i for i in range(rank, self.d) if gram[i][col]), None)
            if pivot is None:
                continue
            gram[rank], gram[pivot] = gram[pivot], gram[rank]
            inv = pow(gram[rank][col], -1, self.p)
            gram[rank] = [(a * inv) % self.p for a in gram[rank]]
            for i in range(self.d):
                if i != rank:
                    a = gram[i][col]
                    gram[i] = [(x - a * y) % self.p
                               for x, y in zip(gram[i], gram[rank])]
            rank += 1
        require(rank == self.d, self.name + ': Frobenius functional')
        return rank


def truncated(p, d):
    def mul(x, y):
        return tuple(sum(x[i] * y[k-i] for i in range(k+1)) % p for k in range(d))
    return Algebra('F_%d[t]/(t^%d)' % (p, d), p, d, mul,
                   (1,) + (0,) * (d-1), lambda x: x[-1])


def matrix2(p):
    def mul(x, y):
        return tuple(sum(x[2*i+k] * y[2*k+j] for k in range(2)) % p
                     for i in range(2) for j in range(2))
    return Algebra('M_2(F_%d)' % p, p, 4, mul, (1, 0, 0, 1), lambda x: x[0]+x[3])


def cycle2(p):
    # Basis e_0,e_1,a_01,a_10 of the two-cycle path algebra modulo paths of length 2.
    def mul(x, y):
        return (x[0]*y[0] % p, x[1]*y[1] % p,
                (x[0]*y[2]+x[2]*y[1]) % p,
                (x[1]*y[3]+x[3]*y[0]) % p)
    return Algebra('two-cycle radical-square-zero/F_%d' % p, p, 4, mul,
                   (1, 1, 0, 0), lambda x: x[2]+x[3])


def direct_product(A, B):
    require(A.p == B.p, 'same prime field')
    def mul(x, y):
        a, b = x[:A.d], x[A.d:]
        c, d = y[:A.d], y[A.d:]
        return (A.vectors[A.mul[A.index[a]][A.index[c]]] +
                B.vectors[B.mul[B.index[b]][B.index[d]]])
    return Algebra('(%s) x (%s)' % (A.name, B.name), A.p, A.d+B.d, mul,
                   A.vectors[A.one]+B.vectors[B.one],
                   lambda x: A.functional(x[:A.d])+B.functional(x[A.d:]))


def check_algebra(A):
    rank = A.validate()
    subspaces = list(A.subspaces())
    require(len(subspaces) == len({s for s, _ in subspaces}), 'RREF uniqueness')
    def is_left(s, basis):
        return all(A.mul[a][x] in s for a in A.basis for x in basis)
    def is_right(s, basis):
        return all(A.mul[x][a] in s for a in A.basis for x in basis)
    left = [s for s, b in subspaces if is_left(s, b)]
    right = [s for s, b in subspaces if is_right(s, b)]
    two = [s for s in left if s in right]
    zero = frozenset({0})
    central = [e for e in A.elements if A.mul[e][e] == e and
               all(A.mul[e][x] == A.mul[x][e] for x in A.elements)]
    generated = {frozenset(A.mul[x][e] for x in A.elements): e for e in central}
    left_lcd = [I for I in left if I & A.right_ann(I) == zero]
    right_lcd = [J for J in right if J & A.left_ann(J) == zero]
    # Fidelity: Theorem 19, all LCD ideals of a finite Frobenius ring are
    # generated by central idempotents; also the size condition of Theorem 15.
    for I in left:
        require(A.left_ann(A.right_ann(I)) == I, 'left double annihilator')
        require(len(I)*len(A.right_ann(I)) == len(A.elements), 'source size condition')
    for J in right:
        require(A.right_ann(A.left_ann(J)) == J, 'right double annihilator')
    require(set(left_lcd) == set(generated) == set(right_lcd),
            A.name + ': central-idempotent LCD classification')
    for I in left_lcd:
        K = A.span(tuple(A.mul[x][a] for x in I for a in A.elements))
        require(K in two, 'IR is a two-sided ideal')
        require(A.left_ann(K) == A.left_ann(I), 'l(IR)=l(I)')
        require(A.right_ann(I) <= A.left_ann(I), 'LCD implies r(I)<=l(I)')
        require(A.left_ann(A.left_ann(K)) <= I <= K, 'contracting sandwich')
        require(K == I == A.left_ann(A.left_ann(I)), 'LCD two-sidedness bridge')
        require(A.left_ann(I) == A.right_ann(I), 'LCD annihilator equality')
        require(A.span(tuple(I | A.right_ann(I))) == A.elements, 'direct-sum spanning')
    # Test the orbit mechanism on every two-sided ideal, without asserting
    # that l^2 fixes all ideals (the non-symmetric example below rules that out).
    for K in two:
        T = lambda s: A.left_ann(A.left_ann(s))
        U = lambda s: A.right_ann(A.right_ann(s))
        require(T(K) in two and U(T(K)) == K and T(U(K)) == K, 'inverse orbit maps')
        if T(K) <= K:
            require(T(K) == K, 'Artinian contracting orbit')
    print('%s: |R|=%d, pairing_rank=%d, subspaces=%d, left=%d, right=%d, two-sided=%d, '
          'left_LCD=%d, right_LCD=%d, central_idempotents=%d; OK' %
          (A.name, len(A.elements), rank, len(subspaces), len(left), len(right), len(two),
           len(left_lcd), len(right_lcd), len(central)))
    return len(subspaces), len(left), len(right), len(left_lcd), len(right_lcd)


def negative_controls():
    A = cycle2(2)
    a = A.span((A.index[(0, 0, 1, 0)],))
    require(A.left_ann(a) != A.right_ann(a), 'must not assert l(I)=r(I) for arbitrary QF I')
    require(A.left_ann(A.left_ann(a)) != a, 'must not assert l^2 fixes every QF ideal')
    print('Sidedness control: in the two-cycle QF algebra, (a_01) has l(I) != r(I) '
          'and l^2(I) != I; OK')
    B = direct_product(truncated(2, 2), truncated(2, 1))
    I = frozenset(B.index[x+(0,)] for x in product(range(2), repeat=2))
    eps = B.index[(0, 1, 0)]
    require(I & B.right_ann(I) == frozenset({0}), 'product-factor ideal is LCD')
    require(eps != 0 and eps in I and B.mul[eps][eps] == 0, 'LCD can contain square-zero ideal')
    # The one-dimensional (eps,0) span is a nonzero square-zero ideal contained in I.
    N = B.span((eps,))
    require(all(B.mul[x][y] == 0 for x in N for y in N), 'N^2=0')
    require(all(B.mul[x][n] in N and B.mul[n][x] in N for x in B.elements for n in N),
            'N is two-sided')
    print('Semiprimeness control: (F_2[eps]/eps^2) x 0 is a proper LCD ideal '
          'containing a nonzero square-zero ideal; OK')


def main():
    print('Source: Dominik Krasula, Free and projective LCD codes, arXiv:2609.15423v1.')
    print('Fidelity: Theorem 19 (finite Frobenius LCD ideals) and Theorem 15 (size condition).')
    print('Settlement checks: both-sided classification and each proposed annihilator/orbit identity.')
    rings = [truncated(2, 1), truncated(2, 2), truncated(2, 3),
             truncated(3, 2), truncated(3, 3), matrix2(2), matrix2(3), cycle2(2), cycle2(3),
             direct_product(truncated(2, 2), truncated(2, 1)),
             direct_product(matrix2(2), truncated(2, 1)),
             direct_product(cycle2(2), truncated(2, 1))]
    totals = [0]*5
    for A in rings:
        values = check_algebra(A)
        totals = [x+y for x, y in zip(totals, values)]
    negative_controls()
    print('TOTAL rings=%d, subspaces=%d, left_ideals=%d, right_ideals=%d, '
          'left_LCD=%d, right_LCD=%d' % (len(rings), *totals))
    print('Scope: exhaustive for the listed finite algebras; no claim of an exhaustive '
          'list of QF rings or a proof for arbitrary/infinite QF rings.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
