#!/usr/bin/env python3
"""Exact small-case checks for Wu--Zhou, arXiv:2609.00533v1, Section 8(iv).

Construction uses increasing coordinate indices, nonzero minors and Cramer's
rule. Gaussian elimination is used only by the independent rank checker.
No ordering of field elements is used by the map. Integer labels in enumeration
are a test harness, not an input to the construction's choices.
"""
from functools import lru_cache
from itertools import combinations, product
from collections import Counter


class GF:
    def __init__(self, degree, modulus):
        self.e, self.q, self.modulus = degree, 1 << degree, modulus
        self.mul = tuple(tuple(self.raw_mul(a, b) for b in range(self.q))
                         for a in range(self.q))
        self.inv = (0,) + tuple(next(b for b in range(1, self.q)
                                    if self.mul[a][b] == 1)
                               for a in range(1, self.q))
        for a, b, c in product(range(self.q), repeat=3):
            assert self.mul[a][b] == self.mul[b][a]
            assert self.mul[self.mul[a][b]][c] == self.mul[a][self.mul[b][c]]
            assert self.mul[a][b ^ c] == self.mul[a][b] ^ self.mul[a][c]

    def raw_mul(self, a, b):
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a & self.q:
                a ^= self.modulus
        return out

    def sum(self, values):
        out = 0
        for x in values:
            out ^= x
        return out

    def dot(self, a, b):
        return self.sum(self.mul[x][y] for x, y in zip(a, b))

    @lru_cache(maxsize=100000)
    def det(self, a):
        if not a:
            return 1
        return self.sum(self.mul[a[0][j]][self.det(tuple(
            tuple(row[k] for k in range(len(a)) if k != j) for row in a[1:]))]
                        for j in range(len(a)) if a[0][j])

    def independent(self, vectors):
        k = len(vectors)
        if not k:
            return True
        n = len(vectors[0])
        if k > n:
            return False
        return any(self.det(tuple(tuple(v[i] for v in vectors) for i in rows))
                   for rows in combinations(range(n), k))

    @lru_cache(maxsize=100000)
    def selected(self, vectors):
        kept, indices = [], []
        for j, v in enumerate(vectors):
            if self.independent(tuple(kept) + (v,)):
                kept.append(v)
                indices.append(j)
        return tuple(indices)

    @lru_cache(maxsize=100000)
    def coordinates_data(self, basis):
        r = len(basis)
        rows = tuple(tuple(v[i] for v in basis) for i in range(len(basis[0]))) if r else ()
        indices = self.selected(rows)
        assert len(indices) == r
        d = tuple(rows[i] for i in indices)
        return indices, d, self.det(d)

    def coordinates(self, basis, v):
        indices, d, determinant = self.coordinates_data(basis)
        assert determinant
        rhs = tuple(v[i] for i in indices)
        return tuple(self.mul[self.inv[determinant]][self.det(tuple(
            tuple(rhs[i] if k == j else row[k] for k in range(len(basis)))
            for i, row in enumerate(d)))] for j in range(len(basis)))

    def lincomb(self, basis, coefficients, n):
        return tuple(self.sum(self.mul[c][v[i]] for c, v in zip(coefficients, basis))
                     for i in range(n))

    @lru_cache(maxsize=100000)
    def kernel(self, a):
        n = len(a)
        columns = transpose(a)
        indices = self.selected(columns)
        basis = tuple(columns[j] for j in indices)
        out = []
        for j in range(n):
            if j not in indices:
                alpha = self.coordinates(basis, columns[j])
                v = [int(i == j) for i in range(n)]
                for i, coefficient in zip(indices, alpha):
                    v[i] ^= coefficient
                v = tuple(v)
                assert all(self.dot(row, v) == 0 for row in a)
                out.append(v)
        assert self.independent(tuple(out))
        return tuple(out)

    @lru_cache(maxsize=100000)
    def adapted(self, image, kernel_matrix):
        columns = transpose(image)
        basis = [columns[j] for j in self.selected(columns)]
        for v in self.kernel(kernel_matrix):
            if self.independent(tuple(basis) + (v,)):
                basis.append(v)
        return tuple(basis)

    def psi(self, x, u, v, inverse=False):
        ex = self.adapted(x, x)
        eu = self.adapted(u, transpose(u))
        if inverse:
            ex, eu = eu, ex
        assert len(ex) == len(eu)
        coordinates = self.coordinates(ex, v)
        assert self.lincomb(ex, coordinates, len(v)) == v
        return self.lincomb(eu, coordinates, len(v))

    @lru_cache(maxsize=100000)
    def phi(self, x):
        if not x:
            return ()
        x0 = truncate(x)
        u = self.phi(x0)
        w = self.psi(x0, u, tuple(row[-1] for row in x[:-1]))
        return extend(u, w, self.sum(w))

    @lru_cache(maxsize=100000)
    def inverse_phi(self, u):
        if not u:
            return ()
        u0, w = truncate(u), tuple(row[-1] for row in u[:-1])
        assert u[-1][-1] == self.sum(w)
        x = self.inverse_phi(u0)
        return extend(x, self.psi(x, u0, w, inverse=True), 0)

    def rank(self, a):
        """Independent Gaussian elimination; never used to choose the map."""
        a = [list(row) for row in a]
        r = 0
        for j in range(len(a[0]) if a else 0):
            pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
            if pivot is None:
                continue
            a[r], a[pivot] = a[pivot], a[r]
            a[r] = [self.mul[self.inv[a[r][j]]][v] for v in a[r]]
            for i in range(r + 1, len(a)):
                c = a[i][j]
                a[i] = [v ^ self.mul[c][w] for v, w in zip(a[i], a[r])]
            r += 1
        return r

    def zero_product(self, a, b):
        return all(self.dot(row, col) == 0 for row in a for col in transpose(b))

    def frobenius(self, a):
        return tuple(tuple(self.mul[v][v] for v in row) for row in a)


def transpose(a):
    return tuple(zip(*a))


def truncate(a):
    return tuple(row[:-1] for row in a[:-1])


def extend(a, v, c):
    return tuple(row + (v[i],) for i, row in enumerate(a)) + ((0,) * len(a) + (c,),)


def triangular(n, q, strict):
    positions = [(i, j) for i in range(n) for j in range(i + int(strict), n)]
    for values in product(range(q), repeat=len(positions)):
        a = [[0] * n for _ in range(n)]
        for (i, j), v in zip(positions, values):
            a[i][j] = v
        yield tuple(map(tuple, a))


def recurrence(q, n):
    previous = [1]
    for size in range(1, n + 1):
        current = [1] + [0] * (size // 2)
        for r in range(1, len(current)):
            current[r] = q ** r * (previous[r] if r < len(previous) else 0)
            current[r] += (q ** (size-r) - q ** (r-1)) * previous[r-1]
        previous = current
    return previous


def check_field(f, max_n, brute_c_n):
    for n in range(max_n + 1):
        bs = {x for x in triangular(n, f.q, True) if f.zero_product(x, x)}
        counts = Counter(f.rank(x) for x in bs)
        actual = [counts[r] for r in range(n // 2 + 1)]
        assert actual == recurrence(f.q, n), (f.q, n, actual)
        image = set()
        for x in bs:
            u = f.phi(x)
            assert f.zero_product(transpose(u), u)
            assert all(u[i][j] == 0 for i in range(n) for j in range(i))
            assert f.rank(x) == f.rank(u)
            assert f.inverse_phi(u) == x
            assert f.phi(f.frobenius(x)) == f.frobenius(u)
            assert not n or truncate(u) == f.phi(truncate(x))
            image.add(u)
        assert len(image) == len(bs)
        if n <= brute_c_n:
            cs = {u for u in triangular(n, f.q, False)
                  if f.zero_product(transpose(u), u)}
            assert image == cs
            assert all(f.phi(f.inverse_phi(u)) == u for u in cs)
            coverage = 'exhaustive B and independently exhaustive C'
        else:
            assert len(image) == sum(recurrence(f.q, n))
            coverage = 'exhaustive B; C cardinality checked against source recurrence'
        print(f'GF({f.q}) n={n}: rank_counts={actual}; {coverage}; '
              'inverse, rank, Gram, Frobenius, truncation OK', flush=True)
    return bs


def check_embedding(small, large, source):
    roots = [a for a in range(large.q) if large.mul[a][a] ^ a ^ 1 == 0]
    assert len(roots) == 2
    for root in roots:
        embedding = (0, 1, root, root ^ 1)
        for a, b in product(range(4), repeat=2):
            assert embedding[a ^ b] == embedding[a] ^ embedding[b]
            assert embedding[small.mul[a][b]] == large.mul[embedding[a]][embedding[b]]
        def mapped(a):
            return tuple(tuple(embedding[v] for v in row) for row in a)
        for x in source:
            assert large.phi(mapped(x)) == mapped(small.phi(x))
        print(f'GF(4) -> GF(16), generator={root}: all {len(source)} '
              'size-4 square-zero matrices commute with embedding', flush=True)


def main():
    f2, f4, f8, f16 = GF(1, 0b11), GF(2, 0b111), GF(3, 0b1011), GF(4, 0b10011)
    print('Exact binary-polynomial field arithmetic; field laws exhaustively verified for q=2,4,8,16.')
    print('Fidelity: Proposition 6.1 and Theorem 6.2 even-characteristic rank recurrence.')
    check_field(f2, 5, 5)
    source = check_field(f4, 4, 4)
    check_field(f8, 3, 3)
    check_embedding(f4, f16, source)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
