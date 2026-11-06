#!/usr/bin/env python3
"""Exact finite checks for Paulista, arXiv:2511.10612v1, Problem 7.1.

No external packages, randomness, floating point, or Lean proofs are used.
A proper coloring and a clique of equal size certify both graph invariants.
"""
from itertools import permutations, product


class Field:
    """F_p, or F_p[t]/(t^2-d) for a nonsquare d; encoded by integers."""

    def __init__(self, p, degree=1):
        assert p in (3, 5, 7, 11) and degree in (1, 2)
        self.p, self.degree, self.q = p, degree, p ** degree
        d = next(x for x in range(1, p) if pow(x, (p - 1) // 2, p) == p - 1)
        self.d = d
        q = self.q
        self.add = [[0] * q for _ in range(q)]
        self.mul = [[0] * q for _ in range(q)]
        for x, y in product(range(q), repeat=2):
            a, b = x % p, x // p
            c, e = y % p, y // p
            self.add[x][y] = (a + c) % p + p * ((b + e) % p)
            self.mul[x][y] = (a * c + d * b * e) % p + p * ((a * e + b * c) % p)
        self.neg = [next(y for y in range(q) if self.add[x][y] == 0) for x in range(q)]
        self.inv = [None] + [next(y for y in range(1, q) if self.mul[x][y] == 1)
                             for x in range(1, q)]
        self.squares = {self.mul[x][x] for x in range(1, q)}
        # Verify the extension really is a field, including associativity/distributivity.
        for x, y, z in product(range(q), repeat=3):
            assert self.add[self.add[x][y]][z] == self.add[x][self.add[y][z]]
            assert self.mul[self.mul[x][y]][z] == self.mul[x][self.mul[y][z]]
            assert self.mul[x][self.add[y][z]] == self.add[self.mul[x][y]][self.mul[x][z]]

    def sub(self, x, y):
        return self.add[x][self.neg[y]]

    def char(self, x):
        return 0 if x == 0 else (1 if x in self.squares else -1)


def witnesses(adj, colors, clique, expected):
    """Check an actual clique and a proper coloring of matching cardinality."""
    assert len(colors) == len(adj)
    assert len(clique) == expected == len(set(clique))
    assert len(set(colors)) == expected
    assert all(0 <= c < expected for c in colors)
    assert all(colors[v] != colors[w] for v in range(len(adj)) for w in adj[v])
    assert all(w in adj[v] for v in clique for w in clique if v != w)


def product_check(multiply, elements, adj, colors, clique, m, expected):
    """Construct G x C_m and independently recompute every commutation pair.

    Removing the *whole* center is checked before comparing with the blow-up.
    """
    elems = list(product(range(len(elements)), range(m)))

    def mul(x, y):
        return multiply(elements[x[0]], elements[y[0]]), (x[1] + y[1]) % m

    central = {x for x in elems if all(mul(x, y) == mul(y, x) for y in elems)}
    base_central = {i for i, x in enumerate(elements)
                    if all(multiply(x, y) == multiply(y, x) for y in elements)}
    assert central == set(product(base_central, range(m)))
    vertices = [x for x in elems if x not in central]
    base_vertices = [i for i in range(len(elements)) if i not in base_central]
    index = {v: i for i, v in enumerate(base_vertices)}
    lifted_colors = {(v, a): colors[index[v]] * m + a for v, a in vertices}
    for i, x in enumerate(vertices):
        for y in vertices[i + 1:]:
            commutes = mul(x, y) == mul(y, x)
            u, v = index[x[0]], index[y[0]]
            assert commutes == (u == v or v in adj[u])
            if commutes:
                assert lifted_colors[x] != lifted_colors[y]
    lifted_clique = [(base_vertices[v], a) for v in clique for a in range(m)]
    assert len(lifted_clique) == expected * m
    assert all(mul(x, y) == mul(y, x) for x in lifted_clique for y in lifted_clique)
    assert len(set(lifted_colors.values())) == expected * m
    return len(central), len(vertices)


def fidelity():
    """Reproduce source Theorems 4.3 and 5.3: S_3 x C_m has omega=chi=2m.

    The theorem numbers are checked against the rendered arXiv v1 source.
    """
    elements = list(permutations(range(3)))
    mul = lambda x, y: tuple(x[y[i]] for i in range(3))
    center = [x for x in elements if all(mul(x, y) == mul(y, x) for y in elements)]
    assert center == [(0, 1, 2)]
    vertices = [x for x in elements if x not in center]
    adj = [{j for j, y in enumerate(vertices) if i != j and mul(x, y) == mul(y, x)}
           for i, x in enumerate(vertices)]
    edge = next((i, j) for i in range(len(adj)) for j in adj[i])
    colors = [0] * len(adj)
    colors[edge[1]] = 1
    witnesses(adj, colors, list(edge), 2)
    for m in range(1, 7):
        z, v = product_check(mul, elements, adj, colors, list(edge), m, 2)
        assert z == m and v == 5 * m
        print(f'FIDELITY S3xC{m}: order={6*m} center={z} vertices={v} omega=chi={2*m}')


def pgl_check(p, degree):
    f = Field(p, degree)
    q, add, mul = f.q, f.add, f.mul
    identity = (1, 0, 0, 1)

    def det(a):
        return f.sub(mul[a[0]][a[3]], mul[a[1]][a[2]])

    def normalize(a):
        scale = f.inv[next(x for x in a if x)]
        return tuple(mul[x][scale] for x in a)

    def raw_product(a, b):
        return (add[mul[a[0]][b[0]]][mul[a[1]][b[2]]],
                add[mul[a[0]][b[1]]][mul[a[1]][b[3]]],
                add[mul[a[2]][b[0]]][mul[a[3]][b[2]]],
                add[mul[a[2]][b[1]]][mul[a[3]][b[3]]])

    def pmul(a, b):
        return normalize(raw_product(a, b))

    elements = sorted({normalize(a) for a in product(range(q), repeat=4) if det(a)})
    assert len(elements) == q * (q * q - 1)
    vertices = [a for a in elements if a != identity]
    adj = [set() for _ in vertices]
    for i, a in enumerate(vertices):
        for j in range(i + 1, len(vertices)):
            if pmul(a, vertices[j]) == pmul(vertices[j], a):
                adj[i].add(j)
                adj[j].add(i)
    assert all(len(n) < len(vertices) - 1 for n in adj)  # Z(PGL_2)=1
    stage, centralizers = [], {0: set(), 1: set(), 2: set()}
    four = 4 % p
    for i, a in enumerate(vertices):
        tr = add[a[0]][a[3]]
        involution = pmul(a, a) == identity
        assert involution == (tr == 0)
        if involution:
            d = raw_product(a, a)[0]
            kind = 1 if f.char(d) == 1 else 2
            expected = 2 * (q - 1 if kind == 1 else q + 1)
        else:
            kind = 0
            discr = f.sub(mul[tr][tr], mul[four][det(a)])
            expected = q if discr == 0 else (q - 1 if f.char(discr) == 1 else q + 1)
        assert len(adj[i]) + 2 == expected
        stage.append(kind)
        centralizers[kind].add(expected)
    for i, kind in enumerate(stage):
        if kind == 1:
            assert sum(stage[j] != 0 for j in adj[i]) == q - 1
        if kind == 2:
            assert sum(stage[j] == 2 for j in adj[i]) == (q + 1) // 2
    # Verify the precise deletion order, not just existence of a q-coloring.
    order = sorted(range(len(vertices)), key=lambda i: (stage[i], i))
    remaining = set(order)
    maxima = [0, 0, 0]
    for i in order:
        later = len(adj[i] & remaining)
        maxima[stage[i]] = max(maxima[stage[i]], later)
        assert later <= q - 1
        remaining.remove(i)
    colors = [-1] * len(vertices)
    for i in reversed(order):
        used = {colors[j] for j in adj[i] if colors[j] >= 0}
        colors[i] = next(c for c in range(q) if c not in used)
    d = next(x for x in range(1, q) if f.char(x) == -1)
    torus = {normalize((a, mul[b][d], b, a))
             for a, b in product(range(q), repeat=2) if a or b}
    assert len(torus) == q + 1
    idx = {a: i for i, a in enumerate(vertices)}
    clique = [idx[a] for a in torus if a != identity]
    witnesses(adj, colors, clique, q)
    charsum = sum(f.char(f.sub(mul[a][a], d)) for a in range(q))
    assert charsum == -1
    projective_counts = [f.char(f.sub(mul[a][a], d)) for a in range(q)] + [1]
    assert projective_counts.count(-1) == (q + 1) // 2
    print(f'PGL2 F{q}: order={len(elements)} center=1 vertices={len(vertices)} '
          f'involutions(split,nonsplit)=({stage.count(1)},{stage.count(2)}) '
          f'deletion_max={tuple(maxima)} character_sum={charsum} omega=chi={q}')
    if q <= 7:
        for m in (3, 5):
            z, v = product_check(pmul, elements, adj, colors, clique, m, q)
            assert z == m and v == (len(elements) - 1) * m
            print(f'PRODUCT PGL2(F{q})xC{m}: center={z} vertices={v} omega=chi={q*m}')


def main():
    fidelity()
    for p, degree in ((3, 1), (5, 1), (7, 1), (3, 2), (11, 1)):
        pgl_check(p, degree)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
