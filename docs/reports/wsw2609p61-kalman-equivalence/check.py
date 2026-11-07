#!/usr/bin/env python3
"""Exact finite tests of the proposed Wang--Shi--Wang Kalman extension.

No floating-point values, external modules, or downloaded program source are used.
The enumeration covers all unary FARL expansions of the explicitly listed reducts,
not all finite residuated lattices. Passing is evidence for a route, not a proof.
"""
from itertools import product


def require(ok, label, *witness):
    if not ok:
        raise AssertionError((label,) + witness)


class Algebra:
    def __init__(self, elements, le, mul, residual=None):
        self.e = tuple(elements)
        self.le = le
        self.mul = mul
        self.zero, self.one = self.e[0], self.e[-1]
        self.meet = self.table(lambda a, b: self.greatest(x for x in self.e if le(x, a) and le(x, b)))
        self.join = self.table(lambda a, b: self.least(x for x in self.e if le(a, x) and le(b, x)))
        self.imp = self.table(residual or (lambda a, b: self.greatest(x for x in self.e if le(mul(a, x), b))))

    def table(self, f):
        d = {(x, y): f(x, y) for x in self.e for y in self.e}
        return lambda x, y: d[x, y]

    def greatest(self, xs):
        xs = tuple(xs)
        found = [x for x in xs if all(self.le(y, x) for y in xs)]
        require(len(found) == 1, 'greatest', xs)
        return found[0]

    def least(self, xs):
        xs = tuple(xs)
        found = [x for x in xs if all(self.le(x, y) for y in xs)]
        require(len(found) == 1, 'least', xs)
        return found[0]

    def check_rl(self):
        e, le, m, r, j, h = self.e, self.le, self.mul, self.imp, self.join, self.meet
        for x in e:
            require(le(x, x), 'order reflexivity', x)
            require(le(self.zero, x) and le(x, self.one), 'bounds', x)
            require(m(x, self.one) == x, 'unit', x)
        for x, y in product(e, repeat=2):
            require(not (le(x, y) and le(y, x)) or x == y, 'order antisymmetry', x, y)
            require(m(x, y) in e and m(x, y) == m(y, x), 'commutativity', x, y)
        for x, y, z in product(e, repeat=3):
            require(not (le(x, y) and le(y, z)) or le(x, z), 'order transitivity', x, y, z)
            require(m(m(x, y), z) == m(x, m(y, z)), 'associativity', x, y, z)
            require(le(m(x, y), z) == le(x, r(y, z)), 'residuation', x, y, z)
            require(le(z, h(x, y)) == (le(z, x) and le(z, y)), 'meet', x, y, z)
            require(le(j(x, y), z) == (le(x, z) and le(y, z)), 'join', x, y, z)

    def distributive(self):
        return all(self.meet(x, self.join(y, z)) == self.join(self.meet(x, y), self.meet(x, z))
                   for x, y, z in product(self.e, repeat=3))


def farl(a, f, g):
    e, le, r, j, h = a.e, a.le, a.imp, a.join, a.meet
    return (all(le(f(x), x) for x in e)
            and all((le(g(x), y) == le(x, f(y))
                     and f(r(g(x), y)) == r(g(x), f(y))
                     and f(j(g(x), y)) == j(g(x), f(y))
                     and g(h(x, g(y))) == h(g(x), g(y))) for x, y in product(e, repeat=2)))


def all_farl_expansions(a):
    # FARL1 and g left adjoint to f imply x <= g(x) and g(0)=0.
    # Every possible such g is visited, and its unique right adjoint is computed.
    e, le = a.e, a.le
    choices = [[a.zero] if x == a.zero else [y for y in e if le(x, y)] for x in e]
    for values in product(*choices):
        gd = dict(zip(e, values))
        g = gd.__getitem__
        if not all(le(g(x), g(y)) for x, y in product(e, repeat=2) if le(x, y)):
            continue
        if not all(g(a.join(x, y)) == a.join(g(x), g(y)) for x, y in product(e, repeat=2)):
            continue
        fd = {y: a.greatest(x for x in e if le(g(x), y)) for y in e}
        f = fd.__getitem__
        if farl(a, f, g):
            yield f, g


def twist(a, f, g):
    e = [(x, y) for x, y in product(a.e, repeat=2) if a.mul(x, y) == a.zero]
    zero, one = (a.zero, a.one), (a.one, a.zero)
    e.remove(zero)
    if one != zero:
        e.remove(one)
    e = [zero] + e + ([one] if one != zero else [])
    le = lambda x, y: a.le(x[0], y[0]) and a.le(y[1], x[1])
    mul = lambda x, y: (a.mul(x[0], y[0]), a.meet(a.imp(x[0], y[1]), a.imp(y[0], x[1])))
    imp = lambda x, y: (a.meet(a.imp(x[0], y[0]), a.imp(y[1], x[1])), a.mul(y[1], x[0]))
    b = Algebra(e, le, mul, imp)
    b.center = (a.zero, a.zero)
    b.neg = lambda x: (x[1], x[0])
    b.f = lambda x: (f(x[0]), g(x[1]))
    b.g = lambda x: b.neg(b.f(b.neg(x)))
    return b


def target_axioms(b):
    b.check_rl()
    e, le, h, j, m, r, n, f, c = b.e, b.le, b.meet, b.join, b.mul, b.imp, b.neg, b.f, b.center
    require(n(c) == c and f(c) == c, 'center')
    signatures = {}
    for x in e:
        require(n(n(x)) == x and f(x) in e and le(f(x), x), 'involution/Fc2', x)
        require(f(h(x, c)) == h(f(x), c), 'Fc5', x)
        require(f(j(x, c)) == j(f(x), c), 'Fc6', x)
        sig = (j(x, c), j(n(x), c))
        require(sig not in signatures, 'S', x, signatures.get(sig))
        signatures[sig] = x
    for x, y in product(e, repeat=2):
        require(le(n(x), n(y)) == le(y, x), 'order reversal', x, y)
        require(r(x, y) == n(m(x, n(y))), 'involutive residual', x, y)
        require(h(m(x, y), c) == j(m(h(x, c), y), m(x, h(y, c))), 'Leibniz', x, y)
        require(f(r(f(x), y)) == r(f(x), f(y)), 'Fc3', x, y)
        require(f(j(f(x), y)) == j(f(x), f(y)), 'Fc4', x, y)
        require(j(h(x, y), c) == h(j(x, c), j(y, c)), 'D', x, y)
        if le(c, x) and le(c, y) and le(m(x, y), c):
            require((x, y) in signatures, 'CK', x, y)


def cone(b):
    e = [x for x in b.e if b.le(b.center, x)]
    a = Algebra(e, b.le, lambda x, y: b.join(b.mul(x, y), b.center), b.imp)
    return a, b.f, b.g


def check_iso(a, b, mapping, unary_a, unary_b, label):
    require(set(mapping) == set(a.e) and set(mapping.values()) == set(b.e)
            and len(set(mapping.values())) == len(a.e), label + '/bijection')
    require(mapping[a.zero] == b.zero and mapping[a.one] == b.one, label + '/constants')
    for op_a, op_b in [(a.meet, b.meet), (a.join, b.join), (a.mul, b.mul), (a.imp, b.imp)]:
        for x, y in product(a.e, repeat=2):
            require(mapping[op_a(x, y)] == op_b(mapping[x], mapping[y]), label + '/binary', x, y)
    for u, v in zip(unary_a, unary_b):
        require(all(mapping[u(x)] == v(mapping[x]) for x in a.e), label + '/unary')


def check_extension(a, f, g):
    require(farl(a, f, g), 'input FARL')
    b = twist(a, f, g)
    target_axioms(b)
    ca, cf, cg = cone(b)
    ca.check_rl()
    require(farl(ca, cf, cg), 'cone FARL')
    rho = {x: (x, a.zero) for x in a.e}
    check_iso(a, ca, rho, (f, g), (cf, cg), 'rho')
    kc = twist(ca, cf, cg)
    sigma = {x: (b.join(x, b.center), b.join(b.neg(x), b.center)) for x in b.e}
    require(sigma[b.center] == kc.center, 'sigma/center')
    check_iso(b, kc, sigma, (b.f, b.neg), (kc.f, kc.neg), 'sigma')
    # Test the multiplication step identified in the candidate's argument.
    for x, y in product(b.e, repeat=2):
        lhs = b.join(b.mul(b.join(x, b.center), b.meet(y, b.center)),
                     b.mul(b.join(y, b.center), b.meet(x, b.center)))
        require(lhs == b.meet(b.mul(x, y), b.center), 'multiplication bridge', x, y)
    return len(b.e)


def chain_reducts(n):
    e = tuple(range(n))
    slots = [(x, y) for x in e[1:-1] for y in e[1:-1] if x <= y]
    for vals in product(*(range(x + 1) for x, y in slots)):
        tab = {(x, y): v for (x, y), v in zip(slots, vals)}
        def mul(x, y):
            if not x or not y:
                return 0
            if x == n - 1:
                return y
            if y == n - 1:
                return x
            return tab[min(x, y), max(x, y)]
        if not all(mul(mul(x, y), z) == mul(x, mul(y, z)) for x, y, z in product(e, repeat=3)):
            continue
        if not all(mul(x, min(y, z)) <= mul(x, max(y, z)) for x, y, z in product(e, repeat=3)):
            continue
        # Freeze the table, so the yielded algebra cannot depend on the next iteration.
        table = {(x, y): mul(x, y) for x, y in product(e, repeat=2)}
        yield Algebra(e, lambda x, y: x <= y, lambda x, y, t=table: t[x, y])


def nondistributive_six():
    # 0 < p,q,r < e < 1; p,q,r are incomparable. Products of two nonunits are 0.
    le = lambda x, y: x == 0 or y == 5 or x == y or (x in (1, 2, 3) and y == 4)
    return Algebra(range(6), le, lambda x, y: y if x == 5 else x if y == 5 else 0)


def main():
    total = 0
    for n in range(1, 6):
        reducts = expansions = 0
        for a in chain_reducts(n):
            a.check_rl()
            reducts += 1
            for f, g in all_farl_expansions(a):
                check_extension(a, f, g)
                expansions += 1
        total += expansions
        print(f'chains n={n}: {reducts} reducts, {expansions} FARL expansions checked')
    boolean = Algebra(range(4), lambda x, y: x & y == x, lambda x, y: x & y)
    boolean.check_rl()
    be = 0
    for f, g in all_farl_expansions(boolean):
        check_extension(boolean, f, g)
        be += 1
    total += be
    print(f'Boolean square: {be} FARL expansions checked')
    nd = nondistributive_six()
    nd.check_rl()
    require(not nd.distributive(), 'non-distributive input')
    require(nd.meet(1, nd.join(2, 3)) == 1 and nd.join(nd.meet(1, 2), nd.meet(1, 3)) == 0,
            'explicit distributivity failure')
    ne = 0
    for f, g in all_farl_expansions(nd):
        require(check_extension(nd, f, g) == 27, '27-element twist')
        ne += 1
    total += ne
    print(f'Non-distributive six-element source: {ne} FARL expansions; all twists have 27 elements')
    # Source fidelity: the displayed L3 x L3 example with diagonal meet/max operators.
    e = tuple(product(range(3), repeat=2))
    le = lambda x, y: all(a <= b for a, b in zip(x, y))
    mul = lambda x, y: tuple(max(0, a + b - 2) for a, b in zip(x, y))
    a = Algebra(e, le, mul)
    a.check_rl()
    f = lambda x: (min(x), min(x))
    g = lambda x: (max(x), max(x))
    neg = lambda x: a.imp(x, a.zero)
    for x in a.e:
        require(neg(neg(x)) == x, 'source FAIRL involutivity', x)
        require(a.le(f(x), x), 'source FAIRL1', x)
        require(g(x) == neg(f(neg(x))), 'source derived g', x)
    for x, y in product(a.e, repeat=2):
        require(a.imp(x, y) == tuple(min(2, 2 - xx + yy) for xx, yy in zip(x, y)),
                'source displayed residual', x, y)
        require(f(a.join(f(x), y)) == a.join(f(x), f(y)), 'source FARL14', x, y)
        require(f(a.imp(f(x), y)) == a.imp(f(x), f(y)), 'source FARL17', x, y)
    require(farl(a, f, g), 'source L3-square example')
    size = check_extension(a, f, g)
    print(f'Source fidelity: L3 x L3, scaled integer coordinates 0,1,2: FAIRL and FARL; twist size={size}; rho and sigma preserve every operation')
    print(f'Total exhaustive family expansions: {total}; plus the source L3-square example')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
