#!/usr/bin/env python3
"""Exact D-type weight multiplicities and atomic coefficients.

Independent implementations: Freudenthal recurrence with Weyl normalization,
Kostant's signed Weyl sum, and direct dominant-interval Mobius inversion.
All coordinates use fractions.Fraction; there is no floating point arithmetic.
Mathematical definitions: Hindi--Schutzer, arXiv:2609.08024v1, Sections 2, 3, 6.
"""
from fractions import Fraction as Q
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations, product
import json


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def dot(x, y):
    return sum((a * b for a, b in zip(x, y)), Q(0))


def scale(k, x):
    return tuple(k * a for a in x)


class D:
    def __init__(self, n):
        assert n >= 4
        self.n = n
        self.zero = (Q(0),) * n
        e = [tuple(Q(i == j) for j in range(n)) for i in range(n)]
        self.simple = tuple(sub(e[i], e[i + 1]) for i in range(n - 1)) + (add(e[-2], e[-1]),)
        self.positive = tuple(sub(e[i], e[j]) for i in range(n) for j in range(i + 1, n)) + tuple(add(e[i], e[j]) for i in range(n) for j in range(i + 1, n))
        self.rho = tuple(Q(n - 1 - i) for i in range(n))
        self.root_vectors = tuple(self.rootcoords(a) for a in self.positive)
        self.complement_vectors = tuple(self.rootcoords(a) for a in self.positive if a[0] != 0 or sum(a) == 2)
        assert len(set(self.root_vectors)) == n * (n - 1)
        assert all(all(v >= 0 and v.denominator == 1 for v in a) for a in self.root_vectors)

    def weight(self, labels):
        assert len(labels) == self.n
        n = self.n
        return tuple(sum((Q(labels[i]) for i in range(j, n - 2)), Q(0)) + Q(labels[-2] + labels[-1], 2) for j in range(n - 1)) + (Q(labels[-1] - labels[-2], 2),)

    def labels(self, x):
        return tuple(dot(x, a) for a in self.simple)

    def rootcoords(self, x):
        c = tuple(sum(x[:i + 1], Q(0)) for i in range(self.n - 2))
        total = sum(x[:-1], Q(0))
        return c + ((total - x[-1]) / 2, (total + x[-1]) / 2)

    def fromroots(self, c):
        return tuple(sum((Q(c[i]) * self.simple[i][j] for i in range(self.n)), Q(0)) for j in range(self.n))

    def below(self, x, y):
        return all(v >= 0 and v.denominator == 1 for v in self.rootcoords(sub(y, x)))

    def dominant(self, x):
        return all(v >= 0 for v in self.labels(x))

    def normalize(self, x):
        y = sorted((abs(v) for v in x), reverse=True)
        if sum(v < 0 for v in x) % 2:
            y[-1] = -y[-1]
        y = tuple(y)
        assert self.dominant(y)
        return y

    def multiplicities(self, lam):
        assert self.dominant(lam)
        bound = dot(lam, lam)

        @cache
        def m(nu):
            nu = self.normalize(nu)
            if not self.below(nu, lam):
                return 0
            if nu == lam:
                return 1
            den = dot(add(lam, self.rho), add(lam, self.rho)) - dot(add(nu, self.rho), add(nu, self.rho))
            assert den > 0
            num = Q(0)
            for a in self.positive:
                v = add(nu, a)
                while dot(v, v) <= bound:
                    num += 2 * dot(v, a) * m(v)
                    v = add(v, a)
            out = num / den
            assert out.denominator == 1 and out >= 0, (lam, nu, out)
            return int(out)
        return m

    def interval(self, mu, lam):
        d = self.rootcoords(sub(lam, mu))
        assert all(v >= 0 and v.denominator == 1 for v in d)
        points = [add(mu, self.fromroots(c)) for c in product(*(range(int(v) + 1) for v in d))]
        points = [v for v in points if self.dominant(v)]
        return sorted(points, key=lambda v: sum(self.rootcoords(sub(v, mu))))

    def atomic(self, mu, lam, m=None):
        if not self.below(mu, lam):
            return 0, [], [], {}
        if m is None:
            m = self.multiplicities(lam)
        interval = self.interval(mu, lam)
        mobius = {}
        for eta in interval:
            mobius[eta] = 1 if eta == mu else -sum(mobius[z] for z in mobius if self.below(z, eta))
        value = sum(mobius[z] * m(z) for z in interval)
        # Independent triangular inversion of the defining character expansion.
        atomic = {}
        for eta in reversed(interval):
            atomic[eta] = m(eta) - sum(atomic[z] for z in atomic if self.below(eta, z))
        assert atomic[mu] == value
        atoms = [z for z in interval if z != mu and not any(t != mu and t != z and self.below(t, z) for t in interval)]
        # Verify the crosscut formula and every join, without assuming a three-atom shape.
        crosscut = 0
        joins = []
        for mask in range(1 << len(atoms)):
            selected = [a for i, a in enumerate(atoms) if mask >> i & 1]
            upper = [z for z in interval if all(self.below(a, z) for a in selected)]
            minima = [z for z in upper if not any(t != z and self.below(t, z) for t in upper)]
            assert len(minima) == 1
            join = minima[0]
            crosscut += (-1) ** len(selected) * m(join)
            joins.append((mask, tuple(int(v) for v in self.rootcoords(sub(join, mu))), m(join)))
        assert crosscut == value
        return value, atoms, joins, mobius

    def graph_partitions(self, beta):
        """Remaining-root partitions as loopless multigraphs after A cancellation."""
        total = sum(beta, Q(0)) / 2
        if total < 0 or total.denominator != 1:
            return 0
        assert total <= 2
        edges = tuple(combinations(range(self.n), 2))
        count = 0
        for chosen in combinations_with_replacement(edges, int(total)):
            deg = [0] * self.n
            for i, j in chosen:
                deg[i] += 1
                deg[j] += 1
            if all(deg[i] >= beta[i] for i in range(1, self.n)):
                count += 1
        return count

    def weyl(self):
        for p in permutations(range(self.n)):
            parity = sum(p[i] > p[j] for i in range(self.n) for j in range(i + 1, self.n)) % 2
            for signs in product((-1, 1), repeat=self.n - 1):
                ss = signs + ((-1) ** sum(s == -1 for s in signs),)
                yield p, ss, (-1) ** parity

    @cache
    def partitions(self, beta, index=0, complement=False):
        if any(v < 0 or Q(v).denominator != 1 for v in beta):
            return 0
        if all(v == 0 for v in beta):
            return 1
        roots = self.complement_vectors if complement else self.root_vectors
        if index == len(roots):
            return 0
        a = roots[index]
        upper = min(int(v // r) for v, r in zip(beta, a) if r)
        return sum(self.partitions(tuple(v - k * r for v, r in zip(beta, a)), index + 1, complement) for k in range(upper + 1))

    def kostant_table(self, lam, targets):
        out = {eta: 0 for eta in targets}
        support = {eta: 0 for eta in targets}
        v = add(lam, self.rho)
        for p, ss, sign in self.weyl():
            wv = tuple(ss[i] * v[p[i]] for i in range(self.n))
            for eta in targets:
                beta = self.rootcoords(sub(sub(wv, self.rho), eta))
                if all(b >= 0 and b.denominator == 1 for b in beta):
                    # Parabolic {2,...,n-1} permutes coordinates 2,...,n.
                    assert p[0] == 0 and all(s == 1 for s in ss), (self.n, lam, eta, p, ss, beta)
                    support[eta] += 1
                    out[eta] += sign * self.partitions(beta)
        return out, support


def as_labels(r, x):
    return tuple(int(v) for v in r.labels(x))


def main():
    # Fidelity: source Proposition 6.8, its two displayed D4 pairs.
    d4 = D(4)
    for ml, ll, expected in [((0, 1, 0, 0), (0, 0, 2, 2), -2), ((1, 1, 0, 0), (1, 0, 2, 2), -1)]:
        value, atoms, joins, mobius = d4.atomic(d4.weight(ml), d4.weight(ll))
        assert value == expected
        print('FIDELITY', json.dumps(dict(type='D4', mu=ml, lam=ll, atomic=value, interval_size=len(mobius))))
    for n in range(4, 9):
        r = D(n)
        for s in range(4):
            ml = (0,) * (n - 3) + (1, 0, s)
            ll = (2,) + (0,) * (n - 2) + (s + 2,)
            pl = (0,) * (n - 1) + (s + 2,)
            mu, lam, prime = map(r.weight, (ml, ll, pl))
            val, atoms, joins, mobius = r.atomic(mu, lam)
            assert 2 * val == n * n - 9 * n + 16, (n, s, val)
            assert len(atoms) == 3
            mult = r.multiplicities(lam)
            assert all(mult(z) == r.partitions(r.rootcoords(sub(lam, z)), complement=True) == r.graph_partitions(sub(lam, z)) for z in mobius)
            expected_multiplicities = {(0,) * n: (3*n*n - 9*n + 8)//2,
                (1,) * (n - 3) + (0, 0, 0): 3*(n - 2),
                (0,) * (n - 2) + (1, 0): n*(n - 1)//2,
                (0,) * (n - 1) + (1,): n*(n - 1)//2,
                (1,) * n: n - 1}
            if s >= 1:
                expected_multiplicities[(1,) * (n - 1) + (0,)] = n - 1
            for c, expected_m in expected_multiplicities.items():
                assert mult(add(mu, r.fromroots(c))) == expected_m
            expected_increments = {(1,) * (n - 3) + (0, 0, 0), (0,) * (n - 2) + (1, 0), (0,) * (n - 1) + (1,)}
            assert {tuple(int(v) for v in r.rootcoords(sub(a, mu))) for a in atoms} == expected_increments
            sharper = r.atomic(mu, prime)[0]
            assert sharper == 0 and sub(prime, mu) == r.simple[-1]
            indices = (n - 4, n - 3, n - 2, n - 1)
            local_mu = tuple(ml[i] for i in indices)
            local_lam = tuple(ll[i] for i in indices)
            assert local_lam == tuple(pl[i] for i in indices) or n == 4
            lv, la, lj, lm = d4.atomic(d4.weight(local_mu), d4.weight(local_lam))
            if n >= 5:
                assert local_mu == (0, 1, 0, s) and local_lam == (0, 0, 0, s + 2)
                assert lv == 0 and len(lm) == 2 and all(d4.multiplicities(d4.weight(local_lam))(v) == 1 for v in lm)
            if n >= 5:
                for outer in permutations((n - 4, n - 2, n - 1)):
                    star = (outer[0], n - 3, outer[1], outer[2])
                    a = d4.weight(tuple(ml[i] for i in star))
                    b = d4.weight(tuple(ll[i] for i in star))
                    assert d4.atomic(a, b)[0] == 0
            detail = dict(rank=n, s=s, atomic=val, local_atomic=lv, prime_atomic=sharper, interval_size=len(mobius), atoms=[as_labels(r,a) for a in atoms], joins=joins)
            if n <= 6 and s <= 1:
                m = r.multiplicities(lam)
                kt, support = r.kostant_table(lam, list(mobius))
                assert all(kt[z] == m(z) for z in mobius)
                detail['kostant_checked_weights'] = len(kt)
                detail['nonzero_mobius_rows'] = [dict(labels=as_labels(r,z), mobius=mobius[z], multiplicity=m(z)) for z in mobius if mobius[z]]
                detail['parabolic_support_terms_max'] = max(support.values())
            print('FAMILY', json.dumps(detail, sort_keys=True))
        r.partitions.cache_clear()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
