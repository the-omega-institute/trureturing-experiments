#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exact finite checks for the actual biased trace and occurrence covariance.

Arithmetic in Q(a), a^2=1-a, uses rational coefficient pairs. Mass stopping
uses integer address weights, including equality, never floating comparisons.
The checks do not prove an infinite spectral asymptotic.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import json
import numpy as np


@dataclass(frozen=True)
class Q:
    x: F = F(0)
    y: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "x", F(self.x))
        object.__setattr__(self, "y", F(self.y))

    @staticmethod
    def coerce(v):
        return v if isinstance(v, Q) else Q(v)

    def __add__(self, v):
        v = Q.coerce(v)
        return Q(self.x + v.x, self.y + v.y)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.x, -self.y)

    def __sub__(self, v):
        return self + -Q.coerce(v)

    def __mul__(self, v):
        v = Q.coerce(v)
        return Q(self.x*v.x + self.y*v.y,
                 self.x*v.y + self.y*v.x - self.y*v.y)

    __rmul__ = __mul__

    def __truediv__(self, v):
        v = Q.coerce(v)
        norm = v.x*v.x - v.x*v.y - v.y*v.y
        if norm == 0:
            raise ZeroDivisionError
        return self * Q((v.x-v.y)/norm, -v.y/norm)

    def sign(self):
        # 2*(x+y*a) = (2*x-y) + y*sqrt(5); compare exact squares.
        u, v = 2*self.x-self.y, self.y
        if v == 0:
            return (u > 0) - (u < 0)
        if u == 0:
            return (v > 0) - (v < 0)
        if (u > 0) == (v > 0):
            return (u > 0) - (u < 0)
        cmp = u*u - 5*v*v
        return ((cmp > 0) - (cmp < 0)) * ((u > 0) - (u < 0))

    def serial(self):
        return [str(self.x), str(self.y)]


ZERO, ONE, A = Q(), Q(1), Q(0, 1)
B = ONE-A
assert A*A == B and A+B == ONE
assert A.sign() > 0 and B.sign() > 0 and (A-B).sign() > 0


@lru_cache(None)
def power(n):
    return ONE if n == 0 else power(n-1)*A


def weight(p):
    return p.count("L") + 2*p.count("R")


def mass(p):
    return power(weight(p))


def words(depth):
    return ["".join(p) for j in range(1, depth+1)
            for p in product("LR", repeat=j)]


def cut(k):
    todo, leaves = [""], []
    while todo:
        p = todo.pop()
        if weight(p) >= k:  # mass(p) <= a^k, equality is a stop
            leaves.append(p)
        else:
            todo.extend([p+"R", p+"L"])
    return sorted(leaves)


def meet_mass(u, v):
    if u.startswith(v):
        return mass(u)
    if v.startswith(u):
        return mass(v)
    return ZERO


def covariance(u, v):
    return meet_mass(u, v) - mass(u)*mass(v)


def cut_checks():
    rows = []
    for k in range(1, 17):
        leaves, t = cut(k), power(k)
        assert sum((mass(p) for p in leaves), ZERO) == ONE
        assert all(weight(p) >= k > weight(p[:-1]) for p in leaves)
        assert all((mass(p)-B*t).sign() > 0 and (t-mass(p)).sign() >= 0
                   for p in leaves)
        assert (len(leaves)*t-ONE).sign() >= 0
        assert (ONE-len(leaves)*B*t).sign() > 0
        equals = sum(weight(p) == k for p in leaves)
        assert equals > 0
        rows.append(dict(k=k, threshold_q_a=t.serial(), cells=len(leaves),
                         exact_threshold_equalities=equals,
                         total_mass_q_a=ONE.serial()))
    return rows


def residual_gram_checks():
    # Test all matrix entries, not a selected vector's accidental cancellation.
    edges = words(5)
    entries, nonzero = 0, 0
    for k in (1, 2, 3, 5):
        leaves = cut(k)
        # Conditional Gram of I-Q directly from cylinder intersections.
        for u in edges:
            for v in edges:
                left = meet_mass(u, v) - sum(
                    (meet_mass(u, p)*meet_mass(v, p)/mass(p) for p in leaves), ZERO)
                # Independent suffix covariance formula for strict descendants.
                right = ZERO
                for p in leaves:
                    if (len(u) > len(p) and len(v) > len(p)
                            and u.startswith(p) and v.startswith(p)):
                        right += mass(p)*covariance(u[len(p):], v[len(p):])
                assert left == right, (k, u, v, left, right)
                nonzero += left != ZERO
                entries += 1
    # Replacing the conditional biased suffix covariance with the fair law
    # changes the variance of the L child on every stopped cylinder.
    fair_variance = Q(F(1, 4))
    assert covariance("L", "L") != fair_variance
    return dict(cuts=[1, 2, 3, 5], edge_support_max_depth=5,
                matrix_entries=entries, nonzero_entries=nonzero,
                biased_child_variance_q_a=covariance("L", "L").serial(),
                fair_child_variance_q_a=fair_variance.serial(),
                fair_substitution_detected=True)


@lru_cache(None)
def nodes(n):
    if n < 2:
        return ("",)
    return ("",) + tuple("L"+p for p in nodes(n-1)) + tuple("R"+p for p in nodes(n-2))


def source_checks():
    rows = []
    candidates = words(7)
    for n in range(13):
        vs = nodes(n)
        for p in candidates:
            reached = n-weight(p) >= 0 and n-weight(p[:-1]) >= 2
            actual = sum(x.startswith(p) for x in vs)
            expected = len(nodes(n-weight(p))) if reached else 0
            assert actual == expected, (n, p, actual, expected)
        rows.append(dict(n=n, nodes=len(vs), tested_addresses=len(candidates)))
    return rows


def finite_inverse_checks():
    rows = []
    for n in range(9):
        vs = nodes(n)
        es = vs[1:]
        size = len(vs)
        lookup = {p: j for j, p in enumerate(vs)}
        # J is path summation, D is the literal unit-edge difference.
        J = np.array([[int(x.startswith(p)) for p in es] for x in vs], dtype=object)
        D = np.zeros((size-1, size), dtype=object)
        for i, p in enumerate(es):
            D[i, lookup[p]] = 1
            D[i, lookup[p[:-1]]] = -1
        L = D.T @ D
        I = np.eye(size, dtype=object)
        S = size*I - np.ones((size, size), dtype=object)  # size*Pi
        SJ = S @ J
        assert np.array_equal(D @ J, np.eye(size-1, dtype=object))
        assert np.array_equal(SJ @ SJ.T @ L, size*S)
        counts = J.sum(axis=0)
        sigma_num = size*(J.T @ J) - np.outer(counts, counts)
        assert np.array_equal(SJ.T @ SJ, size*sigma_num)
        assert np.array_equal(L @ np.ones(size, dtype=object), np.zeros(size, dtype=object))
        # These integer identities are equivalent to
        # C_e C_e^T (size L) = Pi and C_e^T C_e = Sigma_n,
        # with C_e = Pi J/sqrt(size) the normalized Euclidean representation.
        rows.append(dict(n=n, nodes=size, edges=size-1,
                         normalized_inverse_identity=True,
                         centered_covariance_identity=True,
                         constant_zero_mode=True))
    return rows


def main():
    result = dict(arithmetic="Q(a), a^2+a-1=0; exact rational pairs; integer matrices",
                  threshold_rule="stop iff address weight >= k; includes exact equality",
                  mass_cuts=cut_checks(), residual_gram=residual_gram_checks(),
                  finite_sources=source_checks(), finite_inverse=finite_inverse_checks(),
                  limitation="Finite identities only; not proof of infinite spectral growth.")
    output = Path(__file__).with_name("results.json")
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps({"mass_cuts":len(result["mass_cuts"]),
                      "residual_gram":result["residual_gram"],
                      "finite_sources":len(result["finite_sources"]),
                      "finite_inverse":len(result["finite_inverse"]),
                      "result":str(output)},ensure_ascii=False))
    print("ALL_EXACT_FINITE_IDENTITIES_PASSED")


if __name__ == "__main__":
    main()
