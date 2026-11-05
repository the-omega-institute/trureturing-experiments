#!/usr/bin/env python3
"""Exact finite checks of Hanna's self-composition divisibility conjecture.

Sources: OEIS A212392/A212425 COMMENTS; Abbondandolo--Majer,
arXiv:2610.33287v1, Lemma 7.1 and Remark 7.2. No external packages.
Sparse dictionaries represent ordinary polynomials modulo x^(cap+1).
All arithmetic uses Python integers or fractions.Fraction.
"""
from fractions import Fraction
import json


def add(a, b, cap):
    out = dict(a)
    for k, c in b.items():
        if k <= cap:
            out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c and k <= cap}


def scale(a, c):
    return {k: c * v for k, v in a.items() if c * v}


def mul(a, b, cap):
    out = {}
    for i, ci in a.items():
        for j, cj in b.items():
            if i + j <= cap:
                out[i + j] = out.get(i + j, 0) + ci * cj
    return {k: c for k, c in out.items() if c}


def power(a, k, cap):
    out = {0: 1}
    while k:
        if k % 2:
            out = mul(out, a, cap)
        k //= 2
        if k:
            a = mul(a, a, cap)
    return out


def compose(outer, inner, cap):
    assert inner.get(0, 0) == 0
    order = min(inner) if inner else cap + 1
    out = {}
    for k, c in outer.items():
        if k * order <= cap:
            out = add(out, scale(power(inner, k, cap), c), cap)
    return out


def fixed_point(step, start, max_steps):
    current = start
    for iteration in range(1, max_steps + 1):
        nxt = step(current)
        if nxt == current:
            return current, iteration
        current = nxt
    raise AssertionError("finite truncation did not stabilize")


def run_case(b, count):
    q = b * b - 1
    top = b + q * (count - 1)
    cap = top + 1
    x = {1: 1}
    step = lambda f: power(add(x, compose(f, f, cap), cap), b, cap)
    f, iterations = fixed_point(step, {}, count + 3)
    assert f == step(f)
    assert all(k >= b and (k - b) % q == 0 for k in f)
    assert all(isinstance(c, int) and c >= 0 for c in f.values())
    coefficients = [f.get(b + q * j, 0) for j in range(count)]
    assert coefficients[0] == 1

    # Independent, reduced-coordinate equation printed in Lemma 7.1.
    # B(z) = (1 + z B(z)^b B(z^b B(z)^(b^2-1)))^b.
    small_cap = count - 1
    def reduced_step(series):
        inner = mul({b: 1}, power(series, q, small_cap), small_cap)
        product = mul(power(series, b, small_cap),
                      compose(series, inner, small_cap), small_cap)
        return power(add({0: 1}, mul(x, product, small_cap), small_cap),
                     b, small_cap)
    reduced, reduced_iterations = fixed_point(reduced_step, {0: 1}, count + 3)
    assert coefficients == [reduced.get(j, 0) for j in range(count)]

    divisors = [(b - 1) * j + 1 for j in range(count)]
    assert all(a % d == 0 for a, d in zip(coefficients, divisors))
    quotients = [a // d for a, d in zip(coefficients, divisors)]

    # Rational primitive P(0)=0, P'=(b+1) f.
    primitive = {k + 1: Fraction((b + 1) * c, k + 1)
                 for k, c in f.items() if k + 1 <= cap}
    assert all(c.denominator == 1 for c in primitive.values())
    g = add(x, compose(f, f, cap), cap)
    rhs = add(scale(mul(x, f, cap), b + 1),
              scale(power(g, b + 1, cap), -b), cap)
    lhs = add(primitive, scale(compose(primitive, f, cap), -1), cap)
    assert lhs == rhs
    for j, quotient in enumerate(quotients):
        s = (b - 1) * j + 1
        assert primitive[(b + 1) * s] == quotient
    return dict(b=b, supported_coefficients=count, checked_through_degree=cap,
                full_iterations=iterations, reduced_iterations=reduced_iterations,
                coefficients=coefficients, quotients=quotients,
                divisibility=True, integral_primitive=True, coboundary=True,
                independent_coordinate_crosscheck=True)


def main():
    rows = [run_case(b, 16) for b in range(2, 9)]
    # Fidelity: every explicit coefficient in Remark 7.2, with its indexing.
    assert rows[0]['coefficients'][:12] == [
        1, 2, 9, 56, 400, 3096, 25256, 213832, 1861272,
        16552320, 149737632, 1373597892]
    assert rows[1]['coefficients'][:11] == [
        1, 3, 30, 406, 6336, 107415, 1922310, 35739990,
        683593902, 13364444808, 265869803598]
    print('FIDELITY_OK: all 12 b=2 and all 11 b=3 displayed coefficients, Remark 7.2')
    for row in rows:
        print(json.dumps(row, sort_keys=True))
    print('CHECKED: b=2..8; n=1..16; 112 divisibilities; both coordinate constructions;')
    print('CHECKED: integral primitive and P-P(f)=(b+1)xf-b(x+f(f))^(b+1) at every retained degree')
    print('SCOPE: exact finite evidence; no proof for unbounded b or n; no Lean proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
