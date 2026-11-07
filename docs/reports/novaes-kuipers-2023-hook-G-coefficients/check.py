#!/usr/bin/env python3
"""Exact LR-tableau checks of Novaes--Kuipers arXiv:2302.09030v1, Eq. (37).

All coefficients are Python integers. No Sage, floating-point approximation,
closed-form LR evaluation, or claimed infinite proof is used.
"""
from collections import defaultdict
from functools import lru_cache


@lru_cache(None)
def partitions(n, ceiling=None):
    if n == 0:
        return ((),)
    ceiling = min(n, n if ceiling is None else ceiling)
    return tuple((first,) + tail for first in range(ceiling, 0, -1)
                 for tail in partitions(n - first, first))


def contains(outer, inner):
    return len(inner) <= len(outer) and all(a <= b for a, b in zip(inner, outer))


@lru_cache(None)
def lr(left, content, outer):
    """Count semistandard skew tableaux of shape outer/left and given content.

    Read rows top to bottom, each right to left. Every prefix must be a lattice
    word: #i >= #(i+1). Entries are 0-based; content[i] counts entry i.
    """
    if sum(left) + sum(content) != sum(outer) or not contains(outer, left):
        return 0
    cells = tuple((row, col) for row, width in enumerate(outer)
                  for col in range(width - 1, -1, -1)
                  if col >= (left[row] if row < len(left) else 0))
    used = [0] * len(content)
    filled = {}

    def visit(pos):
        if pos == len(cells):
            return int(tuple(used) == content)
        row, col = cells[pos]
        answer = 0
        for label in range(len(content)):
            if used[label] == content[label]:
                continue
            if (row, col + 1) in filled and label > filled[row, col + 1]:
                continue
            if (row - 1, col) in filled and label <= filled[row - 1, col]:
                continue
            used[label] += 1
            if label == 0 or used[label] <= used[label - 1]:
                filled[row, col] = label
                answer += visit(pos + 1)
                del filled[row, col]
            used[label] -= 1
        return answer

    return visit(0)


def hook(arm, leg):
    assert arm >= 1 and leg >= 0
    return (arm,) + (1,) * leg


def G(rho, nu):
    """Eq. (37), returned as {power of t: integer coefficient}, t = sqrt(gamma)."""
    n, m = sum(nu), sum(rho)
    answer = defaultdict(int)
    if m > n:
        return {}
    for k in range(n + 1):
        alpha_size = n - m - k
        if alpha_size < 0:
            continue
        column = (1,) * k
        for mu in partitions(m + k):
            first = lr(column, rho, mu)
            if first:
                for alpha in partitions(alpha_size):
                    second = lr(mu, alpha, nu)
                    answer[k + alpha_size] += (-1) ** k * first * second
    return {degree: coefficient for degree, coefficient in answer.items() if coefficient}


def expected(R, r, A, a):
    if A == R and a >= r and (a - r) % 2 == 0:
        return {a - r: 1}
    if A == R + 1 and a >= r + 1 and (a - r - 1) % 2 == 0:
        return {a + 1 - r: 1}
    return {}


def polynomial(p):
    if not p:
        return '0'
    assert len(p) == 1 and next(iter(p.values())) == 1
    power = next(iter(p))
    return '1' if power == 0 else f't^{power}'


def main():
    nu = (3, 1, 1, 1)
    hooks6 = [hook(R, n - R) for n in range(1, 7) for R in range(1, n + 1)]
    observed = {rho: G(rho, nu) for rho in hooks6 if G(rho, nu)}
    stated = {(3, 1, 1, 1): {0: 1}, (2, 1, 1): {2: 1},
              (3, 1): {2: 1}, (2,): {4: 1}}
    assert observed == stated, (observed, stated)
    print('FIDELITY Section II.6 example, nu=(3,1,1,1):')
    for rho, value in stated.items():
        print(f'  rho={rho}: G={polynomial(value)}')
    cases = nonzero = 0
    by_target_size = defaultdict(lambda: [0, 0])
    hooks8 = [(R, n - R) for n in range(1, 9) for R in range(1, n + 1)]
    for R, r in hooks8:
        for A, a in hooks8:
            got, want = G(hook(R, r), hook(A, a)), expected(R, r, A, a)
            assert got == want, (R, r, A, a, got, want)
            cases += 1
            nonzero += bool(got)
            by_target_size[A + a][0] += 1
            by_target_size[A + a][1] += bool(got)
    print('EXACT_IDENTITY all ordered pairs of nonempty hooks, each size <= 8:')
    print(f'  hooks={len(hooks8)}, ordered_pairs={cases}, nonzero={nonzero}, zero={cases-nonzero}')
    for n, (count, nz) in sorted(by_target_size.items()):
        print(f'  target_size={n}: checked={count}, nonzero={nz}')
    print(f'LR_TABLEAU_CACHE unique_triples={lr.cache_info().currsize}')
    print('All sums and tableau counts used exact integer arithmetic.')
    print('Finite checks do not prove the all-partition assertion.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
