#!/usr/bin/env python3
"""Exact finite algebra checks for GLL arXiv:2609.30257v1, Conjecture 7.2.5.

These computations test necessary homology identities in the proposed proof.
They do not verify any topological map, homotopy, or suspension section.
No third-party code or packages are used.
"""
from collections import Counter
from functools import lru_cache
from itertools import product


def clean(p):
    return Counter({w: c for w, c in p.items() if c})


def add(p, q):
    out = p.copy()
    out.update(q)
    return clean(out)


def mul(p, q):
    out = Counter()
    for u, a in p.items():
        for v, b in q.items():
            out[u + v] += a * b
    return clean(out)


def compositions(n):
    if n == 0:
        yield ()
    else:
        for a in range(1, n + 1):
            for b in compositions(n - a):
                yield (a,) + b


@lru_cache(None)
def ribbon(alpha):
    """Inclusion-exclusion in the free integer algebra on complete generators."""
    out = Counter()
    if not alpha:
        return Counter({(): 1})
    for mask in range(1 << (len(alpha) - 1)):
        beta, a = [], alpha[0]
        for i, b in enumerate(alpha[1:]):
            if mask & (1 << i):
                beta.append(a)
                a = b
            else:
                a += b
        beta.append(a)
        out[tuple(beta)] += (-1) ** (len(alpha) - len(beta))
    return clean(out)


def ribbon_word(alpha):
    return tuple(c for i, a in enumerate(alpha)
                 for c in ((1,) if i else ()) + (0,) * (a - 1))


@lru_cache(None)
def complete_finite(r, alphabet):
    return Counter({w: 1 for w in product(range(1, alphabet + 1), repeat=r)
                    if all(a <= b for a, b in zip(w, w[1:]))})


def evaluate_h(p, alphabet):
    out = Counter()
    for alpha, coefficient in p.items():
        term = Counter({(): 1})
        for a in alpha:
            term = mul(term, complete_finite(a, alphabet))
        for w, c in term.items():
            out[w] += coefficient * c
    return clean(out)


def direct_ribbon(alpha, alphabet):
    s = ribbon_word(alpha)
    return Counter({w: 1 for w in product(range(1, alphabet + 1), repeat=sum(alpha))
                    if all((a <= b if c == 0 else a > b)
                           for a, b, c in zip(w, w[1:], s))})


def fidelity():
    # Section 2.4: R_s R_t = R_{s0t} + R_{s1t}.
    pairs = 0
    for total in range(2, 11):
        for left in range(1, total):
            for a in compositions(left):
                for b in compositions(total - left):
                    fused = a[:-1] + (a[-1] + b[0],) + b[1:]
                    assert mul(ribbon(a), ribbon(b)) == add(ribbon(fused), ribbon(a + b))
                    pairs += 1
    # Independent SSYT/word realization of the inclusion-exclusion convention.
    words = 0
    for degree in range(1, 7):
        for a in compositions(degree):
            for alphabet in range(1, 4):
                assert evaluate_h(ribbon(a), alphabet) == direct_ribbon(a, alphabet)
                words += 1
    print(f'FIDELITY: ribbon product, total degree 2..10: {pairs} pairs')
    print(f'FIDELITY: independent finite-alphabet ribbon realization, degree 1..6: {words} cases')


def subsequences(s, limit):
    seen = {()}
    for c in s:
        seen |= {w + (c,) for w in seen if len(w) < limit}
    return seen


def stable_ranks():
    # Remark 7.2.2: enough alternations imply stability; degrees are half-degrees.
    for bound in range(1, 11):
        s = tuple(i % 2 for i in range(2 * bound + 1))
        ss = subsequences(s, bound - 1)
        for d in range(1, bound + 1):
            rank = sum(len(w) == d - 1 for w in ss)
            tensor_rank = sum(1 for _ in compositions(d))
            assert rank == tensor_rank == 2 ** (d - 1)
    print('FIDELITY/COMPARISON: alternating-source and tensor-algebra ranks agree, half-degree 0..10')


def words_up_to(n, length):
    return [w for k in range(length + 1) for w in product(range(1, n + 1), repeat=k)]


def small_quotient_homology():
    # Cellular H_*(CP^N): degrees 0,2,...,2N. James words delete 0-degree letters.
    # The selected sections below are homology-group sections only.
    stages = tuples = reduced_checks = 0
    for n in range(4):
        for m in range(5):
            for k in range(5):
                target = set(product(words_up_to(n, m), words_up_to(n, k)))
                actual_image = set()
                for a in product(range(n + 1), repeat=m):
                    for b in product(range(n + 1), repeat=k):
                        u, v = tuple(x for x in a if x), tuple(x for x in b if x)
                        assert sum(a) + sum(b) == sum(u) + sum(v)
                        actual_image.add((u, v))
                        tuples += 1
                assert actual_image == target
                for u, v in target:
                    a = u + (0,) * (m - len(u))
                    b = v + (0,) * (k - len(v))
                    assert tuple(x for x in a if x) == u
                    assert tuple(x for x in b if x) == v
                # Reduced homology on product and proposed wedge-smash indexing, degree by degree.
                reduced = Counter(sum(u) + sum(v) for u, v in target if u or v)
                smash = Counter()
                for i in range(m + 1):
                    for j in range(k + 1):
                        if i + j:
                            for w in product(range(1, n + 1), repeat=i + j):
                                smash[sum(w)] += 1
                assert reduced == smash
                stages += 1
                reduced_checks += sum(reduced.values())
    print(f'COMPARISON: N=0..3, m,n=0..4: {stages} stages, {tuples} source tuples')
    print(f'COMPARISON: quotient and wedge-smash reduced graded ranks: {reduced_checks} basis elements')


def coalgebra_and_order():
    # beta_r -> h_r respects Delta h_r = sum_{a+b=r} h_a tensor h_b, h_0=1.
    cases = 0
    for r in range(11):
        cp = Counter((a, r - a) for a in range(r + 1))
        complete = Counter((a, b) for a in range(r + 1) for b in range(r + 1) if a + b == r)
        assert cp == complete
        cases += 1
    # Preserve concatenation and order, not a commutative replacement.
    assert mul(Counter({(1,): 1}), Counter({(2,): 1})) != mul(
        Counter({(2,): 1}), Counter({(1,): 1}))
    for total in range(11):
        basis = list(compositions(total))
        for cut in range(total + 1):
            for a in compositions(cut):
                for b in compositions(total - cut):
                    assert mul(Counter({a: 1}), Counter({b: 1})) == Counter({a + b: 1})
                    cases += 1
        assert len(basis) == (1 if total == 0 else 2 ** (total - 1))
    print(f'COMPARISON: complete-generator coproduct and ordered tensor product: {cases} cases')


def main():
    fidelity()
    stable_ranks()
    small_quotient_homology()
    coalgebra_and_order()
    print('SCOPE: exact finite necessary algebra checks; no homotopies or topological equivalences certified')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
