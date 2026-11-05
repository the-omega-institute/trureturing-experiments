#!/usr/bin/env python3
"""Exact finite checks for Moreto--Navarro, arXiv:2301.03917v2, Conjecture 2.6.

No floating point or external algebra package is used. Roots of unity are
represented in Z[X]/Phi_(3^a)(X); all reductions use integer coefficients.
Presentations, rather than isomorphism classes, are enumerated.
"""
from collections import Counter, deque
from itertools import product
from math import prod


def partitions(n, largest=None):
    if n == 0:
        yield ()
        return
    for k in range(min(n, largest or n), 0, -1):
        for tail in partitions(n - k, k):
            yield (k,) + tail


class Abelian:
    def __init__(self, moduli):
        self.moduli = tuple(moduli)
        self.elements = tuple(product(*(range(q) for q in moduli)))
        self.zero = tuple(0 for _ in moduli)
        self.order = prod(moduli)
        self.exponent = max(moduli)
        self.basis = tuple(tuple(int(i == j) for i in range(len(moduli)))
                           for j in range(len(moduli)))

    def add(self, x, y):
        return tuple((a + b) % q for a, b, q in zip(x, y, self.moduli))

    def neg(self, x):
        return tuple(-a % q for a, q in zip(x, self.moduli))

    def sub(self, x, y):
        return self.add(x, self.neg(y))

    def scale(self, k, x):
        return tuple(k * a % q for a, q in zip(x, self.moduli))

    def hom(self, columns):
        return {x: tuple(sum(a * c[i] for a, c in zip(x, columns)) % q
                         for i, q in enumerate(self.moduli)) for x in self.elements}

    def exponent_of_value(self, character, x):
        # lambda(x) = zeta_exponent ^ this integer; all linear characters occur.
        return sum(k * a * (self.exponent // q)
                   for k, a, q in zip(character, x, self.moduli)) % self.exponent


def cyclotomic_reduce(exponents, q):
    """Reduce sum X^e exactly modulo Phi_q, q a positive power of three."""
    assert q >= 3
    third = q // 3
    coefficients = Counter()
    for exponent in exponents:
        e = exponent % q
        if e < 2 * third:
            coefficients[e] += 1
        else:
            j = e - 2 * third
            coefficients[j] -= 1
            coefficients[j + third] -= 1
    return {e: c for e, c in coefficients.items() if c}


def actions(A):
    """All automorphisms alpha of this labelled abelian group with alpha^3=1."""
    choices = [tuple(x for x in A.elements if A.scale(q, x) == A.zero)
               for q in A.moduli]
    for columns in product(*choices):
        alpha = A.hom(columns)
        if len(set(alpha.values())) != A.order:
            continue
        if all(alpha[alpha[alpha[x]]] == x for x in A.elements):
            yield alpha


def structure(A, alpha):
    delta = {u: A.sub(alpha[u], u) for u in A.elements}
    kernel = {u for u in A.elements if delta[u] == A.zero}
    image = set(delta.values())
    c = 1
    while image != {A.zero}:
        image = {delta[u] for u in image}
        c += 1
        assert c <= A.order
    return delta, kernel, c


def character_checks(A, alpha, delta, center, nilclass, n):
    """Every center-faithful lambda, not one selected successful character."""
    assert nilclass == n - 1 and len(center) == 3
    ker2 = {u for u in A.elements if delta[delta[u]] == A.zero}
    assert len(ker2) == 9
    desired = ker2 - center
    count = 0
    first_zero = None
    for lam in A.elements:
        if all(A.exponent_of_value(lam, z) == 0 for z in center):
            continue
        count += 1
        orbit = [tuple(A.exponent_of_value(lam, a[u]) for u in A.elements)
                 for a in ({u: u for u in A.elements}, alpha,
                           {u: alpha[alpha[u]] for u in A.elements})]
        assert len(set(orbit)) == 3
        zero = set()
        norm_terms = Counter()
        for j, u in enumerate(A.elements):
            exps = tuple(row[j] for row in orbit)
            if not cyclotomic_reduce(exps, A.exponent):
                zero.add(u)
            norm_terms.update((a - b) % A.exponent for a in exps for b in exps)
        # Sum |chi(u)|^2 over U equals |G|; outside U the induced trace is zero.
        norm = cyclotomic_reduce((e for e, mult in norm_terms.items()
                                  for _ in range(mult)), A.exponent)
        assert norm == {0: 3 * A.order}
        assert zero == desired
        assert len(zero) == 6
        assert 2 * A.order + len(zero) == 3**n - 3**(n - 1) + 6
        if first_zero is None:
            first_zero = zero
    assert count > 0
    return count, first_zero


def extension_checks(A, alpha, t, zero_in_U, n):
    """Concrete G=U.<g>, g u g^-1=alpha(u), g^3=t; t need not be zero."""
    assert alpha[t] == t
    powers = [{u: u for u in A.elements}, alpha,
              {u: alpha[alpha[u]] for u in A.elements}]

    def mul(x, y):
        u, i = x
        v, j = y
        return (A.add(A.add(u, powers[i][v]), A.scale((i + j) // 3, t)),
                (i + j) % 3)

    elements = tuple(product(A.elements, range(3)))
    identity = (A.zero, 0)
    inverses = {}
    for x in elements:
        u, i = x
        j = -i % 3
        v = powers[j][A.neg(A.add(u, A.scale((i + j) // 3, t)))]
        inverses[x] = (v, j)
        assert mul(x, inverses[x]) == identity == mul(inverses[x], x)
    # The cocycle identity and the action law ensure associativity, checked
    # independently on all carry indices and all vectors, without cubic G loops.
    for i, j, k in product(range(3), repeat=3):
        left = A.add(A.scale((i + j) // 3, t),
                     A.scale(((i + j) % 3 + k) // 3, t))
        right = A.add(powers[i][A.scale((j + k) // 3, t)],
                      A.scale((i + (j + k) % 3) // 3, t))
        assert left == right
    for i, j in product(range(3), repeat=2):
        assert all(powers[i][powers[j][u]] == powers[(i + j) % 3][u]
                   for u in A.elements)
    generators = [(b, 0) for b in A.basis] + [(A.zero, 1)]
    unseen = set(elements)
    zero_classes = 0
    zero_elements = 0
    while unseen:
        x = min(unseen)
        orbit, queue = {x}, deque([x])
        while queue:
            y = queue.popleft()
            for s in generators:
                z = mul(mul(s, y), inverses[s])
                if z not in orbit:
                    orbit.add(z)
                    queue.append(z)
        unseen.difference_update(orbit)
        membership = {i != 0 or u in zero_in_U for u, i in orbit}
        assert len(membership) == 1
        if membership == {True}:
            zero_classes += 1
            zero_elements += len(orbit)
    # Source Theorem 2.7 with p=3: equality in elements iff eight zero classes.
    assert zero_classes == 8
    assert zero_elements == 3**n - 3**(n - 1) + 6
    norm_image = {A.add(A.add(u, alpha[u]), alpha[alpha[u]]) for u in A.elements}
    return t not in norm_image


def literal_boundary():
    # C9: all nine ordinary irreducibles are lambda_k(j)=zeta_9^(kj).
    # Every value is nonzero; this uses exact cyclotomic reduction as above.
    for k, j in product(range(9), repeat=2):
        assert cyclotomic_reduce([k * j], 9)
    G = set(range(9))
    center = {x for x in G if all((x + y) % 9 == (y + x) % 9 for y in G)}
    assert center == G and len(G) > 1  # class exactly 1 = n-1 for n=2
    U = {0, 3, 6}
    assert all((x + y) % 9 in U for x, y in product(U, repeat=2))
    # Check maximality directly: adjoining any element outside U generates C9.
    assert all({(k * x + u) % 9 for k in range(9) for u in U} == G for x in G - U)
    target = 3**2 - 3**(2 - 1) + 6
    assert target == 12 > len(G)
    print('LITERAL_BOUNDARY C9: order=9 class=1 abelian_maximal_order=3 target=12; literal equivalence is false')


def main():
    literal_boundary()
    total_actions = total_extensions = total_characters = nonsplit = 0
    for n in (3, 4):
        found = extensions = chars = ns = 0
        for partition in partitions(n - 1):
            A = Abelian(tuple(3**k for k in partition))
            for alpha in actions(A):
                delta, center, c = structure(A, alpha)
                if c != n - 1:
                    continue
                found += 1
                char_count, zeros = character_checks(A, alpha, delta, center, c, n)
                chars += char_count
                for t in sorted(center):
                    extensions += 1
                    ns += extension_checks(A, alpha, t, zeros, n)
        assert found and ns
        print(f'EXHAUSTIVE n={n}: labelled_actions={found} extensions={extensions} '
              f'nonsplit_relative_to_U={ns} center_faithful_characters_per_action_total={chars}; '
              'zero_set=ker(delta^2)-ker(delta), internal_zeros=6, zero_classes=8')
        total_actions += found
        total_extensions += extensions
        total_characters += chars
        nonsplit += ns
    # Larger exact examples: additive Z[omega]/(omega-1)^m, alpha=omega.
    for m in (4, 5, 6):
        k = m // 2
        if m % 2 == 0:
            A = Abelian((3**k, 3**k))
            alpha = {u: (-u[1] % A.moduli[0], (u[0] - u[1]) % A.moduli[1])
                     for u in A.elements}
        else:
            A = Abelian((3**(k + 1), 3**k))
            alpha = {u: ((u[0] - 3*u[1]) % A.moduli[0], (u[0] - 2*u[1]) % A.moduli[1])
                     for u in A.elements}
        assert all(alpha[alpha[alpha[u]]] == u for u in A.elements)
        delta, center, c = structure(A, alpha)
        count, zeros = character_checks(A, alpha, delta, center, c, m + 1)
        norm_image = {A.add(A.add(u, alpha[u]), alpha[alpha[u]]) for u in A.elements}
        assert norm_image == {A.zero}
        # t=0 is split; either nonzero center t gives a nonsplit extension.
        assert len(center - norm_image) == 2
        total_characters += count
        print(f'EISENSTEIN n={m+1}: U_moduli={A.moduli} class={c} '
              f'center_faithful_characters={count} twist_parameters=3 '
              f'internal_zeros={len(zeros)} element_zeros={2*A.order+6}; '
              'all three twists have the same induced traces')
    print(f'SUMMARY small_labelled_actions={total_actions} small_extensions={total_extensions} '
          f'small_nonsplit={nonsplit} total_center_faithful_characters={total_characters}')
    print('FIDELITY source Theorem 2.7: checked eight zero classes versus 3^n-3^(n-1)+6 zero elements in every enumerated small extension')
    print('SCOPE all abelian-kernel cyclic-extension presentations at orders 27 and 81; selected Eisenstein actions at orders 243,729,2187; no universal proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
