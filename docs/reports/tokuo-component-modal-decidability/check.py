#!/usr/bin/env python3
"""Exact finite orthoframe and filtration checks; not a proof of decidability.

Source: Kenji Tokuo, Component Modalities of Quantum Logic, arXiv:2607.14757v1,
Sections 2 and 4. All sets, graphs and truth values are integer bit masks.
The test spaces and success conditions are fixed in main(). No randomness.
"""
from itertools import combinations, product


def frames(n):
    pairs = list(combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        rows = [1 << i for i in range(n)]
        for k, (i, j) in enumerate(pairs):
            if bits & (1 << k):
                rows[i] |= 1 << j
                rows[j] |= 1 << i
        yield tuple(rows)


def perp(rows, x):
    return sum(1 << w for w, neighbors in enumerate(rows) if not neighbors & x)


def components(rows):
    unseen = (1 << len(rows)) - 1
    result = []
    while unseen:
        reached = unseen & -unseen
        while True:
            expanded = reached
            for w, neighbors in enumerate(rows):
                if reached & (1 << w):
                    expanded |= neighbors
            if expanded == reached:
                break
            reached = expanded
        result.append(reached)
        unseen &= ~reached
    return tuple(result)


def box_component(rows, x):
    return sum(c for c in components(rows) if c & x == c)


def eval_formula(rows, values, formula):
    op = formula[0]
    if op == 'p':
        return values[formula[1]]
    if op == '&':
        return eval_formula(rows, values, formula[1]) & eval_formula(rows, values, formula[2])
    if op == '~':
        return perp(rows, eval_formula(rows, values, formula[1]))
    assert op == 'B'
    return box_component(rows, eval_formula(rows, values, formula[1]))


def subformulas(f):
    result = {f}
    if f[0] in ('~', 'B'):
        result |= subformulas(f[1])
    elif f[0] == '&':
        result |= subformulas(f[1]) | subformulas(f[2])
    return result


def atoms(f):
    return {x[1] for x in subformulas(f) if x[0] == 'p'}


def quotient(rows, values, formula):
    sub = subformulas(formula)
    occurring = atoms(formula)
    sigma = sorted(sub | {('~', ('p', p)) for p in occurring}, key=repr)
    truths = [eval_formula(rows, values, f) for f in sigma]
    signatures = [tuple(bool(x & (1 << w)) for x in truths) for w in range(len(rows))]
    names = {s: i for i, s in enumerate(sorted(set(signatures)))}
    q = tuple(names[s] for s in signatures)
    size = len(names)
    image_rows = [0] * size
    image_values = [0] * len(values)
    for i, neighbors in enumerate(rows):
        for j in range(len(rows)):
            if neighbors & (1 << j):
                image_rows[q[i]] |= 1 << q[j]
        for p, x in enumerate(values):
            if x & (1 << i):
                image_values[p] |= 1 << q[i]
    # Atoms absent from the input are set to the stable empty set.
    image_values = tuple(x if p in occurring else 0 for p, x in enumerate(image_values))
    assert size <= 1 << (len(sub) + len(occurring))
    assert all(image_rows[i] & (1 << i) for i in range(size))
    assert all(bool(image_rows[i] & (1 << j)) == bool(image_rows[j] & (1 << i))
               for i in range(size) for j in range(size))
    assert len(components(image_rows)) == 1
    for x in image_values:
        assert perp(image_rows, perp(image_rows, x)) == x
    for f in sub:
        old = eval_formula(rows, values, f)
        new = eval_formula(image_rows, image_values, f)
        assert all(bool(old & (1 << w)) == bool(new & (1 << q[w]))
                   for w in range(len(rows))), (rows, values, formula, f)
    return size < len(rows)


def formula_space():
    previous = {('p', 0), ('p', 1)}
    for _ in range(2):
        previous = ({('p', 0), ('p', 1)} |
                    {(op, f) for op in ('~', 'B') for f in previous} |
                    {('&', f, g) for f in previous for g in previous})
    return sorted(previous, key=repr)


def main():
    # Fidelity: source's Modal Alternatives across Components example (Section 2).
    rows = (1, 2)
    values = (1,)
    p = ('p', 0)
    bx = eval_formula(rows, values, ('B', p))
    nbx = eval_formula(rows, values, ('~', ('B', p)))
    assert (bx, nbx) == (1, 2) and bx | nbx == 3 and bx != 3 and nbx != 3
    print('FIDELITY Section 2 modal alternatives: box(p)={0}; neg(box(p))={1}; local MEM=true; neither conclusion globally true')

    # Fidelity: Theorem 4.2's component approximation identities, all graphs n<=5.
    graph_count = stable_count = 0
    for n in range(1, 6):
        top = (1 << n) - 1
        for rows in frames(n):
            graph_count += 1
            cs = components(rows)
            for x in range(1 << n):
                assert x & ~perp(rows, perp(rows, x)) == 0
                if perp(rows, perp(rows, x)) != x:
                    continue
                stable_count += 1
                bx = box_component(rows, x)
                dx = perp(rows, box_component(rows, perp(rows, x)))
                assert bx == sum(c for c in cs if c & x == c)
                assert dx == sum(c for c in cs if c & x)
                assert bx & ~x == 0 and x & ~dx == 0
                assert perp(rows, perp(rows, bx)) == bx
                assert bx | perp(rows, bx) == top
    print(f'FIDELITY Theorem 4.2: frames={graph_count}; stable_sets={stable_count}; all n=1..5 component interior/cover identities hold')

    # Settlement check: every connected graph n<=4, every stable two-atom valuation,
    # every formula of modal/negation/conjunction depth<=2 over those atoms.
    formulas = formula_space()
    model_count = check_count = strict_count = 0
    for n in range(1, 5):
        for rows in frames(n):
            if len(components(rows)) != 1:
                continue
            stable_sets = [x for x in range(1 << n) if perp(rows, perp(rows, x)) == x]
            for values in product(stable_sets, repeat=2):
                model_count += 1
                for formula in formulas:
                    strict_count += quotient(rows, values, formula)
                    check_count += 1
    assert strict_count > 0
    print(f'FILTRATION: formulas={len(formulas)}; models={model_count}; quotients={check_count}; strict_quotients={strict_count}; all truth/stability/connectivity/bound checks hold')
    print('SCOPE: finite exact checks only; canonical completeness and the unbounded filtration theorem remain proof obligations')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
