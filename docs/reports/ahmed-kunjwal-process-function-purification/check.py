#!/usr/bin/env python3
"""Exact small-instance checks for Ahmed--Kunjwal Conjecture VII.1.

All arithmetic is over integers/rationals (SymPy exact matrices); no floating
point value is used as evidence.
"""
from itertools import product
from sympy import Matrix, Rational, I, zeros, eye, simplify


def tuples(dims):
    return list(product(*[range(d) for d in dims]))


def local_maps(input_dims, output_dims):
    choices = [list(product(range(output_dims[k]), repeat=input_dims[k])) for k in range(len(input_dims))]
    for fs in product(*choices):
        yield fs


def apply_local(fs, i):
    return tuple(fs[k][i[k]] for k in range(len(fs)))


def is_process(omega, input_dims, output_dims):
    inputs = tuples(input_dims)
    outputs = tuples(output_dims)
    for fs in local_maps(input_dims, output_dims):
        fixed = []
        for i in inputs:
            o = apply_local(fs, i)
            if omega[o] == i:
                fixed.append(i)
        if len(fixed) != 1:
            return False
    return True


def process_functions(input_dims, output_dims):
    outputs = tuples(output_dims)
    inputs = tuples(input_dims)
    for values in product(inputs, repeat=len(outputs)):
        omega = dict(zip(outputs, values))
        if is_process(omega, input_dims, output_dims):
            yield omega


def add_mod(x, y, dims):
    return tuple((x[k] + y[k]) % dims[k] for k in range(len(dims)))


def sub_mod(x, y, dims):
    return tuple((x[k] - y[k]) % dims[k] for k in range(len(dims)))


def reversible_extension(omega, input_dims, output_dims):
    def u(o, s):
        return add_mod(s, omega[o], input_dims), o
    def inverse(i, o):
        return o, sub_mod(i, omega[o], input_dims)
    return u, inverse


def check_bijection(omega, input_dims, output_dims):
    outputs = tuples(output_dims)
    inputs = tuples(input_dims)
    u, inv = reversible_extension(omega, input_dims, output_dims)
    domain = [(o, s) for o in outputs for s in inputs]
    codomain = [(i, o) for i in inputs for o in outputs]
    assert {u(*x) for x in domain} == set(codomain)
    assert all(inv(*u(*x)) == x for x in domain)
    assert all(u(*inv(*y)) == y for y in codomain)


def purification_check(omega, input_dims, output_dims):
    outputs = tuples(output_dims)
    zero = tuple(0 for _ in input_dims)
    u, _ = reversible_extension(omega, input_dims, output_dims)
    for o in outputs:
        assert u(o, zero) == (omega[o], o)
        for i in tuples(input_dims):
            assert ((i == omega[o]) == (u(o, zero) == (i, o)))


def L_omega(omega, matrices, input_dims, output_dims):
    outputs = tuples(output_dims)
    total = 0
    for o in outputs:
        term = 1
        i = omega[o]
        for k in range(len(input_dims)):
            term *= matrices[k][o[k], i[k]]
        total += term
    return simplify(total)


def column_sum_matrix(rows, cols, entries):
    m = Matrix(entries)
    sums = [simplify(sum(m[r, c] for r in range(rows))) for c in range(cols)]
    assert all(s == sums[0] for s in sums)
    return m, sums[0]


def normalization_check(omega, input_dims, output_dims):
    # Exact complex matrices with constant column sums, including non-real and
    # negative entries.  This is the source's nontrivial multilinear step.
    mats = []
    cs = []
    for k, (di, do) in enumerate(zip(input_dims, output_dims)):
        if (di, do) == (2, 2):
            m, c = column_sum_matrix(2, 2, [[1 + I, 2], [3, 2 + I]])
        elif (di, do) == (1, 1):
            m, c = column_sum_matrix(1, 1, [[2 - I]])
        else:
            # A rectangular exact test: every column has the same sum.
            rows = [[Rational(r + 1) + I * Rational(c + 1) for c in range(di)] for r in range(do)]
            m, c = column_sum_matrix(do, di, rows)
        mats.append(m)
        cs.append(c)
    got = L_omega(omega, mats, input_dims, output_dims)
    expected = 1
    for c in cs:
        expected *= c
    assert simplify(got - expected) == 0, (omega, got, expected)


def linked_kraus(omega, kraus_families, input_dims, output_dims):
    """The no-ancilla linked Kraus matrices from Appendix A, Eq. (33)."""
    outputs = tuples(output_dims)
    inputs = tuples(input_dims)
    out = []
    for choices in product(*[range(len(Ks)) for Ks in kraus_families]):
        K = zeros(len(outputs), len(inputs))
        for ro, o in enumerate(outputs):
            for ci, s in enumerate(inputs):
                value = 1
                i = add_mod(s, omega[o], input_dims)
                for k, Ks in enumerate(kraus_families):
                    value *= Ks[choices[k]][o[k], i[k]]
                K[ro, ci] = simplify(value)
        out.append(K)
    return out


def kraus_completeness(Ks):
    total = zeros(Ks[0].cols, Ks[0].cols)
    for K in Ks:
        total += K.conjugate().T * K
    return total


def validity_check(omega, input_dims, output_dims):
    # Identity/unitary and exact dephasing channels are CPTP; the latter tests
    # off-diagonal Gram coefficients rather than only deterministic maps.
    local_families = []
    for di, do in zip(input_dims, output_dims):
        assert di == do
        if di == 1:
            local_families.append([Matrix([[1]])])
        elif di == 2:
            local_families.append([Matrix([[1, 0], [0, 0]]), Matrix([[0, 0], [0, 1]])])
        else:
            local_families.append([eye(di)])
    for fam in local_families:
        assert kraus_completeness(fam) == eye(fam[0].cols)
    linked = linked_kraus(omega, local_families, input_dims, output_dims)
    assert kraus_completeness(linked) == eye(len(tuples(input_dims)))


def run():
    # Fidelity: the AF/BW three-party binary process and the explicit reversible
    # extension in Appendix A are checked on every basis state.
    def afbw(o):
        o1, o2, o3 = o
        return ((o2 ^ 1) * o3, (o3 ^ 1) * o1, (o1 ^ 1) * o2)
    af = {o: afbw(o) for o in tuples((2, 2, 2))}
    assert is_process(af, (2, 2, 2), (2, 2, 2))
    check_bijection(af, (2, 2, 2), (2, 2, 2))
    purification_check(af, (2, 2, 2), (2, 2, 2))

    checked = 0
    cases = [((1, 1), (1, 1), None), ((2, 2), (2, 2), None)]
    # Exhaustive enumeration is feasible for one and two parties.  For three
    # binary parties use the published AF/BW process plus several distinct
    # causal process functions; this keeps the check exact without claiming a
    # finite enumeration is a proof of the universal statement.
    causal3 = []
    for mode in range(4):
        def make_causal(mode):
            def w(o):
                o1, o2, _ = o
                if mode == 0:
                    return (0, 0, 0)
                if mode == 1:
                    return (0, o1, 0)
                if mode == 2:
                    return (0, o1, o1 ^ o2)
                return (0, o1, o2)
            return w
        w = make_causal(mode)
        causal3.append({o: w(o) for o in tuples((2, 2, 2))})
    pfs3 = [af] + causal3
    assert all(is_process(w, (2, 2, 2), (2, 2, 2)) for w in pfs3)
    cases.append(((2, 2, 2), (2, 2, 2), pfs3))
    for dims in cases:
        input_dims, output_dims, supplied = dims
        pfs = list(process_functions(input_dims, output_dims)) if supplied is None else supplied
        assert pfs
        for omega in pfs:
            check_bijection(omega, input_dims, output_dims)
            purification_check(omega, input_dims, output_dims)
            normalization_check(omega, input_dims, output_dims)
            validity_check(omega, input_dims, output_dims)
            checked += 1
        print(f"checked process functions for dimensions {input_dims}->{output_dims}: {len(pfs)}")

    # A rectangular/singleton-alphabet fidelity case: one party with a
    # singleton output and a singleton input is still covered exactly.
    singleton = {(0,): (0,)}
    check_bijection(singleton, (1,), (1,))
    purification_check(singleton, (1,), (1,))
    print(f"checked total process functions: {checked + 1}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    run()
