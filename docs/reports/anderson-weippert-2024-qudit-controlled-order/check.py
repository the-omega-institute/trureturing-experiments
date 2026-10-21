"""Numerical check of the conventions for Anderson-Weippert (arXiv:2410.04711v3, Section 4).

Pauli group: c * X^a Z^b on n qudits of prime dimension p, |c| = 1.
Clifford hierarchy: CH_1 = Pauli, CH_k = {V unitary : V P V^dagger in CH_{k-1} for every Pauli P}.
Controlled gate: C(U) = sum_j |j><j| (x) U_j with U_1 = U and U_j = I otherwise (control first).

The script computes the least level k <= KMAX with C(U) in CH_k, and the multiplicative order of U, for qubit gates
whose levels the paper states (Section 3: C(Z), C(S), C(T), C(H); Corollary 2.4.1: U of order not a power of 2 is
excluded) and for qutrit gates (the open problem: order a power of 3).
"""
import itertools
import sys

import numpy as np

TOL = 1e-9


def digits(p, n):
    return list(itertools.product(range(p), repeat=n))


def weyl(p, n, a, b):
    w = np.exp(2j * np.pi / p)
    basis = digits(p, n)
    index = {v: i for i, v in enumerate(basis)}
    m = np.zeros((p ** n, p ** n), dtype=complex)
    for j in basis:
        i = tuple((j[k] + a[k]) % p for k in range(n))
        m[index[i], index[j]] = w ** (sum(b[k] * j[k] for k in range(n)) % p)
    return m


def paulis(p, n):
    return [weyl(p, n, a, b) for a in digits(p, n) for b in digits(p, n)]


def is_pauli(m, plist):
    for w in plist:
        k = np.argmax(np.abs(w[:, 0]))
        c = m[k, 0] / w[k, 0]
        if abs(abs(c) - 1) < 1e-7 and np.allclose(m, c * w, atol=1e-7):
            return True
    return False


def in_ch(v, k, plist):
    if k == 1:
        return is_pauli(v, plist)
    if not np.allclose(v @ v.conj().T, np.eye(v.shape[0]), atol=TOL):
        return False
    vd = v.conj().T
    return all(in_ch(v @ q @ vd, k - 1, plist) for q in plist)


def controlled(p, u):
    d = u.shape[0]
    blocks = [u if j == 1 else np.eye(d) for j in range(p)]
    m = np.zeros((p * d, p * d), dtype=complex)
    for j, blk in enumerate(blocks):
        m[j * d:(j + 1) * d, j * d:(j + 1) * d] = blk
    return m


def level(p, u, kmax):
    c = controlled(p, u)
    n = int(round(np.log(c.shape[0]) / np.log(p)))
    plist = paulis(p, n)
    for k in range(1, kmax + 1):
        if in_ch(c, k, plist):
            return k
    return None


def order(u, cap=200):
    m = np.eye(u.shape[0], dtype=complex)
    for r in range(1, cap + 1):
        m = m @ u
        if np.allclose(m, np.eye(u.shape[0]), atol=1e-8):
            return r
    return None


def is_power(x, p):
    if x is None:
        return False
    while x % p == 0:
        x //= p
    return x == 1


def main():
    s2 = 1 / np.sqrt(2)
    H = np.array([[s2, s2], [s2, -s2]], dtype=complex)
    S = np.diag([1, 1j])
    qubit = [
        ("X", np.array([[0, 1], [1, 0]], dtype=complex), 2),
        ("Z", np.diag([1, -1]).astype(complex), 2),
        ("S", S, 3),
        ("T", np.diag([1, np.exp(1j * np.pi / 4)]), 4),
        ("H", H, 3),
        ("iX", 1j * np.array([[0, 1], [1, 0]], dtype=complex), 2),
        ("SH", S @ H, None),
        ("Rz(1)", np.diag([1, np.exp(1j)]), None),
    ]
    w3 = np.exp(2j * np.pi / 3)
    qutrit = [
        ("Z3", np.diag([1, w3, w3 ** 2]), None),
        ("X3", np.roll(np.eye(3), 1, axis=0).astype(complex), None),
        ("diag(1,1,e^{2pi i/9})", np.diag([1, 1, np.exp(2j * np.pi / 9)]), None),
        ("e^{2pi i/9} I", np.exp(2j * np.pi / 9) * np.eye(3), None),
        ("-I", -np.eye(3, dtype=complex), None),
        ("diag(1,1,-1)", np.diag([1, 1, -1]).astype(complex), None),
    ]
    bad = 0
    print("p=2 (qubit), control + 1 target, KMAX=4")
    for name, u, expected in qubit:
        lv = level(2, u, 4)
        o = order(u)
        ok = lv == expected
        bad += not ok
        print(f"  C({name}): level={lv} expected={expected} order(U)={o} power_of_2={is_power(o, 2)} {'OK' if ok else 'MISMATCH'}")
    print("p=3 (qutrit), control + 1 target, KMAX=3")
    for name, u, _ in qutrit:
        lv = level(3, u, 3)
        o = order(u)
        consistent = lv is None or is_power(o, 3)
        bad += not consistent
        print(f"  C({name}): level={lv} order(U)={o} power_of_3={is_power(o, 3)} claim_consistent={consistent}")
    print(f"EXIT={1 if bad else 0}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
