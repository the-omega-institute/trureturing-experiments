#!/usr/bin/env python3
"""Qiu-Song-Chen, Phys. Rev. A 111 (2025) 022407 (arXiv:2410.15253v1), Sec. V, Conjecture 21:
K_{AD:BC}(F_4) = 2 ebits for the four-qubit Fredkin gate (swap C,D iff A = B = 1), with Definition 4:
sup over product pure inputs psi_i in H_{A_i R_i} of the entanglement entropy across the cut.
Numerical and exact checks (not a proof) of the counterexample family
  |a_e>_A |0>_{R_A} (x) |a_e>_B |0>_{R_B} (x) |Phi>_{C R_C} (x) |Phi>_{D R_D},  |a_e> = sqrt(e)|0> + sqrt(1-e)|1>:
  A. the reduced spectrum on {A R_A, D R_D} equals the closed form ((1+e+d)/4)^2, ((1+e-d)/4)^2, z, z, z,
     d = sqrt(1+10e-7e^2), z = (1-e)^2/4, and its entropy exceeds 2 bits for 0 < e <= 1/256 (e = 1/257 included);
  B. at e = 1/257 the coefficient matrix is rational: C = (1/257)[[1,16],[16,128]] plus three entries 128/257;
  C. the controlled-swap with n - 2 >= 3 controls split across the cut reproduces the same spectrum (n = 5, 6).
Exit 0 and final line `bad= 0` when every check agrees."""
import itertools, math, sys
import numpy as np
import mpmath as mp

mp.mp.dps = 60
bad = 0

def state_vector(controls_left, controls_right, eps_left, eps_right):
    """Qubits ordered: left controls (each with a reference), right controls (with references), C, R_C, D, R_D.
    Each control is a product sqrt(1-p)|0>+sqrt(p)|1> on its own qubit, reference |0>."""
    def ctrl(p):
        return np.array([math.sqrt(1 - p), math.sqrt(p)])
    phi = np.array([1, 0, 0, 1]) / math.sqrt(2)          # |Phi> on (target, reference)
    parts = []
    for p in eps_left + eps_right:
        parts.append(np.kron(ctrl(p), np.array([1.0, 0.0])))
    v = parts[0]
    for x in parts[1:]:
        v = np.kron(v, x)
    v = np.kron(v, np.kron(phi, phi))                    # C R_C D R_D
    return v

def apply_cswap(v, ncontrols):
    """Swap C and D when every control qubit is 1. Qubit layout: 2*ncontrols control/ref qubits, then C, R_C, D, R_D."""
    nq = 2 * ncontrols + 4
    out = np.zeros_like(v)
    for idx in range(2 ** nq):
        bits = [(idx >> (nq - 1 - k)) & 1 for k in range(nq)]
        if all(bits[2 * j] == 1 for j in range(ncontrols)):
            b = bits[:]
            c, d = 2 * ncontrols, 2 * ncontrols + 2
            b[c], b[d] = bits[d], bits[c]
            j = sum(bit << (nq - 1 - k) for k, bit in enumerate(b))
            out[j] += v[idx]
        else:
            out[idx] += v[idx]
    return out

def cut_spectrum(v, ncontrols, left_controls):
    """Left side: the first `left_controls` control parties (with refs) and D R_D; right: the rest and C R_C."""
    nq = 2 * ncontrols + 4
    T = v.reshape([2] * nq)
    left = []
    for j in range(left_controls):
        left += [2 * j, 2 * j + 1]
    left += [2 * ncontrols + 2, 2 * ncontrols + 3]
    right = [k for k in range(nq) if k not in left]
    M = np.transpose(T, left + right).reshape(2 ** len(left), 2 ** len(right))
    s = np.linalg.svd(M, compute_uv=False)
    return sorted([x * x for x in s if x * x > 1e-30], reverse=True)

def closed(e):
    e = mp.mpf(e)
    d = mp.sqrt(1 + 10 * e - 7 * e ** 2)
    z = (1 - e) ** 2 / 4
    return [((1 + e + d) / 4) ** 2, ((1 + e - d) / 4) ** 2, z, z, z]

def bits(sp):
    return -sum(x * mp.log(x, 2) for x in sp if x > 0)

# A: four-qubit Fredkin, controls A (left) and B (right)
for e in [mp.mpf(1) / 257, mp.mpf(1) / 256, mp.mpf(1) / 100, mp.mpf(1) / 50, mp.mpf(1) / 10]:
    v = state_vector(1, 1, [1 - float(e)], [1 - float(e)])   # control amplitude on |1> squared = 1 - e
    w = apply_cswap(v, 2)
    sp = cut_spectrum(w, 2, 1)
    cl = sorted(closed(e), reverse=True)
    err = max(abs(mp.mpf(a) - b) for a, b in zip(sp, cl)) if len(sp) == 5 else 1
    S = bits(cl)
    ok = err < 1e-12 and (S > 2 if e <= mp.mpf(1) / 256 else True)
    print(f"e={mp.nstr(e,6)} rank={len(sp)} max|numeric-closed|={float(err):.2e} entropy_bits={mp.nstr(S,25)} {'ok' if ok else 'BAD'}")
    if not ok:
        bad += 1
    if e <= mp.mpf(1) / 256 and not (S - 2 > e ** 2 / mp.log(2)):
        bad += 1; print('margin below e^2/ln2', e)

# B: rational coefficient matrix at e = 1/257
e = mp.mpf(1) / 257
C = mp.matrix([[e, mp.sqrt(e * (1 - e))], [mp.sqrt(e * (1 - e)), (1 - e) / 2]])
target = mp.matrix([[1, 16], [16, 128]]) / 257
if mp.norm(C - target) > mp.mpf(10) ** -50 or abs((1 - e) / 2 - mp.mpf(128) / 257) > mp.mpf(10) ** -50:
    bad += 1; print('B rational matrix mismatch')
else:
    print('B: C_e at e=1/257 equals (1/257)[[1,16],[16,128]]; scalar blocks 128/257')

# C: more controls, split across the cut, each control with probability (1-e)^(1/r) of being one
for (l, r) in [(1, 2), (2, 2)]:
    e = 1 / 257
    pl = [(1 - e) ** (1 / l)] * l
    pr = [(1 - e) ** (1 / r)] * r
    v = state_vector(l, r, pl, pr)
    w = apply_cswap(v, l + r)
    sp = cut_spectrum(w, l + r, l)
    cl = sorted(closed(e), reverse=True)
    err = max(abs(mp.mpf(a) - b) for a, b in zip(sp, cl)) if len(sp) == 5 else 1
    print(f"C: n={l + r + 2} controls split {l}+{r}: rank={len(sp)} max|numeric-closed|={float(err):.2e}")
    if err > 1e-10:
        bad += 1
print('bad=', bad)
sys.exit(0 if bad == 0 else 1)
