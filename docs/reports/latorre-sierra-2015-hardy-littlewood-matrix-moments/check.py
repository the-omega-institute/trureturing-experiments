#!/usr/bin/env python3
"""Checks for Latorre-Sierra, arXiv:1403.4765v2 / QIC 15 (2015), conjecture (68):
2^(-ms) Tr(C_m^s) -> prod_{p>2} (1 + (p-1)^(1-2s)),  s >= 2,
C_m: d = 2^(m-1), (C_m)_ii = 0, (C_m)_ij = C(2|i-j|), C(2h) = 2 prod_{p>2}(1 - 1/(p-1)^2) prod_{p>2, p|h} (p-1)/(p-2).

1. Local-factor identity (exact rationals): for each odd prime p, a_p^s [ (p + 1/(p-2))^s + (p-1)(1/(p-2))^s ] / p^s
   = 1 + (p-1)^(1-2s), with a_p = p(p-2)/(p-1)^2; this is the per-prime factor of the fixed-cutoff limit
   (K_p = J + (r_p - 1) I on Z/pZ).
2. Numerics (float64): the ratio 2^(-ms) Tr(C_m^s) / prod approaches 1 monotonically for s = 2..5, m = 6..12.
Prints ALL_OK on success.
"""
from fractions import Fraction
import numpy as np

def primes(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(3, n + 1) if s[i]]

ok = True
for p in primes(200):
    for s in range(2, 7):
        ap = Fraction(p * (p - 2), (p - 1) ** 2); e = Fraction(1, p - 2)
        lhs = ap ** s * ((p + e) ** s + (p - 1) * e ** s) / Fraction(p) ** s
        if lhs != 1 + Fraction(1, (p - 1) ** (2 * s - 1)):
            print('LOCAL_FAIL', p, s); ok = False
print('local-factor identity ok: odd p < 200, s = 2..6')

P = primes(2_000_000)
logA = np.log(2.0) + float(np.sum(np.log1p(-1.0 / (np.array(P, dtype=float) - 1) ** 2)))
A = np.exp(logA)
def target(s):
    pa = np.array(P, dtype=float)
    return float(np.exp(np.sum(np.log1p((pa - 1) ** (1 - 2 * s)))))
def Cvals(dmax):
    F = np.ones(dmax)
    for p in P:
        if p >= dmax: break
        F[p::p] *= (p - 1) / (p - 2)
    return A * F   # C(2h) = A F(h), h >= 1
ratios = {s: [] for s in range(2, 6)}
for m in range(6, 13):
    d = 2 ** (m - 1); Cv = Cvals(d)
    idx = np.abs(np.subtract.outer(np.arange(d), np.arange(d)))
    C = np.where(idx == 0, 0.0, Cv[idx])
    ev = np.linalg.eigvalsh(C)
    for s in range(2, 6):
        ratios[s].append(float(np.sum(ev ** s)) / 2.0 ** (m * s) / target(s))
    print('m=%2d ' % m + ' '.join('s=%d:%.6f' % (s, ratios[s][-1]) for s in range(2, 6)), flush=True)
for s in range(2, 6):
    r = ratios[s]
    if not (abs(r[-1] - 1) < abs(r[0] - 1) and abs(r[-1] - 1) < 0.05):
        print('NOT_APPROACHING', s, r); ok = False
print('ALL_OK' if ok else 'FAIL')
