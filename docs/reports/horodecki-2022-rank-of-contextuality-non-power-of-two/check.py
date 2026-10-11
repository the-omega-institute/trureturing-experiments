#!/usr/bin/env python3
"""Horodecki, Zhou, Stankiewicz, Salazar, Horodecki, Raussendorf, Horodecki, Ramanathan, Tyhurst, "The rank of
contextuality", arXiv:2205.10307v1; New J. Phys. 25 (2023) 073003. Section IX asks for "a quantum behavior having the
rank of contextuality different from the power of 2". Definition III.1: RC(P) is the least size of a set S of
noncontextual behaviors such that every context X has some N in S with N(.|X) = P(.|X).

Exact checks (integer and Fraction arithmetic only):
  1. Peres-Mermin square P on C^4, state I/4: Born probabilities (1 + delta s t u)/8 in every row and column.
  2. P is contextual: every one of the 2^9 global assignments gives positive weight to some outcome of zero P-probability.
  3. For each of the six contexts c a noncontextual model N_c (uniform over 16 assignments) equals P on the five other
     contexts and differs from P on c.  {N_R1, N_R2} covers P, so RC(P) = 2.
  4. Q = P (x) P on C^16, state I/16, 18 observables A (x) I and I (x) B, 36 contexts: Born probabilities equal the
     products P(.|S) P(.|T) (direct 16x16 computation).
  5. M_j = N_Rj (x) N_Rj, j = 1, 2, 3, cover all 36 contexts of Q, so RC(Q) <= 3.  The lower bound RC(Q) >= 3 is the
     marginal argument (each model's copy-1 and copy-2 marginals are noncontextual, step 2), which needs no enumeration.
  6. For k = 2..5 the constant products N_c^{(x)k}, c ranging over k+1 contexts, cover every one of the 6^k contexts of
     P^{(x)k}: RC(P^{(x)k}) <= k + 1, matching the marginal lower bound k + 1."""
import sys, itertools
from fractions import Fraction as F
import numpy as np

bad = 0
def report(label, ok):
    global bad
    print(label, ok)
    bad += 0 if ok else 1

I2 = np.eye(2, dtype=complex); X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]]); Z = np.array([[1, 0], [0, -1]], dtype=complex)
k2 = np.kron
grid = [[k2(Z, I2), k2(I2, Z), k2(Z, Z)], [k2(I2, X), k2(X, I2), k2(X, X)], [k2(Z, X), k2(X, Z), k2(Y, Y)]]
obs = [grid[r][c] for r in range(3) for c in range(3)]          # observable v = 3r + c
contexts = [[3 * r + c for c in range(3)] for r in range(3)] + [[3 * r + c for r in range(3)] for c in range(3)]
names = ['R1', 'R2', 'R3', 'C1', 'C2', 'C3']
target = [0, 0, 0, 0, 0, 1]                                     # parity of outcome bits (eigenvalue (-1)^b)

def born(ops, state_dim, bits):
    """tr(rho * prod (I + (-1)^b O)) / 2^len with rho = I / state_dim; exact because entries are Gaussian integers."""
    M = np.eye(state_dim, dtype=complex)
    for O, b in zip(ops, bits):
        M = M @ (np.eye(state_dim) + (-1) ** b * O)
    t = np.trace(M)
    num = int(round(t.real)); assert abs(t.real - num) < 1e-9 and abs(t.imag) < 1e-9
    return F(num, state_dim * 2 ** len(ops))

# 1. Born probabilities and commutation
comm = all(np.allclose(obs[u] @ obs[v], obs[v] @ obs[u]) for ctx in contexts for u in ctx for v in ctx)
herm = all(np.allclose(O, O.conj().T) and np.allclose(O @ O, np.eye(4)) for O in obs)
P = {}
ok = comm and herm
for x, ctx in enumerate(contexts):
    for bits in itertools.product((0, 1), repeat=3):
        p = born([obs[v] for v in ctx], 4, bits)
        P[x, bits] = p
        ok &= p == (F(1, 4) if sum(bits) % 2 == target[x] else 0)
report('1. PM Born probabilities are 1/4 on parity-consistent outcomes, 0 otherwise; commuting +-1 observables:', ok)

# 2. contextuality
ok = all(any(P[x, tuple(g[v] for v in ctx)] == 0 for x, ctx in enumerate(contexts))
         for g in itertools.product((0, 1), repeat=9))
report('2. every global assignment of the 9 observables hits a zero-probability outcome (P contextual):', ok)

# 3. single-context-defect models
assignments = list(itertools.product((0, 1), repeat=9))
def parities(g):
    return [sum(g[v] for v in ctx) % 2 for ctx in contexts]
def model(defect):
    want = [t ^ (1 if x == defect else 0) for x, t in enumerate(target)]
    supp = [g for g in assignments if parities(g) == want]
    return {g: F(1, len(supp)) for g in supp}, len(supp)
def marg(mu, x):
    out = {}
    for g, w in mu.items():
        key = tuple(g[v] for v in contexts[x]); out[key] = out.get(key, 0) + w
    return {bits: out.get(bits, F(0)) for bits in itertools.product((0, 1), repeat=3)}
Nmods = {}
ok = True
for c in range(6):
    mu, n = model(c); Nmods[c] = mu
    agree = [all(marg(mu, x)[b] == P[x, b] for b in itertools.product((0, 1), repeat=3)) for x in range(6)]
    ok &= n == 16 and agree == [x != c for x in range(6)]
    print(f'   N_{names[c]}: support {n}, agrees with P exactly on', [names[x] for x in range(6) if agree[x]])
report('3. each N_c is a probability distribution on 16 assignments agreeing with P except on context c:', ok)
report('3. {N_R1, N_R2} covers P (RC(P) = 2):',
       all(any(all(marg(Nmods[c], x)[b] == P[x, b] for b in itertools.product((0, 1), repeat=3)) for c in (0, 1))
           for x in range(6)))

# 4. Q = P (x) P from 16x16 Born probabilities
I4 = np.eye(4)
obsQ = [k2(O, I4) for O in obs] + [k2(I4, O) for O in obs]
ok = True
for S, T in itertools.product(range(6), repeat=2):
    ops = [obsQ[v] for v in contexts[S]] + [obsQ[9 + v] for v in contexts[T]]
    for bits in itertools.product((0, 1), repeat=6):
        ok &= born(ops, 16, bits) == P[S, bits[:3]] * P[T, bits[3:]]
ok &= all(np.allclose(obsQ[u] @ obsQ[v], obsQ[v] @ obsQ[u])
          for S, T in itertools.product(range(6), repeat=2)
          for u in contexts[S] + [9 + w for w in contexts[T]] for v in contexts[S] + [9 + w for w in contexts[T]])
report('4. Q(.|S,T) = P(.|S) P(.|T) for all 36 contexts x 64 outcomes (16x16 Born rule):', ok)

# 5. M_j = N_Rj (x) N_Rj covers Q; product marginals agree iff both factors agree
def prod_agrees(c1, c2, S, T):
    m1, m2 = marg(Nmods[c1], S), marg(Nmods[c2], T)
    return all(m1[a] * m2[b] == P[S, a] * P[T, b] for a in itertools.product((0, 1), repeat=3)
               for b in itertools.product((0, 1), repeat=3))
cover = {(S, T): [j for j in range(3) if prod_agrees(j, j, S, T)] for S, T in itertools.product(range(6), repeat=2)}
ok = all(cover[st] for st in cover)
uncovered_by_two = {pair: [st for st in cover if not set(cover[st]) & set(pair)] for pair in itertools.combinations(range(3), 2)}
print('   contexts missed by each pair of the three models:',
      {f'M{a+1},M{b+1}': [(names[s], names[t]) for s, t in v] for (a, b), v in uncovered_by_two.items()})
report('5. {M_1, M_2, M_3} covers all 36 contexts of Q (RC(Q) <= 3):', ok)

# 6. constant-word covers of P^{(x)k}, k = 2..5, via the per-context agreement table of step 3
agreeP = {(c, x): x != c for c in range(6) for x in range(6)}
for k in range(2, 6):
    ok = all(any(all(agreeP[c, x] for x in word) for c in range(k + 1))
             for word in itertools.product(range(6), repeat=k))
    report(f'6. k={k}: the {k+1} models N_c^(x){k} cover all {6**k} contexts of P^(x){k} (RC <= {k+1}):', ok)
print('bad=', bad)
sys.exit(1 if bad else 0)
