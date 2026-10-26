"""Liang-Sagan Question 8.6(a): is the middle order restricted to Pin_n(sigma) a distributive lattice?

Definitions follow arXiv:2408.02782v2, Section 5 and Question 8.6:
  pword(pi)  = values pi_j at internal peaks pi_{j-1} < pi_j > pi_{j+1}, left to right;
  iota_v(pi) = number of w < v that appear to the right of v in pi;
  pi <= tau  iff iota(pi) <= iota(tau) componentwise (the middle order).

Usage: python3 ls-pinnacle.py [NMAX]   (default NMAX = 6; NMAX = 7 takes about 10 minutes)
Checks:
  1. the family A_n, B_n (n = 6..9) lies in Pin_n((n-2, n)) and has no common upper bound there;
  2. for n <= NMAX, lists the sigma for which (Pin_n(sigma), <=) is not a lattice.
Exits nonzero if any expected reading fails.
"""
import itertools, collections, sys

def pword(p):
    return tuple(p[j] for j in range(1, len(p) - 1) if p[j - 1] < p[j] > p[j + 1])

def iota(p):
    pos = {v: i for i, v in enumerate(p)}
    return tuple(sum(1 for w in range(1, v) if pos[w] > pos[v]) for v in range(1, len(p) + 1))

def le(a, b):
    return all(x <= y for x, y in zip(a, b))

def family(n):
    A = (1, n - 2, 2, n) + tuple(range(3, n - 2)) + (n - 1,)
    B = (n - 1, 1, n - 2, 2, n) + tuple(range(3, n - 2))
    return A, B

def is_lattice(I):
    for a in I:
        for b in I:
            ubs = [c for c in I if le(a, c) and le(b, c)]
            if sum(1 for c in ubs if not any(d != c and le(d, c) for d in ubs)) != 1:
                return False
            lbs = [c for c in I if le(c, a) and le(c, b)]
            if sum(1 for c in lbs if not any(d != c and le(c, d) for d in lbs)) != 1:
                return False
    return True

ok = True
assert iota((4, 1, 5, 6, 3, 2)) == (0, 0, 1, 3, 2, 2)  # paper's worked example 415632
for n in range(6, 10):
    A, B = family(n)
    ia, ib = iota(A), iota(B)
    good = pword(A) == pword(B) == (n - 2, n)
    ubs = sum(1 for p in itertools.permutations(range(1, n + 1))
              if pword(p) == (n - 2, n) and le(ia, iota(p)) and le(ib, iota(p)))
    print(f"n={n} A={A} B={B} iota_n(A)={ia[n-1]} iota_(n-1)(B)={ib[n-2]} common_upper_bounds={ubs}")
    ok &= good and ubs == 0
nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
for n in range(1, nmax + 1):
    g = collections.defaultdict(list)
    for p in itertools.permutations(range(1, n + 1)):
        g[pword(p)].append(iota(p))
    bad = sorted(s for s, I in g.items() if not is_lattice(I))
    print(f"n={n} classes={len(g)} nonlattice={len(bad)} {bad}")
    if n <= 5:
        ok &= not bad
sys.exit(0 if ok else 1)
