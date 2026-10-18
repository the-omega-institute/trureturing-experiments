"""For every p coprime to 506: a word f_p over {A, B} with f_p(ab) = christoffel(p, 506), found by peeling
A^{-1} (ab -> b, a -> a) or B^{-1} (ab -> a, b -> b). By the factor identities L(c) = 3 + |f_p|; the minimum is 17.
Writes christoffel_certificates_506.tsv (p, f_p outer-to-inner, |f_p|, quotient sum)."""
import os, sys
from math import gcd
here = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(here, 'verify506.py')).read().split("def main")[0])

def inv_A(w):  # images a -> a, b -> ab
    out, i = [], 0
    while i < len(w):
        if w[i] == 'a' and i + 1 < len(w) and w[i + 1] == 'b': out.append('b'); i += 2
        elif w[i] == 'a': out.append('a'); i += 1
        else: return None
    return ''.join(out)

def inv_B(w):  # images a -> ab, b -> b
    out, i = [], 0
    while i < len(w):
        if w[i] == 'a' and i + 1 < len(w) and w[i + 1] == 'b': out.append('a'); i += 2
        elif w[i] == 'b': out.append('b'); i += 1
        else: return None
    return ''.join(out)

def peel(w):
    seq = []
    while w != 'ab':
        a, b = inv_A(w), inv_B(w)
        # exactly one inverse keeps a word containing both letters of length >= 2
        cands = [(n, x) for n, x in (('A', a), ('B', b)) if x is not None and 'a' in x and 'b' in x]
        if len(cands) != 1: return None
        n, x = cands[0]; seq.append(n); w = x
    return ''.join(seq)

A = {"a": "a", "b": "ab"}; B = {"a": "ab", "b": "b"}
def apply(seq, x='ab'):
    for c in reversed(seq): x = ''.join((A if c == 'A' else B)[ch] for ch in x)
    return x

rows, bad = [], 0
for p in range(1, 506):
    if gcd(p, 506) != 1: continue
    c = christoffel(p, 506); f = peel(c)
    ok = f is not None and apply(f) == c and len(f) + 3 == qsum(p, 506) + 1 and L_fast(c) == 3 + len(f)
    bad += not ok; rows.append((p, f or '-', len(f) if f else -1, qsum(p, 506)))
open(os.path.join(here, 'christoffel_certificates_506.tsv'), 'w').write(''.join(f"{p}\t{f}\t{l}\t{q}\n" for p, f, l, q in rows))
print(f"certificates: {len(rows)}; failures: {bad}; min |f_p| = {min(r[2] for r in rows)} (so min L = {3 + min(r[2] for r in rows)}); max |f_p| = {max(r[2] for r in rows)}")
print("RESULT", "PASS" if bad == 0 and len(rows) == 220 and min(r[2] for r in rows) == 14 else "FAIL")
sys.exit(0 if bad == 0 else 1)
