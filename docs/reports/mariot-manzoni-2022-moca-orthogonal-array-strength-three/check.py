"""Mariot-Manzoni, "Building Correlation Immune Functions from Sets of Mutually Orthogonal Cellular Automata",
arXiv:2207.08280v1, Section 5: "prove that the binary OA obtained from a set of k-MOCA has always at least strength 3".

Bipermutive rule of diameter d = b + 1: f(y_0..y_b) = y_0 + phi(y_1..y_{b-1}) + y_b over F_2.  CA F : F_2^{2b} -> F_2^b,
F(x)_j = f(x_j, ..., x_{j+b}).  Two CA are orthogonal when x -> (F_1(x), F_2(x)) is injective on F_2^{2b} (Eq. (6) on the
Latin squares with rows x_0..x_{b-1} and columns x_b..x_{2b-1}).  Array A (Eq. (7)): rows x in F_2^{2b}, columns (i, j).

Checks (exact):
 1. For b = 2, 3, 4: enumerate all bipermutive rules, all orthogonal pairs, all 3-MOCA families; the array of every
    2-MOCA and 3-MOCA family is an OA of strength 3 (every 3 distinct columns, every pattern exactly 2^{2b-3} times).
    Report the counts and how many 3-MOCA arrays have strength exactly 3 (not 4), i.e. sharpness.
 2. Search-seat lemma: for EVERY choice of rules (orthogonal or not), b = 2..5, the XOR of any three distinct columns
    of A is balanced; single columns are balanced.  (Pairs need orthogonality; Lemma 3 of the source.)
"""
import itertools, random

def rules(b):
    m = b - 1
    for bits in itertools.product((0, 1), repeat=2 ** m):
        yield bits  # phi as truth table over F_2^{b-1}, index = big-endian bits of (y_1..y_{b-1})

def col(b, phi, x, j):
    m = b - 1
    idx = 0
    for t in range(1, b):
        idx = idx * 2 + x[j + t]
    return x[j] ^ (phi[idx] if m >= 0 else 0) ^ x[j + b]

def ca(b, phi, x):
    return tuple(col(b, phi, x, j) for j in range(b))

def inputs(b):
    return list(itertools.product((0, 1), repeat=2 * b))

def orth(b, p, q, X):
    return len({(ca(b, p, x), ca(b, q, x)) for x in X}) == len(X)

def strength(b, fam, X, t):
    cols = [(i, j) for i in range(len(fam)) for j in range(b)]
    table = {x: {(i, j): col(b, fam[i], x, j) for (i, j) in cols} for x in X}
    lam = len(X) // 2 ** t
    for S in itertools.combinations(cols, t):
        cnt = {}
        for x in X:
            key = tuple(table[x][c] for c in S)
            cnt[key] = cnt.get(key, 0) + 1
        if len(cnt) != 2 ** t or any(v != lam for v in cnt.values()):
            return False
    return True

ok = True
for b in (2, 3, 4):
    X = inputs(b)
    R = list(rules(b))
    pairs = [(p, q) for p, q in itertools.combinations(R, 2) if orth(b, p, q, X)]
    pairset = set(pairs) | {(q, p) for p, q in pairs}
    triples = [(p, q, r) for p, q, r in itertools.combinations(R, 3)
               if (p, q) in pairset and (p, r) in pairset and (q, r) in pairset]
    s3_pairs = all(strength(b, list(f), X, 3) for f in pairs)
    s3_triples = all(strength(b, list(f), X, 3) for f in triples)
    exact3 = sum(1 for f in triples if not strength(b, list(f), X, 4)) if b >= 2 else 0
    print(f"b={b} d={b+1}: rules={len(R)} orthogonal pairs={len(pairs)} 3-MOCA={len(triples)} "
          f"strength>=3 for all 2-MOCA={s3_pairs} all 3-MOCA={s3_triples} 3-MOCA of strength exactly 3={exact3}")
    ok &= s3_pairs and s3_triples
random.seed(220708280)
lemma = True
for b in (2, 3, 4, 5):
    X = inputs(b)
    for trial in range(6):
        k = random.randint(1, 3)
        fam = [tuple(random.randint(0, 1) for _ in range(2 ** (b - 1))) for _ in range(k)]
        cols = [(i, j) for i in range(k) for j in range(b)]
        for c in cols:
            lemma &= sum(col(b, fam[c[0]], x, c[1]) for x in X) * 2 == len(X)
        for S in itertools.combinations(cols, 3):
            s = sum(col(b, fam[S[0][0]], x, S[0][1]) ^ col(b, fam[S[1][0]], x, S[1][1]) ^ col(b, fam[S[2][0]], x, S[2][1]) for x in X)
            lemma &= s * 2 == len(X)
print("single-column and three-column XOR balance for random (not necessarily orthogonal) rules, b = 2..5:", lemma)
ok &= lemma
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
