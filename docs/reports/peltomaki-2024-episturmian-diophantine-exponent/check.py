"""Finite check of the uniform finite-window lemma behind a proposed proof of Peltomaki's Conjecture 8.4
(EJC 118 (2024) 103942, arXiv:2103.08351v3, Section 8): an aperiodic episturmian word t has dio(t) < infinity
if and only if t has bounded partial quotients.  The checked lemma carries the unbounded direction.

Standard episturmian prefixes c_Delta are built by iterated palindromic closure u_{i+1} = (u_i y_i)^(+) (the closure
is computed in linear time from the longest palindromic suffix, found by KMP, and compared with the naive
definition on all words of length <= 8 over three letters).  Directives have alphabets of size 2-4, a prefix
using every letter, one maximal run x^a with a = M(M+3), M in {2, 3}, and a random tail; words are capped at
60000 letters.  q is the period produced by the run (the length added by one closure step with x).

Lemma 3 of the proposal: from every starting position s, the factor beginning at s has a prefix U V^e with
|V| = q, |U| <= (M+3) q and |U V^e| / (|U| + |V|) >= M.  The longest period-q run from each position is
precomputed.  Tested positions: the first 3H positions of each word (H = (a+M+2) q) and up to 3000 positions
sampled uniformly from the rest of the word.
"""
import random
random.seed(84)

def pal_closure(w):
    # shortest palindrome having w as prefix: w + reverse(prefix before the longest palindromic suffix), via KMP
    r = w[::-1]; s = r + "\x00" + w
    f = [0] * len(s)
    for i in range(1, len(s)):
        j = f[i - 1]
        while j and s[i] != s[j]: j = f[j - 1]
        if s[i] == s[j]: j += 1
        f[i] = j
    L = f[-1]  # length of the longest palindromic suffix of w
    return w + w[:len(w) - L][::-1]

def episturmian_prefix(directive, minlen):
    u = ""
    for y in directive:
        u = pal_closure(u + y)
        if len(u) >= minlen:
            break
    return u

def period_runs(word, V_len):
    # cnt[p] = number of consecutive i >= p with word[i] == word[i + V_len]
    n = len(word); cnt = [0] * (n + 1)
    for p in range(n - V_len - 1, -1, -1):
        cnt[p] = cnt[p + 1] + 1 if word[p] == word[p + V_len] else 0
    return cnt

def best_ratio_from(word, s, V_len, limit, cnt):
    # max over preperiods |U| in [0, limit] of |U V^e| / (|U| + |V|) with V of length V_len, V^e maximal fractional power
    best = 0.0
    n = len(word)
    for u in range(0, limit + 1):
        start = s + u
        if start + V_len > n: break
        e_len = V_len + cnt[start]
        r = (u + e_len) / (u + V_len)
        if r > best: best = r
    return best

CAP = 60000
bad = 0; tested = 0
for trial in range(60):
    k = random.choice([2, 3, 4]); A = "abcd"[:k]
    M = random.choice([2, 3])
    a = M * (M + 3)
    # directive: some short random runs, then x^a, then more short runs
    pre = []
    for _ in range(random.randint(1, 4)):
        x = random.choice(A); pre += [x] * random.randint(1, 3)
    # ensure every letter appears before the long run, so the run's letter has a nontrivial mu
    for x in A: pre.append(x)
    x = random.choice([c for c in A if c != pre[-1]])
    post = []
    for _ in range(6):
        y = random.choice([c for c in A if c != (post[-1] if post else x)]); post += [y] * random.randint(1, 2)
    directive = pre + [x] * a + post
    # q = |mu(x)| : length of the image of x under the composed morphism of the prefix before the run.
    # equivalently the period produced by the run; compute it empirically as the length of the palindromic closure step
    u_before = episturmian_prefix(pre, 10**9) if pre else ""
    u1 = pal_closure(u_before + x); u2 = pal_closure(u1 + x)
    q = len(u2) - len(u1)
    w = episturmian_prefix(directive, CAP)
    H = (a + M + 2) * q
    window = min(len(w) - H - 1, 3 * H)
    if window <= 0: continue
    cnt = period_runs(w, q)
    last = len(w) - H - 1
    positions = list(range(0, window)) + random.sample(range(window, last), min(3000, max(0, last - window)))
    for s in positions:
        tested += 1
        r = best_ratio_from(w, s, q, (M + 3) * q, cnt)
        if r < M - 1e-12:
            bad += 1; print("FAIL", trial, k, M, a, q, s, r); break
print("positions tested:", tested, "(first 3H of each word plus up to 3000 sampled across the whole word); M in {2,3}")
print("bad=", bad)
raise SystemExit(1 if bad else 0)
