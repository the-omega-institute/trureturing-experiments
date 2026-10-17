"""Vivion arXiv:2509.11172v2, Question 3: does b_w^2 = p_w together with b_w^1 != p_w force a binary word to be
1-balanced? Family w = 1^m (0^a 1)^omega, a, m >= 2: b^2 = p for all lengths checked, b^1(2) = 3 != 4 = p(2),
and the factors 00, 11 (more generally 0^h, 1^h, h <= min(a, m)) violate 1-balance.
Factors of length n are read from a prefix long enough to contain every factor (the word is eventually periodic
with period a + 1 after position m). The k-binomial class is the vector of subword-occurrence counts of all words
of length <= k. Control: the Thue-Morse word has b^2(4) < p(4) (the source's example 0110 ~2 1001).
"""
import itertools, sys

def binom(u, v):  # number of occurrences of v as a (scattered) subword of u
    dp = [1] + [0] * len(v)
    for c in u:
        for j in range(len(v), 0, -1):
            if v[j - 1] == c: dp[j] += dp[j - 1]
    return dp[len(v)]

def sig(u, k):
    return tuple(binom(u, ''.join(x)) for L in range(1, k + 1) for x in itertools.product('01', repeat=L))

def factors(prefix, n):
    return {prefix[i:i + n] for i in range(len(prefix) - n + 1)}

def main():
    fails = []
    N = 36
    for a in range(2, 6):
        for m in range(2, 6):
            prefix = '1' * m + ('0' * a + '1') * (2 * N // (a + 1) + 8)
            ok = True
            for n in range(1, N + 1):
                F = factors(prefix, n)
                p, b2 = len(F), len({sig(u, 2) for u in F})
                if b2 != p: ok = False
            F2 = factors(prefix, 2); b1_2 = len({sig(u, 1) for u in F2})
            unbal = any(abs(u.count('1') - v.count('1')) > 1 for n in range(1, 12) for u in factors(prefix, n) for v in factors(prefix, n))
            h = min(a, m)
            witness = ('0' * h in factors(prefix, h)) and ('1' * h in factors(prefix, h))
            print(f"a={a} m={m}: b2 = p for n <= {N}: {ok}; p(2) = {len(F2)}, b1(2) = {b1_2}; not 1-balanced: {unbal}; 0^{h}, 1^{h} occur: {witness}")
            if not (ok and len(F2) == 4 and b1_2 == 3 and unbal and witness): fails.append((a, m))
    # control: Thue–Morse
    tm = ''.join(str(bin(i).count('1') % 2) for i in range(4096))
    F4 = factors(tm, 4); c = (len(F4), len({sig(u, 2) for u in F4}))
    print(f"control Thue-Morse: p(4) = {c[0]}, b2(4) = {c[1]} (expected b2 < p)")
    if not c[1] < c[0]: fails.append("control")
    print("RESULT", "PASS" if not fails else f"FAIL {fails}")
    sys.exit(0 if not fails else 1)

if __name__ == "__main__":
    main()
