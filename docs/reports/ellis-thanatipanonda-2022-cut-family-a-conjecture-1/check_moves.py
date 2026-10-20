"""Check the explicit moves of the search seat's reachability lemma for Ellis-Thanatipanonda Family A.

For b = 2c and every n <= N, F(n) = 2*floor((n-1)/b) + ((n-1) mod 2). Each listed move must use a cut d in {1, 3, b}
(d + 1 positive parts summing to n), and the xor of F over its parts must cover every value below F(n).
"""
import sys

N = 400


def F(n, b):
    return 2 * ((n - 1) // b) + ((n - 1) % 2)


def moves(n, b):
    c = b // 2
    q, r = divmod(n - 1, b)
    r += 1
    eps = (r - 1) % 2
    out = []
    for x in range(q):
        out.append([b * x + r] + [q - x] * b)
    if r > 1:
        for x in range(q):
            out.append([b * x + r - 1, 1, c * (q - x), c * (q - x)])
        if eps == 1:
            out.append([n - 1, 1])
    else:
        for x in range(q - 1):
            out.append([b * (x + 1), 1, c * (q - x - 1), c * (q - x - 1)])
        if q >= 1:
            out.append([b * q, 1])
    return out


ok = True
for b in (4, 6, 8, 10, 12):
    for n in range(1, N + 1):
        got = set()
        for m in moves(n, b):
            if sum(m) != n or min(m) < 1 or len(m) - 1 not in (1, 3, b):
                print(f"bad move b={b} n={n} {m}")
                ok = False
            v = 0
            for h in m:
                v ^= F(h, b)
            got.add(v)
        if not set(range(F(n, b))) <= got or F(n, b) in got:
            print(f"coverage fails b={b} n={n}")
            ok = False
    print(f"b={b}: moves valid and cover [0, F(n)) without F(n) for n <= {N}: {ok}")
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
