"""Checks for inequality (4) of Gupta, arXiv:2606.22614v2, Section 4.

disp(phi) = min_K sum_v |phi(v) - K| over integers K; D(m, n) = disp of (i, j) -> i + j on the 1-indexed grid.
(a) D(2, n) = ceil(n^2 / 2) (Proposition 4.3) and D(m, m) = (m^3 - m) / 3 (the square diagonal after Theorem 4.1).
(b) Exhaustive: every integer 1-Lipschitz function on G_{m,n} with phi(1,1) = 0 has disp <= D(m, n), for the grids the
    source enumerates (up to 4 x 4 and 3 x 5) and the two-row grids up to 2 x 7.
(c) Monotone functions (Proposition 4.4 reduces (4) to them) on larger grids: 3 x 6, 4 x 5, 5 x 5.
"""
import sys


def disp(values):
    v = sorted(values)
    med = v[len(v) // 2]
    return sum(abs(x - med) for x in v)


def D(m, n):
    return disp([i + j for i in range(1, m + 1) for j in range(1, n + 1)])


def lipschitz_functions(m, n):
    cells = [(i, j) for i in range(m) for j in range(n)]
    phi = {}

    def rec(k):
        if k == len(cells):
            yield list(phi.values())
            return
        i, j = cells[k]
        if k == 0:
            candidates = [0]
        else:
            lo, hi = -10 ** 9, 10 ** 9
            for (a, b) in ((i - 1, j), (i, j - 1)):
                if (a, b) in phi:
                    lo = max(lo, phi[(a, b)] - 1)
                    hi = min(hi, phi[(a, b)] + 1)
            candidates = range(lo, hi + 1)
        for c in candidates:
            phi[(i, j)] = c
            yield from rec(k + 1)
            del phi[(i, j)]

    yield from rec(0)


def monotone_functions(m, n):
    cells = [(i, j) for i in range(m) for j in range(n)]
    phi = {}

    def rec(k):
        if k == len(cells):
            yield list(phi.values())
            return
        i, j = cells[k]
        if k == 0:
            candidates = [0]
        else:
            lo, hi = -10 ** 9, 10 ** 9
            for (a, b) in ((i - 1, j), (i, j - 1)):
                if (a, b) in phi:
                    lo = max(lo, phi[(a, b)])
                    hi = min(hi, phi[(a, b)] + 1)
            candidates = range(lo, hi + 1)
        for c in candidates:
            phi[(i, j)] = c
            yield from rec(k + 1)
            del phi[(i, j)]

    yield from rec(0)


def main():
    bad = 0
    for n in range(2, 12):
        ok = D(2, n) == -(-n * n // 2)
        bad += not ok
        print(f"D(2,{n})={D(2, n)} ceil(n^2/2)={-(-n * n // 2)} {'OK' if ok else 'MISMATCH'}")
    for m in range(1, 9):
        ok = D(m, m) == (m ** 3 - m) // 3
        bad += not ok
        print(f"D({m},{m})={D(m, m)} (m^3-m)/3={(m ** 3 - m) // 3} {'OK' if ok else 'MISMATCH'}")
    for (m, n) in [(2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (3, 3), (3, 4), (3, 5), (4, 4)]:
        count = 0
        worst = 0
        for vals in lipschitz_functions(m, n):
            count += 1
            worst = max(worst, disp(vals))
        ok = worst == D(m, n)
        bad += not ok
        print(f"lipschitz {m}x{n}: functions={count} max_disp={worst} D={D(m, n)} {'OK' if ok else 'MISMATCH'}")
    for (m, n) in [(3, 6), (4, 5), (5, 5)]:
        count = 0
        worst = 0
        for vals in monotone_functions(m, n):
            count += 1
            worst = max(worst, disp(vals))
        ok = worst == D(m, n)
        bad += not ok
        print(f"monotone {m}x{n}: functions={count} max_disp={worst} D={D(m, n)} {'OK' if ok else 'MISMATCH'}")
    print(f"EXIT={1 if bad else 0}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
