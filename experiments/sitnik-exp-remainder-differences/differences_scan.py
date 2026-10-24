"""Forward differences in n of g_n(x) = R_{n-1}(x)/R_n(x), R_n(x) = e^x - sum_{k<=n} x^k/k!.
For each order k = 1..12 report the first (x, n) on the grid with a negative difference,
then the witness k = 3, n = 1, x = 1 and the k = 3 rows for x in {0.001, 1, 5}."""
import mpmath as mp

mp.mp.dps = 60


def R(n, x):
    return x ** (n + 1) / mp.factorial(n + 1) * mp.hyp1f1(1, n + 2, x)


def g(n, x):
    return R(n - 1, x) / R(n, x)


def D(k, n, x):
    return mp.fsum((-1) ** (k - j) * mp.binomial(k, j) * g(n + j, x) for j in range(k + 1))


grid = [0.001, 0.5, 1, 2, 3, 5, 8, 12, 20, 30, 50]
for k in range(1, 13):
    found = None
    for x in grid:
        for n in range(1, 40):
            v = D(k, n, mp.mpf(x))
            if v < 0:
                found = (x, n, mp.nstr(v, 5))
                break
        if found:
            break
    print(f"k={k} first negative (x, n, value) = {found}")
e = mp.e
s = [1, 2, mp.mpf(5) / 2, mp.mpf(8) / 3, mp.mpf(65) / 24]
G = [(e - s[m - 1]) / (e - s[m]) for m in range(1, 5)]
print("witness k=3 n=1 x=1:", mp.nstr(G[3] - 3 * G[2] + 3 * G[1] - G[0], 12), "direct:", mp.nstr(D(3, 1, mp.mpf(1)), 12))
for x in [0.001, 1, 5]:
    print(f"k=3 x={x}:", [mp.nstr(D(3, n, mp.mpf(x)), 5) for n in range(1, 9)])
