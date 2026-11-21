from fractions import Fraction
from math import prod
def divisors(fs):
    ds = [1]
    for p, a in fs:
        ds = [d * p**e for d in ds for e in range(a + 1)]
    return sorted(ds)
ds = divisors([(2, 3), (3, 3), (5, 1), (7, 1)])
h = sum((Fraction(1, e-d) for i, d in enumerate(ds)
         for e in ds[i+1:]), Fraction()) / len(ds)
assert Fraction(99275449, 10**8) < h < Fraction(99275450, 10**8)
print(len(ds), float(h))
for fs in [[(2, 10), (3, 4), (5, 3), (7, 4)],
           [(2, 10), (3, 4), (5, 4), (7, 3)]]:
    ds = divisors(fs)
    print(prod(p**a for p, a in fs), len(ds),
          sum(10**10 // (e-d) for i, d in enumerate(ds) for e in ds[i+1:]))
