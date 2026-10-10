from math import comb
from fractions import Fraction
P = {15: 1, 7: 15}
Q = {8: 15, 0: 1}
def local_order(a):
    return next(j for j in range(16)
                if sum(c * comb(k,j) * a**(k-j) for k,c in P.items() if k >= j)
                - a * sum(c * comb(k,j) * a**(k-j) for k,c in Q.items() if k >= j))
assert sum(c for k,c in Q.items()) == 16
assert sum(c*(-1)**k for k,c in Q.items()) == 16
orders = (min(P), min(15-k for k in Q), local_order(1), local_order(-1))
assert orders == (7, 7, 3, 3)
print("source [[15,1,3]]; points (0, infinity, 1, -1); orders", orders)
