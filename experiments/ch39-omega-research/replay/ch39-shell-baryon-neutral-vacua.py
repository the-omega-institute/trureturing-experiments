"""Lowest tensor degree of a color/electrically neutral direction breaking the parent element."""
from fractions import Fraction as F
from itertools import combinations_with_replacement, product


# A color irrep is encoded by 0 (trivial), 1 (3), -1 (3-bar).
# Each weight is (color, electric charge, turns of g_B).
def singlet(q):
    return (0,F(q),F(0))


def dual(block):
    return {(-c,-q,(-g)%1) for c,q,g in block}


weak = {singlet(F(1,2)),singlet(F(-1,2))}
adj_six = {(1,F(-1,3),F(17,18)),(-1,F(1,3),F(1,18))}
T = adj_six | weak | {singlet(0),singlet(1),singlet(-1)}
D = T | {singlet(0)}
E = {(1,F(1,6),F(1,9)),singlet(F(-1,2))}
big = E | D
letters = (big,dual(big),D,dual(D),weak)


def dangerous_pair(a,b):
    return any(c+h==0 and q+v==0 and (g+k)%1
               for (c,q,g),(h,v,k) in product(a,b))


assert all(not any(c==0 and q==0 and g for c,q,g in a) for a in letters)
assert not any(dangerous_pair(a,b) for a,b in combinations_with_replacement(letters,2))
new_vev = [(1,F(1,6),F(1,9)),(-1,F(1,3),F(1,18)),singlet(F(-1,2))]
assert all(v in a for v,a in zip(new_vev,(big,D,weak)))
assert sum(v[1] for v in new_vev)==0 and sum(v[2] for v in new_vev)==F(1,6)
assert [n for n in range(18) if (F(n,6)%1)==0] == [0,6,12]
# On the six-block, compare the sixth power with the square of the color center.
for c,q,g in big|D:
    assert (6*g-F(2*c,3))%1==0
print("no color- and charge-neutral breaking weight at tensor degree one or two")
print("F x D x right-doublet color contraction: Y=0, g_B charge 1/6 turn")
print("preserved powers of g_B are 0,6,12 modulo 18; sixth power equals squared color center")
