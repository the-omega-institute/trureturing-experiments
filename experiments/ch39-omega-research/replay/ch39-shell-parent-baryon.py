"""Order, light charges, axion character, and neutral mass directions of baryon triality."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd


# Entries (turns of the gauge element, internal line degree).
phase = F(1,12)
color_lepton = [F(1,36)]*3 + [-phase]  # exp(i*pi*(B-L)/6) on the four
E = [((a+phase)%1,1) for a in color_lepton for i in range(3)]
C = [(F(0),-2)]*3
Cstar = [(F(0),2)]*3
six = [((a+b)%1,0) for a,b in combinations(color_lepton,2)]
shared = six + [(F(0),0)]*8
large = E+C+C+shared
small = C+Cstar+shared
assert (len(large),len(small)) == (32,20)
for branch in (large,small):
    assert sum(q for q,d in branch)%1 == 0 and sum(d for q,d in branch) == 0


def order(entries):
    n = 1
    for q,_ in entries:
        n = n*q.denominator//gcd(n,q.denominator)
    return n


assert order(large+small)==18 and order(E)==9
assert [q for q,_ in E] == [F(1,9)]*9 + [F(0)]*3
assert [q for q,_ in six] == [F(1,18),F(1,18),F(17,18),F(1,18),F(17,18),F(17,18)]
axion = sum(q*d for q,d in large)-sum(q*d for q,d in small)
assert axion == 1
# The two possible spectator masses have phase zero, even before the SU4 part is applied.
assert large[12:15] == small[:3] and large[18:] == small[6:]
assert [((-q)%1,-d) for q,d in small[3:6]] == large[15:18]
# Three family left doublets: nine quark doublets of rotation 1/9 and three neutral leptons.
assert sum(q for q,d in E) == 1
assert 3*F(1,9)%1 == F(1,3) and 6*F(1,9)%1 == F(2,3)
assert (F(-1,2)+F(1,2))==0 and E[-1]==(0,1)
print("order in parent",order(large+small),"order on physical quarks",order(E))
print("integrated tensor character, in turns",axion,"weak instanton phase in turns",sum(q for q,d in E))
print("charges of QQQL, squared, cubed in turns",*[k*F(1,3)%1 for k in (1,2,3)])
