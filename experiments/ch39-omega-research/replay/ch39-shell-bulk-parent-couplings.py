"""A common-coupling bulk parent and the face-energy and perturbativity checks."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path


old = run_path(str(Path(__file__).with_name("ch39-shell-bulk-parent-spectrum.py")))
energies = old["energies"]

# (dimension, central charge, SU4, family SU3, left SU2, right SU2).
E = (12, 1, 3, 4, 0, 0)
C = (3, -2, 0, 1, 0, 0)        # B tensored with inverse J
Cstar = (3, 2, 0, 1, 0, 0)
six = (6, 0, 2, 0, 0, 0)        # exterior square of A tensored with J
lr = (4, 0, 0, 0, 2, 2)
left, right = (2, 0, 0, 0, 1, 0), (2, 0, 0, 0, 0, 1)
shared = (six, lr, left, right)
branches = ((E, C, C) + shared, (C, Cstar) + shared)
levels = [tuple(map(sum, zip(*(block[2:] for block in branch)))) for branch in branches]
central = [sum(Q(d*q*q, 144) for d,q,*_ in branch) for branch in branches]
sizes = [sum(b[0] for b in branch) for branch in branches]
assert sizes == [32,20] and [sum(b[0]*b[1] for b in branch) for branch in branches] == [0,0]
assert levels == [(5,6,3,3), (2,2,3,3)] and central == [Q(1,4),Q(1,6)]
assert tuple(a+b+c for a,b,c in zip(*levels, (0,0,1,1))) == (7,8,7,7)
assert tuple(a-b for a,b in zip(*levels)) == (3,4,0,0)
for cap, selected in ((True,-3),(False,-1)):
    for i in (0,1):
        w4,w3,_,_ = levels[i]
        values = energies(cap,w4,w3,central[i])
        winners = [k for k,e in values.items() if e == min(values.values())]
        expected = [-3,1] if i == 1 and cap else [selected]
        assert winners == expected

budget = sum(size*row[0] for size,row in zip(sizes,levels))
equal_weight_cost = max(sizes)*7
assert (budget,equal_weight_cost) == (200,224)
assert 224/(4*Q(3141,1000)*Q(406,10)) < Q(44,100)
assert 911/(4*Q(22,7)*41) > Q(176,100)
print("gauge group SU(32) x SU(20) x SU(2)L x SU(2)R, common-coupling levels 7,8,7,7")
print("bare anomaly factor Y = c2(20)-c2(32)-p-l-r, integral on the parent")
print("net three at every strictly positive pair of internal weights")
print("adjoint numerator at equal weights",equal_weight_cost,"optimal numerator",budget)
