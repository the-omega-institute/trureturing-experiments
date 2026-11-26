"""Integrated weak anomalies of bulk spectators on the degree-twelve determinant background."""
from fractions import Fraction as F
from math import gcd
from functools import reduce

# q is measured on the primitive central circle (e^{-3i theta},e^{4i theta}).
# Each sector is (multiplicity from S(U4 x U3), q, 6D chirality).
original=(12,1,1)
fund4=(4,-3,1)
fund3=(3,4,-1)
def index(block):
    d,q,sgn=block
    return d*q*sgn
def anomaly(block):
    return F(index(block)*block[1],2)
bundle=[original,fund4,fund3]
assert [index(r) for r in bundle]==[12,-12,-12]
assert [anomaly(r) for r in bundle]==[6,18,-24]
assert sum(anomaly(r) for r in bundle)==0
dual=[(d,-q,-s) for d,q,s in bundle]
assert sum(anomaly(r) for r in dual)==0
assert reduce(gcd,(36,48,144))==12
# The charge-12 line alone gives index 12, and contributes 72 for a weak doublet.
assert anomaly((1,12,1))==72
print("left integrated anomalies",*[anomaly(r) for r in bundle],"right",*[anomaly(r) for r in dual])
print("new net internal multiplicities per weak side",*[abs(index(r)) for r in bundle[1:]])
