"""Exact discrete charges and a conservative three-family clock-breaking bound."""
from fractions import Fraction as Q
from itertools import permutations, product
from collections import Counter
from math import sqrt


def parity(seq):
    return (-1)**sum(a>b for i,a in enumerate(seq) for b in seq[i+1:])


# Use distinct Grassmann generators for each color, weak and spinor component.
# Both weak and both spinor epsilon tensors have two indices.
four_fermion = Counter()
doublets = list(permutations((0, 1)))
for colors in permutations((0,1,2)):
    for weak1, weak2, spin1, spin2 in product(doublets, repeat=4):
        a,b,c = colors
        i,j = weak1
        k,l = weak2
        m,n = spin1
        r,s = spin2
        word = ((a,i,m),(b,j,n),(c,k,r),(3,l,s))
        coefficient = parity(colors)*parity(weak1)*parity(weak2)*parity(spin1)*parity(spin2)
        four_fermion[tuple(sorted(word))] += coefficient*parity(word)
nonzero = sorted((w,v) for w,v in four_fermion.items() if v)
assert nonzero and 3*Q(1,6)-Q(1,2) == 0
print("QQQL single-family nonzero terms",len(nonzero),"example",nonzero[0])

# Left-handed fields. Each entry lists hypercharge and baryon number.
charges = {"Q":(Q(1,6),Q(1,3)), "u":(Q(-2,3),Q(-1,3)),
           "d":(Q(1,3),Q(-1,3)), "L":(Q(-1,2),0),
           "e":(1,0),"n":(0,0),"H":(Q(1,2),0)}
x = {k:2*y-b for k,(y,b) in charges.items()}
assert all(Q(v).denominator==1 for v in x.values())
assert (x["Q"]+x["u"]+x["H"]) % 3 == 0
assert (x["Q"]+x["d"]-x["H"]) % 3 == 0
assert (x["L"]+x["e"]-x["H"]) % 3 == 0
assert (x["L"]+x["n"]+x["H"]) % 3 == 0
assert (2*(x["L"]+x["H"]))%3 == 0 and (3*x["Q"]+x["L"])%3 != 0
print("2Y-B charges modulo three", {k:int(z%3) for k,z in x.items()})

# Sufficient same-scale inputs. Strange and down have the uncertainty from
# the gauge input. Use conservative bounds instead of floating-point fits.
ac = Q(".00311")**2/(Q(".8616")**2-Q(".00000615")**2)
bs = Q(".000271")**2/(Q(".01401")**2-Q(".0000137")**2)
comm_sq = 2*((1-ac)*(Q(".035")*Q(".98")-bs))**2
assert comm_sq>Q(".0479")**2 and 3*comm_sq/4>Q(".0414")**2
print("derived commutator and cyclic symmetry-breaking lower bounds",
      sqrt(float(comm_sq)),sqrt(float(3*comm_sq/4)))

# Enumerate every family-independent Z3 charge with one Higgs and a
# Majorana Weinberg term. Yukawa invariance fixes the conjugate singlets.
for quark,higgs in product(range(3), repeat=2):
    lepton = -higgs
    z = dict(Q=quark, H=higgs, L=lepton, u=-quark-higgs, d=-quark+higgs,
             e=-lepton+higgs, n=-lepton-higgs)
    assert all((z[k]-higgs*x[k]-quark*3*charges[k][1])%3==0 for k in z)
    assert ((3*quark+lepton)%3 != 0) == (higgs != 0)

# Add a continuous gauged X = 3B - 9 L_tau and one neutral singlet
# conjugate neutrino of X charge +9 (the other singlet neutrinos may be neutral).
# Each tuple is multiplicity, hypercharge, X, 2T_color with weak multiplicity,
# and 2T_weak with color multiplicity.
rows = []
for generation in range(3):
    tau = int(generation == 2)
    rows.extend([(6,Q(1,6),1,2,3), (3,Q(-2,3),-1,1,0), (3,Q(1,3),-1,1,0),
                 (2,Q(-1,2),-9*tau,0,1), (1,Q(1),9*tau,0,0), (1,Q(0),9*tau,0,0)])
anomalies = [sum(tc*xq for _,y,xq,tc,tw in rows), sum(tw*xq for _,y,xq,tc,tw in rows),
             sum(m*y*y*xq for m,y,xq,_,_ in rows), sum(m*y*xq*xq for m,y,xq,_,_ in rows),
             sum(m*xq**3 for m,y,xq,_,_ in rows), sum(m*xq for m,y,xq,_,_ in rows)]
assert anomalies == [0]*6
print("continuous X anomaly sums", list(map(int,anomalies)))
