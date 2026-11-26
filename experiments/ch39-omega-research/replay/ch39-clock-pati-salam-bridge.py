"""Threshold matching and one-loop beta coefficients of the finite semisimple sector."""
from fractions import Fraction as F
from flint import arb

# Dimension, index for each simple factor, using T(fundamental)=1/2.
dim_adj = (8,15,3,3,3,8)  # X, color SU(4), weak L, R, spectator S, family
cas_adj = (3,4,2,2,2,3)
one, two, three, four = (1,0), (2,F(1,2)), (3,F(1,2)), (4,F(1,2))
spin1, spin32, spin2 = (3,2), (4,5), (5,10)
adjfamily = (8,3)
fields = (
    (F(2,3),(one,four,two,one,one,three)),  # chiral left
    (F(2,3),(one,four,one,two,one,three)),  # conjugate right
    (F(1,3),(one,four,one,two,one,three)),  # PS breaking
    (F(1,3),(three,one,one,one,spin2,one)),
    (F(1,3),(one,one,one,spin1,spin32,one)),
    (F(1,3),(one,one,two,two,one,one)),  # complex electroweak bidoublet
    (2*F(1,6),(one,one,one,one,one,adjfamily)),
)
from math import prod
beta = [-F(11,3)*c for c in cas_adj]
for weight, reps in fields:
    dimension = prod(d for d,t in reps)
    for index,(d,t) in enumerate(reps):
        beta[index] += weight*t*(dimension//d)
assert beta == [-F(61,6),-F(29,3),-3,F(5,3),F(23,3),-F(10,3)]
ratios = (8,1,1,6,F(8,3),1)  # g_i^2/g_U^2
norms = (F(16,3),F(2,3),0,1,F(4,9),0)
assert sum(F(w)/r for w,r in zip(norms,ratios)) == F(5,3)
def ball(q):
    q = F(q)
    return arb(q.numerator)/q.denominator
g2max = 4*arb.pi()/arb("40.6142")
def effective_charge(reps):
    return sum(F(t * ad, d)*r for (d,t),ad,r in zip(reps,dim_adj,ratios))
charges = [effective_charge(reps) for _,reps in fields]
assert max(charges) == F(80,3)
assert ball(max(charges))*g2max < arb("8.252")
for b,r in zip(beta,ratios):
    assert 1-ball(b*r)*g2max*arb(559).log()/(8*arb.pi()**2) > 0
running = [ball(r)*g2max/(1-ball(b*r)*g2max*arb(559).log()/(8*arb.pi()**2)) for b,r in zip(beta,ratios)]
assert max(sum(ball(F(t*ad,d))*g for (d,t),ad,g in zip(reps,dim_adj,running))
           for _,reps in fields) < 4*arb.pi()
print("finite-sector beta",beta,"maximum breaker charge at Q",(ball(max(charges))*g2max).upper().str(13))
print("squared couplings of finite sector at 559Q on upper input edge",
      [g.upper().str(10) for g in running])
