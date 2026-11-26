"""Exact rational weak weights with identical derivatives through order seven and positive transfer."""
from fractions import Fraction as F
from math import comb

def size(j):
    return 2*j+1
def curvature(j):
    return F(4*j*(j+1),3)
def multiply_s(p):
    # s=(3/8)(1-phi_1), where phi_j is the normalized integer-spin character.
    out=[F(0)]*(len(p)+1)
    for j,weight in enumerate(p):
        out[j] += F(3,8)*weight
        for ell in range(abs(1-j),j+2):
            out[ell] -= weight*F(size(ell),8*size(j))
    return out
def moment(p,k):
    return sum(weight*sum((2*t)**(2*k) for t in range(-j,j+1))/size(j) for j,weight in enumerate(p))
def blocked(p):
    coefficients = [weight**2/size(j)**2 for j,weight in enumerate(p)]
    return sum(curvature(j)*w for j,w in enumerate(coefficients))/sum(coefficients)

f=[F(1)]
for _ in range(4):
    f=multiply_s(f)
base=[F(0),F(0),F(1,100),F(1,300),F(1,500)]
base[1]=(1-sum(w*curvature(j) for j,w in enumerate(base)))/curvature(1)
base[0]=1-sum(base)
a,b=([w+sign*t/8 for w,t in zip(base,f)] for sign in (-1,1))
assert min(a+b)>0 and min(a[0],b[0])>F(1,2)
assert all(moment(a,k)==moment(b,k) for k in range(4))
assert moment(a,1)==1 and moment(a,4)!=moment(b,4)
assert blocked(a)!=blocked(b)
print('exact blocked curvatures',blocked(a),blocked(b))
print('common microscopic curvature and absolute blocked gap',moment(a,1),float(abs(blocked(b)-blocked(a))))
