"""An all-representation congruence for the charge-two clock lift and its residual Z9."""
from itertools import product

mod=729
def low(y,h):
    tr=h%3
    a0=(y-4*h)%36
    v=y-a0
    return v**3-6*a0*v*v+144*tr*tr*v
def odd(y,h):
    return (low(y,h)-low(-y,-h))*pow(2,-1,mod)%mod

fields=[(6,1,1,0),(3,-4,-1,0),(3,2,-1,0),
        (2,-3,0,1),(1,6,0,-1),(1,0,0,-1)]
def invariant(charges):
    return sum(3*dim*odd(y,h) for (dim,y,_,_),h in zip(fields,charges))%mod

def check():
    for p,q in product(range(11),repeat=2):
        weight=[]
        for i in range(q,p+q+1):
            for j in range(q+1):
                weight.extend([3*(i+j)-2*(p+2*q)]*(i-j+1))
        rank=len(weight)
        tau=weight[0]%3
        assert sum(weight)==0 and (tau==0 or rank%3==0)
        H2=sum(h*h for h in weight)
        H3=sum(h**3 for h in weight)
        assert (H2+rank*tau*tau)%9==0
        for a in (-73,-36,-18,-1,0,1,17,36,74):
            dim=3 if (a+4*tau)%3 else 1
            direct=dim*sum(odd(a+4*h,h) for h in weight)
            assert (direct-32*dim*(2*H3-3*a*H2))%mod==0
    baryon=[invariant([b+k*y for _,y,b,_ in fields]) for k in range(9)]
    lepton=invariant([3*l for _,_,_,l in fields])
    print("baryon residual invariant across integer hypercharge representatives",baryon)
    print("lepton residual invariant",lepton)
    candidates=[k for k in range(9) if (1+k)%3==0]
    assert candidates==[2,5,8]
    assert all(baryon[k]==486 for k in candidates) and lepton==0
if __name__=="__main__":
    check()
