"""Lowest-dimensional local violation of baryon triality's continuous parent baryon number."""
from itertools import permutations
from collections import defaultdict


def monomial(species,color,spin):
    # Two spin components for each antiquark color, followed by the charged antilepton.
    bit=12+spin if species=="e" else 6*(species=="d")+2*color+spin
    return {1<<bit:1}


def wedge(a,b):
    product=defaultdict(int)
    for i,x in a.items():
        for j,y in b.items():
            if i&j:
                continue
            sign=1
            bits=i
            while bits:
                t=bits&-bits
                if ((j&(t-1)).bit_count()&1):
                    sign=-sign
                bits-=t
            product[i|j]+=sign*x*y
    return {k:v for k,v in product.items() if v}


def sum_terms(terms):
    total=defaultdict(int)
    for weight,poly in terms:
        for k,v in poly.items():
            total[k]+=weight*v
    return {k:v for k,v in total.items() if v}


def baryon(last,spin):
    terms=[]
    for a,b,c in permutations(range(3)):
        parity=(-1)**((a>b)+(a>c)+(b>c))
        for v in (0,1):
            ud=wedge(monomial("u",a,v),monomial("d",b,1-v))
            terms.append((parity*(-1)**v,wedge(ud,monomial(last,c,spin))))
    return sum_terms(terms)


if __name__=="__main__":
    neutron=[baryon("d",j) for j in range(2)]
    proton=[baryon("u",j) for j in range(2)]
    nn=sum_terms([((-1)**j,wedge(neutron[j],neutron[1-j])) for j in range(2)])
    pe=sum_terms([((-1)**j,wedge(proton[j],monomial("e",0,1-j))) for j in range(2)])
    result=wedge(nn,pe)
    assert result and all(mask.bit_count()==10 for mask in result)
    assert all((mask&63).bit_count()==4 and ((mask>>6)&63).bit_count()==5
               and (mask>>12).bit_count()==1 for mask in result)
    assert 4*(-4)+5*2+6==0    # integer hypercharges 6Y
    assert 9*(-1)+9==0        # X=3B-9 L_tau for the tau singlet
    print("surviving exterior monomials, largest coefficient, one bitmask/coefficient",
          len(result),max(abs(v) for v in result.values()),min(result.items()))
