"""识别源在 Pati–Salam 公共中心商中的最小完整表示。"""
from fractions import Fraction
from math import comb
from itertools import product


def dimension(a, b, c):
    return (a+1)*(b+1)*(c+1)*(a+b+2)*(b+c+2)*(a+b+c+3)//12


def indices(a,b,c,p,left):
    dim = dimension(a,b,c)*left*(p+1)
    casimir = Fraction(3*a*a+4*b*b+3*c*c+4*a*b+2*a*c+4*b*c+12*a+16*b+12*c,8)
    return dim, dim*casimir/15, Fraction(dim*p*(p+2),12)


def weight(a,b,c,charge,triplet):
    total = a+2*b+3*c
    if (charge+3*total)%4:
        return False
    boxes = (charge+3*total)//4
    if triplet:
        if (boxes-1)%3:
            return False
        k = (boxes-1)//3
        return 0<=k<=c and c<=k<=b+c and b+c<=k+1<=a+b+c
    return b+2*c <= boxes <= total


def minima(n, triplet):
    charge, p, left = 3*n+1,n,2
    limit, u, v = indices(1,0,n,p,left)
    # Either optimum is inside this colour-dimension budget, by T_R alone.
    size_cap = max(Fraction(limit,left*(p+1)), (u+v)*12/(left*(p+1)*p*(p+2)))
    stop = 0
    while comb(stop+3,3) <= size_cap:
        stop += 1
    assert all(dimension(*labels)>size_cap for labels in ((stop,0,0),(0,stop,0),(0,0,stop)))
    rows = []
    for a,b,c in product(range(stop), repeat=3):
        if dimension(a,b,c)>size_cap or not weight(a,b,c,charge,triplet):
            continue
        assert (a+2*b+3*c+p+left-1)%2 == 0
        dim,t4,tr = indices(a,b,c,p,left)
        rows.append((dim,t4+tr,(a,b,c),t4,tr))
    winners = [min(row[i] for row in rows) for i in (0,1)]
    best = [row for row in rows if row[0]==winners[0] or row[1]==winners[1]]
    assert len(best)==1 and best[0][:2]==tuple(winners)
    return best[0]


def verify():
    for n in (9,10):
        any_weight = minima(n,False)
        fundamental = minima(n,True)
        # A colour singlet has D = 3(c-a) with b=0 by U(3) interlacing.
        rank = 3*n+1
        neutral = indices(0,0,rank,3*n,2)
        assert any_weight[2] == (0,0,n+3) and fundamental[2] == (1,0,n)
        assert neutral[0] == 2*(3*n+1)*comb(rank+3,3)
        # At every larger positive allowed multiple (>=6), even a left
        # singlet exceeds the 3-fold multiplet in each gauge index and size.
        larger = indices(0,0,2*rank,6*n,1)
        assert all(a>b for a,b in zip(larger,neutral,strict=True))
        print("n",n,"primitive dimension/T4/TR",any_weight[0],*any_weight[3:],
              "colour triplet",fundamental[0],*fundamental[3:],
              "threefold singlet",*neutral,flush=True)


if __name__ == "__main__":
    verify()
