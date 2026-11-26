"""Sharp one-loop splitting bounds for complete mirror generations."""
from fractions import Fraction as F
from itertools import product
from flint import arb, ctx


def scalar(x):
    x = F(x)
    return arb(x.numerator)/x.denominator


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def verify():
    ctx.prec = 160
    rows = [(F(1,15),1,F(2,3)),(F(8,15),0,F(1,3)),
            (F(2,15),0,F(1,3)),(F(1,5),F(1,3),0),(F(2,5),0,0)]
    beta = (F(41,10),F(-19,6),-7)
    normal = (115,-333,218)
    separator = (45,-107,62)
    assert dot(normal,beta) == 0 and sum(normal) == sum(separator) == 0
    assert [dot(normal,b) for b in rows] == [-180,134,88,-88,46]
    for on in product(range(2),repeat=5):
        v = [beta[j]+sum(bit*row[j] for bit,row in zip(on,rows)) for j in range(3)]
        assert dot(separator,v) >= 0
        assert (dot(separator,v) == 0) == (on==(1,0,0,1,0))
    entries = list(map(scalar,(".467453",".63811","1.0583")))
    def initial(factors):
        return [4*arb.pi()/(g*scalar(t))**2 for g,t in zip(entries,factors)]
    for factors in product(("0.999","1.001"),repeat=3):
        a = initial(factors)
        assert -dot(separator,a)/214 > scalar(".02565")
    a = initial((1,1,1))
    minimax = -dot(separator,a)/214
    assert scalar(".08719933") < minimax < scalar(".08719934")
    h = -2*arb.pi()*dot(normal,a)/268
    end = (2*arb.pi()*(a[0]-a[1])+scalar("16/15")*h)/scalar("109/15")
    assert h/2 < end < h
    assert scalar("6.69943e14") < 1000*end.exp() < scalar("6.69944e14")
    for k,lo,hi in ((2,1577686,1577687),(3,13552,13553)):
        assert lo < (h/k).exp() < hi
        active = [k*(h/k),0,0,k*(h/k),0]
        inverse = [a[j]-(scalar(beta[j])*end
                            +sum(x*scalar(b[j]) for x,b in zip(active,rows)))/(2*arb.pi())
                   for j in range(3)]
        assert (inverse[0]-inverse[1]).contains(0) and (inverse[1]-inverse[2]).contains(0)
    print("Sharp one-copy mismatch",minimax.str(12),"three-copy minimum mass ratio",(h/3).exp().str(12))


if __name__ == "__main__":
    verify()
