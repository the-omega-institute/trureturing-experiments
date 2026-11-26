"""Four exact arc readings obstruct periodic positive Gram transport; three can fit."""
from fractions import Fraction as Q
from flint import arb, acb, ctx
from sympy import Matrix, Rational, factor, sqrt, symbols


def main():
    c = Q(25,26)
    h = [Q(1), c]
    for _ in range(2):
        h.append(2*c*h[-1]-h[-2])
    k = [(1+x)/2 for x in h]
    gram = Matrix(4,4,lambda i,j:k[abs(i-j)])
    v = Matrix([-13,38,-38,13])
    assert gram*v == Matrix.zeros(4,1)
    assert gram[:3,:3].det()>0 and gram.rank()==3
    x = symbols("x")
    assert factor(sum(v[j]*x**j for j in range(4)))==(x-1)*(13*x**2-25*x+13)
    # A period-24 realization of the first three readings on a five-dimensional
    # positive Gram space; a rotation and its square are given by radicals.
    t = (sqrt(6)+sqrt(2))/4
    s = (sqrt(6)-sqrt(2))/4
    R = Matrix([[t,-s],[s,t]])
    assert (R.T*R).applyfunc(lambda q:q.expand()) == Matrix.eye(2)
    R2 = (R*R).applyfunc(lambda q:q.expand())
    R3 = (R2*R).applyfunc(lambda q:q.expand())
    R6 = (R3*R3).applyfunc(lambda q:q.expand())
    assert R6 == Matrix([[0,-1],[1,0]]) and R6**4==Matrix.eye(2)
    a,b = 1-t,1-R2[0,0]
    delta = Rational(1,26)
    # E[X]=delta/2, E[X^2]=delta^2/2 on X=0,a,b.
    ma = delta*(b-delta)/(2*a*(b-a))
    mb = delta*(delta-a)/(2*b*(b-a))
    assert a>delta/2 and a<delta and b>delta
    assert ma>0 and mb>0 and (1-ma-mb).simplify()>0
    assert (ma*a+mb*b-delta/2).simplify()==0
    assert (ma*a*a+mb*b*b-delta**2/2).simplify()==0
    # The period-24 witness is PSD but takes negative values at some other
    # lags. For a globally nonnegative fidelity use a finer finite cycle.
    assert (1-2*ma).simplify()<0
    ctx.prec = 256
    cosine,sine = arb(0),arb(1)
    for _ in range(7):
        cosine = ((1+cosine)/2).sqrt()
        sine = sine/(2*cosine)
    z = acb(cosine,sine)
    f,g = z**22,z**23
    a,b = 1-f.real,1-g.real
    gap = arb(1)/26
    ma = gap*(b-gap)/(2*a*(b-a))
    mb = gap*(gap-a)/(2*b*(b-a))
    assert ma>0 and mb>0 and 1-ma-mb>0
    for j in range(512):
        fidelity = 1-ma-mb+ma*(f**j).real+mb*(g**j).real
        assert fidelity>arb("0.001") and (j==0 or fidelity<1)
    # A wholly rational finite periodic kernel at positive tolerance.
    length,horizon = 127,3
    reals = [Q(1),c]
    for _ in range(length-1):
        reals.append(2*c*reals[-1]-reals[-2])
    fidelities = [(1+v)/2 for v in reals]
    cyclic = [sum(fidelities[abs((j+n)%length-j)] for j in range(length))/length
              for n in range(horizon+1)]
    for n,f in enumerate(cyclic):
        assert f == ((length-n)*fidelities[n]+n*fidelities[length-n])/length
        assert abs(f-fidelities[n]) <= Q(n,length)
    print("Exact lags 0..3:", k, "; four-state Gram null vector",list(v))
    print("Three-state Gram determinant",gram[:3,:3].det(),"; period 24 matches lags 0,1,2.")


if __name__=="__main__":
    main()
