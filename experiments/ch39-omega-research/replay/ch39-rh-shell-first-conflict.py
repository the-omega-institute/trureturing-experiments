"""Locate the first inconsistent consecutive RH sample for the doubled shell clock."""
from math import comb, prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx

directory=Path(__file__).parent
read=run_path(str(directory/"ch39-finite-nbb-coherent-four-windows.py"))["readings"]
positive=run_path(str(directory/"ch39-rh-doubled-clock-resolution.py"))["positive_bernstein"]
support = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25, 26, 28, 29, 31, 34, 38, 42, 53, 63, 71, 85, 104, 107, 128, 166, 182, 207, 243, 422, 989, 2693]
anchors = [0, 1, 2, 5, 3, 12, 7, 4, 30, 6, 18, 68, 9, 8, 23, 11, 139, 15, 10, 43, 14, 269, 20, 36, 13, 90, 16, 26, 53, 17, 522, 21, 33, 19, 110, 25, 48, 189, 22, 28, 61, 39]
values = [-61657476612, 62213032658, 53648962878, 7985478015, -8914350015, 16372370422, -29071494509, -37262451856, 4345503910, -51752256491, 50920647772, -55491463489, -44319909615, -16974927978, -29991601438, 9829994927, -58848203348, 35296670215, -46461860314, -23958650330, -22930861090, -59603966720, 68719476736, -32774950484, 36178765046, -57370600452, -6871054039, -3168421538, -54136371490, -32141901816, -61175328776, -3307841777, 14112371808, 21209315921, -60579426929, -35455588146, -58183964424, -59172306445, -10890453266, -52959588823, -51927216528, 6266464269]
phases = [-3704581000, -4174627109, -1899926099, -59087227541]


def bernstein_test(poly, degree):
    coefficients=[sum(poly[j]*comb(i,j)/comb(degree,j) for j in range(i+1)) for i in range(degree+1)]
    positive(coefficients,16)


def main():
    ctx.prec,ctx.cap=1792,3
    levels=read(42)
    logs=[arb(2).log(),arb(3).log(),(arb(3)/2).log(),arb(6).log()]
    def shell(g):
        return [(2*g*l).cos() for l in logs]
    def column(g):
        return [1/(1+j*j/g**2) for j in range(41)]+shell(g)
    low=arb_mat([column(acb.zeta_zero(i).imag) for i in support]).transpose()
    weights=low.solve(arb_mat([[f] for f in levels[:41]]+[[0]]*4)).entries()
    assert min(weights)>arb(".00128") and sum(weights).contains(1)
    nodes=[1/acb.zeta_zero(i).imag**2 if i else arb(0) for i in anchors]
    h=[arb(k)/2**36 for k in values]
    c=[arb(k)/2**36 for k in phases]
    degree=len(nodes)-1
    V=arb_mat([[1/(1+j*j*z) for z in nodes] for j in range(degree+1)])
    coeff=V.transpose().solve(arb_mat([[v] for v in h])).entries()
    target=sum(a*f for a,f in zip(coeff,levels))
    gap=arb(".00005")
    assert target>arb(".00020") and sum(map(abs,c))<2
    cost=sum(map(abs,coeff[1:]))
    assert cost<arb(2)**210
    assert 2*arb(2)**-16+cost*arb(2)**-224<target-gap
    # The two independent phase angles share the opposite and sum angles.
    # Maximize over beta first, then bound the square root in X=cos(alpha).
    x=arb_poly([0,1])
    X=2*x-1
    torus=arb(".8944")
    c1,c2,c3,c4=c
    assert torus>abs(c1)
    bound=(torus-c1*X)**2-(c2+(c3+c4)*X)**2-(c3-c4)**2*(1-X**2)
    bernstein_test(bound,2)

    with ctx.workprec(256):
        # Recompute the indexed actual heights at the precision needed to
        # check the pointwise bound. Index 3200 begins the continuous tail.
        gs=[r.imag for r in acb.zeta_zeros(1,3200)]
        scaled=[]
        for i,z in enumerate(nodes):
            w=prod((1+j*j*z for j in range(1,degree+1)),start=arb(1))
            lag=prod((z-other for k,other in enumerate(nodes) if k!=i),start=arb(1))
            scaled.append(h[i]*w/lag)
        known=dict(zip(anchors,h))
        for i,g in enumerate(gs,1):
            if i in known:
                r=known[i]
            else:
                z=1/g**2
                factors=[z-o for o in nodes]
                w=prod((1+j*j*z for j in range(1,degree+1)),start=arb(1))
                r=prod(factors)*sum(v/f for v,f in zip(scaled,factors))/w
            assert r+sum(v*u for v,u in zip(c,shell(g)))<gap
        end=1/gs[-1]**2
    # For z=end*x the common denominator has positive sign. Test the
    # polynomial numerator of gap - torus - R on the entire unit interval.
    denom=prod((1+j*j*end*x for j in range(1,degree+1)),start=arb_poly([1]))
    nom=arb_poly([0])
    for i,d in enumerate(scaled):
        nom+=d*prod((end*x-o for k,o in enumerate(nodes) if k!=i),start=arb_poly([1]))
    bernstein_test((gap-torus)*denom-nom,degree)
    print("Offsets 0..40 are jointly feasible, min weight",min(weights).str(12))
    print("Offset 41 first breaks the joint doubled shell. Target",target.str(12),"pointwise upper .00005.")
    print("log2 resolvent coefficient cost",(cost.log()/arb(2).log()).str(12))


if __name__=="__main__":
    main()
