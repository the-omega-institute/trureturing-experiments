"""Finite separation of the shell phase at 2 log(6) from 45 RH resolvent readings."""
from math import comb, prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx

read = run_path(str(Path(__file__).with_name("ch39-finite-nbb-coherent-four-windows.py")))["readings"]
# A positive index selects that zero (one-based below), 0 and negative
# indices select multiples of the tail endpoint divided by 64.
anchors = [0, 1, 2, 5, 3, 12, 7, 4, 30, 6, 18, 68, 9, 8, 23, 11, 139, 15, 10, 43, 14, 269, 20, 36, 13, 90, 16, 26, 53, 17, -36, 21, 33, 19, 110, 25, 48, 189, 22, 28, 61, 39, -13, 24, 78]
values = [-1077602819, 994558093, 1071505394, 227719942, -98105971, 374780121, -560295027, -643996656, -284989526, -989990697, 859085670, -1106848664, -778157751, -260923710, -503036655, 261711770, -1129559053, 705311114, -817600663, -286771120, -367051574, -1079754414, 1073741824, -469411971, 618441427, -1128449154, -53596090, 9071048, -1069854741, -548175477, -1078507828, 7676766, 380910004, 475335466, -1070615403, -626116400, -1004722699, -1040648910, -245234688, -467221578, -933237520, 97031554, -1076565619, 40893035, -1059302299]
cosine = -1073797225


def positive_bernstein(coefficients, depth):
    if all(c>0 for c in coefficients):
        return
    assert depth>0
    row = coefficients
    left,right = [row[0]],[row[-1]]
    while len(row)>1:
        row = [(u+v)/2 for u,v in zip(row,row[1:])]
        left.append(row[0])
        right.append(row[-1])
    positive_bernstein(left,depth-1)
    positive_bernstein(right[::-1],depth-1)


def main():
    ctx.prec,ctx.cap=1792,3
    degree=len(anchors)-1
    # Zeros only need enough accuracy for the interpolation and strict signs.
    with ctx.workprec(384):
        gs=[r.imag for r in acb.zeta_zeros(1,360)]
    ze=[1/g**2 for g in gs]
    end=ze[-1]
    nodes=[ze[i-1] if i>0 else -i*end/64 for i in anchors]
    values_at_nodes=[arb(k)/2**30 for k in values]
    c=arb(cosine)/2**30
    # For the continuous tail interval use the common polynomial numerator.
    x=arb_poly([0,1])
    denominator=prod((1+j*j*end*x for j in range(1,degree+1)),start=arb_poly([1]))
    numerator=arb_poly([0])
    for i,node in enumerate(nodes):
        factors=[(end*x-other)*(1/(node-other)) for k,other in enumerate(nodes) if k!=i]
        numerator+=values_at_nodes[i]*prod((1+j*j*node for j in range(1,degree+1)),start=arb(1))*prod(factors,start=arb_poly([1]))
    # At the RH target use partial fractions, with higher precision for the Cauchy matrix.
    with ctx.workprec(1792):
        ref_nodes=[1/acb.zeta_zero(i).imag**2 if i>0 else -i/(64*acb.zeta_zero(360).imag**2) for i in anchors]
        V=arb_mat([[1/(1+j*j*z) for z in ref_nodes] for j in range(degree+1)])
        coeff=V.transpose().solve(arb_mat([[v] for v in values_at_nodes])).entries()
        target=sum(a*f for a,f in zip(coeff,read(degree+1)))
    assert target>arb(".005")
    assert sum(abs(a) for a in coeff[1:])<arb(2)**224
    # At actual zeros evaluate Lagrange cardinal products without polynomial cancellation.
    scales = [prod((1+j*j*node for j in range(1,degree+1)),start=arb(1)) for node in nodes]
    for g,z in zip(gs,ze):
        r=arb(0)
        w=prod((1+j*j*z for j in range(1,degree+1)),start=arb(1))
        for i,node in enumerate(nodes):
            basis=prod(((z-o)/(node-o) for k,o in enumerate(nodes) if i!=k),start=arb(1))
            r+=values_at_nodes[i]*scales[i]*basis/w
        assert r+c*(2*arb(6).log()*g).cos()<arb(".002")
    # For all larger heights use |cos|<=1 and verify the numerator of
    # .002-|c|-R is positive on the *entire* unit interval after z=end*x.
    poly=(arb(".002")-abs(c))*denominator-numerator
    bernstein=[sum(poly[k]*comb(i,k)/comb(degree,k) for k in range(i+1)) for i in range(degree+1)]
    positive_bernstein(bernstein,14)
    print("At integer offsets 0..44 the doubled-clock shell cannot match all RH readings.")
    print("Target",target.str(12),"pointwise upper bound .002; cosine coefficient",c.str(12))
    print("log2 sum of absolute nonconstant partial fractions",(sum(abs(a) for a in coeff[1:]).log()/arb(2).log()).str(10))


if __name__=="__main__":
    main()
