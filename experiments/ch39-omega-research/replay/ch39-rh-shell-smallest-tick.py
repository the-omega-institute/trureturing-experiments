"""For 42 exact RH samples the smallest common integer shell tick is 5; its first feasible zero prefix is 1279."""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, arb_poly, arb_mat, acb, ctx

prior=run_path(str(Path(__file__).with_name("ch39-rh-shell-first-conflict.py")))
support = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 28, 29, 32, 34, 38, 39, 40, 41, 48, 54, 58, 68, 69, 91, 98, 106, 146, 158, 278, 291, 312, 1279]
witnesses = [(3, 26, [0, 1, 2, 5, 3, 12, 7, 4, 30, 6, 18, 68, 9, 8, 23, 11, 139, 15, 10, 43, 14, 269, 20, 36, 13, 90, 16, 26, 53, 17, 522], [-67444330, 56299435, 67108864, 30083529, -53590760, -16575738, 67108864, 66334545, -106357988, 37794843, -38197066, -67617105, 61137774, 61439731, -67323087, -27211279, -67845164, -19321090, -58740159, -44504465, 65036466, -67004268, -9884310, -75663768, 8944742, -68622282, -50896097, 47099529, -69842211, -67429761, -68606396, 0, 0, -130652, -67313678]), (4, 26, [0, 1, 2, 5, 3, 12, 7, 4, 30, 6, 18, 68, 9, 8, 23, 11, 139, 15, 10, 43, 14, 269, 20], [-67813221, 39825007, 60138745, 66525661, 65019999, 61101515, 62465996, -40860504, -76735344, 19780428, -70359726, -67019392, 52427895, -26844427, -50466794, 30693996, -67598014, -23455164, -38547543, -64231560, 16496731, -66936230, -57251704, -292527, 240960, 67108864, -64931]), (5, 44, [0, 1, 2, 5, 3, 12, 7, 4, 30, 6, 18, 68, 9, 8, 23, 11, 139, 15, 10, 43, 14, 269, 20, 36, 13, 90, 16, 26, 53, 17, 522, 21, 33, 19, 110, 25, 48, 189, 22, 28, 61, 39], [17592186044416, 293263386643, 515431803493, 124620032943, -355379316884, -163354252681, -81863473212, 8309886589, -73007823700, -375592839456, -257395944566, -272769617535, -361218327222, -259531163393, -1710692676562, -182248267003, -620269079469, -222006360766, 585331553115, 166624886570, -14338044982, -731149627182, 576377560803, -785704116642, -379086198510, -392604264093, -78563080717, -2187085283452, -128208650073, 133678273440, -8309107502860, 354800925468, -175802047294, 134965368898, -958287097583, -95379809240, -237355737615, -1667765926560, 326727328570, -285261105144, -407067016163, -92510566826, 107591999353, 266103574645, 22196855373, -242206864474])]


def main():
    ctx.prec,ctx.cap=1536,3
    levels=prior["read"](42)
    logs=[arb(t).log() for t in (2,3,arb(3)/2,6)]
    def phase(g,k):
        return [(g*k*l).cos() for l in logs]
    def col(i):
        g=acb.zeta_zero(i).imag
        return [1/(1+j*j/g**2) for j in range(42)]+phase(g,5)
    M=arb_mat([col(i) for i in support]).transpose()
    weights=M.solve(arb_mat([[f] for f in levels]+[[0]]*4)).entries()
    assert min(weights)>arb(".00079") and sum(weights).contains(1) and max(support)==1279
    with ctx.workprec(256):
        heights=[r.imag for r in acb.zeta_zeros(1,3200)]
    for k,bits,indices,values in witnesses:
        h=[arb(v)/2**bits for v in values[:-4]]
        c=[arb(v)/2**bits for v in values[-4:]]
        deg=len(indices)-1
        nodes=[1/acb.zeta_zero(i).imag**2 if i else arb(0) for i in indices]
        W=lambda z:prod((1+j*j*z for j in range(1,deg+1)),start=arb(1))
        V=arb_mat([[1/(1+j*j*z) for j in range(deg+1)] for z in nodes])
        a=V.solve(arb_mat([[v] for v in h])).entries()
        target=sum(v*f for v,f in zip(a,levels))
        bound=arb(".0001")
        assert target>2*bound
        factors=[v*W(z)/prod((z-o for j,o in enumerate(nodes) if j!=i),start=arb(1))
                 for i,(v,z) in enumerate(zip(h,nodes))]
        known=dict(zip(indices,h))
        with ctx.workprec(256):
            for i,g in enumerate(heights[:1278] if k==5 else heights,1):
                z=1/g**2
                if i in known:
                    r=known[i]
                else:
                    diffs=[z-o for o in nodes]
                    r=prod(diffs)*sum(v/d for v,d in zip(factors,diffs))/W(z)
                assert r+sum(v*u for v,u in zip(c,phase(g,k)))<bound
        if k!=5:
            end=1/heights[-1]**2
            x=arb_poly([0,1])
            D=prod((1+j*j*end*x for j in range(1,deg+1)),start=arb_poly([1]))
            N=arb_poly([0])
            for i,v in enumerate(factors):
                N+=v*prod((end*x-o for j,o in enumerate(nodes) if j!=i),start=arb_poly([1]))
            prior["bernstein_test"]((bound-sum(map(abs,c)))*D-N,deg)
        print("Tick",k,"separator at offsets 0..",deg,"target",target.str(8),"upper bound .0001",
              "on first 1278 heights" if k==5 else "on every actual height")
    print("Tick 5 gives all 42 readings and the full shell on 46 positive zero pairs, min weight",min(weights).str(10))


if __name__=="__main__":
    main()
