"""The attained smallest discrepancy at offset 41 under the 41 preceding readings and shell orthogonality."""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx

prior=run_path(str(Path(__file__).with_name("ch39-rh-shell-first-conflict.py")))
support=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
         25,26,29,31,34,39,42,49,53,61,63,75,82,104,113,153,167,304,319,417,902,2693)


def main():
    ctx.prec,ctx.cap=1792,3
    degrees=range(42)
    levels=prior["read"](42)
    gs=[acb.zeta_zero(i).imag for i in support]
    logs=[(arb(t).log()*2) for t in (2,3,arb(3)/2,6)]
    def phase(g):
        return [(k*g).cos() for k in logs]
    def moments(g):
        return [1/(1+j*j/g**2) for j in degrees]
    rows=[moments(g) for g in gs]
    mat=arb_mat([r[:-1]+phase(g) for r,g in zip(rows,gs)])
    w=mat.transpose().solve(arb_mat([[v] for v in levels[:-1]]+[[0]]*4)).entries()
    assert min(w)>arb(".002") and sum(w).contains(1)
    discrepancy=sum(p*r[-1] for p,r in zip(w,rows))-levels[-1]
    assert discrepancy>0
    # The unique functional interpolating f_41 at these 45 columns is a
    # minorant on every other actual-height column. Prove the sign below.
    dual=mat.solve(arb_mat([[r[-1]] for r in rows])).entries()
    scale=arb(10)**45
    coeff=[scale*a for a in dual[:-4]]+[-scale]
    c=[scale*v for v in dual[-4:]]
    denom=lambda z:prod((1+j*j*z for j in degrees),start=arb(1))
    rational=lambda z:sum(a/(1+j*j*z) for j,a in enumerate(coeff))
    # Direct partial fractions here are evaluated at full precision once
    # per interpolation node; cardinal products then need fewer bits.
    an=prior["anchors"]
    nodes=[1/acb.zeta_zero(i).imag**2 if i else arb(0) for i in an]
    heights_to_r=dict(zip(an,map(rational,nodes)))
    scaling=[]
    for i,z in enumerate(nodes):
        ds=prod((z-o for k,o in enumerate(nodes) if k!=i),start=arb(1))
        scaling.append(heights_to_r[an[i]]*denom(z)/ds)
    with ctx.workprec(256):
        low=[r.imag for r in acb.zeta_zeros(1,3200)]
        for i,g in enumerate(low,1):
            if i in support:
                continue
            z=1/g**2
            if i in heights_to_r:
                val=heights_to_r[i]
            else:
                diffs=[z-o for o in nodes]
                val=prod(diffs)*sum(v/t for v,t in zip(scaling,diffs))/denom(z)
            assert val+sum(v*t for v,t in zip(c,phase(g)))<0
        end=1/low[-1]**2
    x=arb_poly([0,1])
    X=2*x-1
    c1,c2,c3,c4=c
    # A slight upward rounding of the value at alpha=0. The
    # Bernstein comparison covers every alpha and beta.
    torus=c1+abs(c2+c3+c4)+arb(".00001")
    assert torus>abs(c1)
    prior["bernstein_test"]((torus-c1*X)**2-(c2+(c3+c4)*X)**2-(c3-c4)**2*(1-X**2),2)
    W=prod((1+j*j*end*x for j in degrees),start=arb_poly([1]))
    P=arb_poly([0])
    for i,v in enumerate(scaling):
        P+=v*prod((end*x-o for k,o in enumerate(nodes) if k!=i),start=arb_poly([1]))
    prior["bernstein_test"](-torus*W-P,41)
    assert discrepancy>arb("9.8e-50") and discrepancy<arb("1.0e-49")
    # In the exposed simplex, move the targets along the two worst-case
    # error directions. Positivity at both endpoints certifies feasibility
    # all along their convex segment; duality gives necessity.
    shell_norm=sum(abs(v) for v in dual[-4:])
    earlier_norm=sum(abs(v) for v in dual[1:-4])
    directions=[int(v.sgn().unique_fmpz()) for v in dual]
    assert all(v in (-1,1) for v in directions)
    tau=discrepancy/shell_norm
    eps=discrepancy/(1+earlier_norm)
    def adjust(e,t):
        desired=[levels[0]]+[f-s*e for f,s in zip(levels[1:-1],directions[1:-4])]
        desired+=[-s*t for s in directions[-4:]]
        new=mat.transpose().solve(arb_mat([[v] for v in desired])).entries()
        assert min(new)>0 and sum(new).contains(1)
        assert (sum(v*r[-1] for v,r in zip(new,rows))-levels[-1]-e).contains(0)
        return min(new)
    print("Minimum max shell error with exact arithmetic",tau.str(13),"min support weight",adjust(arb(0),tau).str(7))
    print("Minimum max arithmetic error with exact shell",eps.str(13),"min support weight",adjust(eps,arb(0)).str(7))
    print("Bits of uniform arithmetic precision at this boundary",(-eps.log()/arb(2).log()).str(10))
    # Which binary place of absolute additive error reaches the gap?
    print("Minimum offset-41 excess",discrepancy.str(17),"bits",(-discrepancy.log()/arb(2).log()).str(10))
    print("The unique minimizer uses 45 actual heights with all weights >",min(w).str(10))


if __name__=="__main__":
    main()
