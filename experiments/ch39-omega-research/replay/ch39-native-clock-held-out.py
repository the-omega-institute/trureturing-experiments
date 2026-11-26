"""Test a held-out RH reading and certify a 50-reading positive realization of the same generator."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def main(cases=None, minimum=None, ticks=5):
    root=Path(__file__).resolve().parent
    old=run_path(str(root/"ch39-rh-shell-smallest-tick.py"))
    law=run_path(str(root/"ch39-two-prime-height-law.py"))
    words=run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    ctx.prec,ctx.cap=2048,3
    targets=old["prior"]["read"](max(length for _,length in cases)+1 if cases else 51)
    extended=[*range(1,27),28,29,30,32,34,36,39,43,48,54,58,62,68,79,91,106,127,146,204,276,340,375,462,699,1828,3133,3170,3191]
    a=law["rate_numerators"](ticks)
    for support,length in (cases or [(old["support"],42),(extended,50)]):
        n=len(support)
        assert n==length+4
        heights=[acb.zeta_zero(i).imag for i in support]
        def row(j):
            return [1/(1+(j/g)**2) for g in heights]
        phase=[[] for _ in range(4)]
        for g in heights:
            h=int((g*2**288).floor().unique_fmpz())
            u,v=[-43235114615559591407865398948523*((h*s+2**311)//2**312) for s in a]
            for dest,k in zip(phase,(u,v,v-u,v+u)):
                lo,hi=words["power"](k)[0]
                dest.append((arb(lo)/words["DEN"]).union(arb(hi)/words["DEN"]))
        rows=[row(j) for j in range(length)]+phase
        mid=arb_mat([[v.mid() for v in r] for r in rows])
        inv=mid.inv()
        error=inv*arb_mat(rows)-arb_mat([[int(i==j) for j in range(n)] for i in range(n)])
        norm=max(sum(abs(error[i,j]).upper() for j in range(n)).upper() for i in range(n))
        assert norm<arb(2)**-200
        guess=(inv*arb_mat([[v] for v in targets[:length]]+[[0]]*4)).entries()
        change=norm*max(abs(x).upper() for x in guess)/(1-norm)
        assert min(x.lower() for x in guess)-change>(arb(".00069") if minimum is None else minimum)
        prediction=sum(w*t for w,t in zip(guess,row(length)))-targets[length]
        gap=prediction+arb(0,n*change)
        if length==42:
            assert arb("-3.891835e-49")<gap<arb("-3.891834e-49")
        print("least weight >", (min(x.lower() for x in guess)-change).lower().str(10), "exact nodes",length,"support",n,"next predicted minus RH",gap.str(12))


if __name__=="__main__":
    main()

