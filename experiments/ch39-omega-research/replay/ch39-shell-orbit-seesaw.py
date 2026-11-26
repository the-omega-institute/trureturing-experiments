"""Check the mass-gap and reactor-angle tests along the interpolation from averaged to marked orbits."""
from itertools import permutations
from pathlib import Path
from runpy import run_path
from flint import arb,arb_mat,acb_mat

m=run_path(str(Path(__file__).with_name("ch39-shell-four-seesaw.py")))
ind=range(3)


def eigen(M):
    es=sorted((x.real for x in acb_mat(M).eig()),key=lambda x:float(x.mid()))
    assert es[0]<es[1]<es[2]
    return es


def main():
    for n,trial in ((25,".85"),(255,".678"),(719,".662")):
        marked=sum(m["matrices"](n),arb_mat(3,3))/4
        T=arb_mat([[sum(marked[p[i],p[j]] for p in permutations(ind))/6 for j in ind] for i in ind])
        D=marked-T
        norm=max(sum(abs(D[i,j]) for j in ind).upper() for i in ind)
        # Verify a normal-ordering mass ratio in a narrow central data range.
        ev=eigen(T+D*arb(trial))
        sq=[1/t**2 for t in ev[::-1]]
        example=(sq[1]-sq[0])/(sq[2]-sq[0])
        assert arb(".028")<example<arb(".031")
        leaves=0
        def check(lo,hi):
            nonlocal leaves
            mid=(lo+hi)/2
            err=(hi-lo)*norm/2
            ev=[(t-err).union(t+err) for t in eigen(T+D*mid)]
            if ev[0]>0:
                squares=[1/t**2 for t in ev[::-1]]
                if squares[2]>squares[0]:
                    gaps=[(squares[i+1]-squares[i])/(squares[2]-squares[0]) for i in (0,1)]
                    H=T+D*lo.union(hi)
                    H2=H*H
                    permitted=False
                    for gap,isolated in zip(gaps,(0,2)):
                        if gap<arb(".02655") or gap>arb(".034"):
                            continue
                        a,b=[x for j,x in enumerate(ev) if j!=isolated]
                        denom=(ev[isolated]-a)*(ev[isolated]-b)
                        if denom.contains(0):
                            permitted=True
                            break
                        pp=[(H2[i,i]-(a+b)*H[i,i]+a*b)/denom for i in ind]
                        if any(not (x>arb(".025")) for x in pp):
                            permitted=True
                            break
                    if not permitted:
                        leaves+=1
                        return
            assert hi-lo>arb(1)/2**22
            check(lo,mid)
            check(mid,hi)
        check(arb(0),arb(1))
        print("zero",n,"eta",trial,"mass ratio",example.str(13),"simultaneous reactor exclusion, boxes",leaves)


if __name__=="__main__":
    main()
