"""Global one-dimensional fits, using integral difference constraints for every ordering."""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

fixed=run_path(str(Path(__file__).with_name("ch39-rh-fixed-clock-crossings.py")))
points=[p for p in product(range(5),repeat=3) if sum(x*x for x in p)==17]
labels=[[17+sum(a*b for a,b in zip(p,q)) for q in points] for p in points]


def has_fit(low,high):
    def extend(order,dist,remaining):
        if not remaining:
            return True
        positions=range(len(order))
        for v in remaining:
            upper=[high[labels[p][v]] for p in order]
            lower=[low[labels[p][v]] for p in order]
            arrive=[min(dist[i][j]+upper[j] for j in positions) for i in positions]
            if any(a<b for a,b in zip(arrive,lower)):
                continue
            depart=[min(dist[j][i]-lower[j] for j in positions) for i in positions]
            distances=[[min(dist[i][j],arrive[i]+depart[j]) for j in positions]+[arrive[i]]
                       for i in positions]+[depart+[0]]
            if extend(order+[v],distances,remaining-{v}):
                return True
        return False
    # Acting by all coordinate permutations leaves one representative per
    # possible orbit of the leftmost point (six outer points or three inner).
    return any(extend([v],[[0]],set(range(9))-{v}) for v in (0,4))


def main():
    ctx.prec,ctx.cap=1536,3
    z,weights,_=fixed["shell"]["quad"](21)
    keys=set(n for row in labels for n in row if n!=34)
    orders=([0,2,1,3,4,5,6,8,7],[0,2,1,4,3,5,6,8,7])
    witnesses=(["1.107201683672",".988448692086",".563650686737",".424798005350",
                ".424798005349",".424798005348",".424798005347",".138852681387","0"],
               ["40.116188640930","35.808945885667","20.880586616813","20.880586616812",
                "14.928359268860","14.928359268859","14.928359268858","5.952227347958","0"])
    for j,(lo,hi) in enumerate([(".22948510774",".22948510779"),(".17988096371",".17988096377")]):
        e=arb(lo)
        def read(length):
            t=length**2
            if j==0:
                return 1/(1+t)
            return sum(w/(1+n*t) for w,n in zip(weights,z))
        def length_for(value):
            if value>=1:
                return arb(0)
            a,b=arb(0),arb(128)
            assert read(b)<value<1
            for _ in range(82):
                mid=(a+b)/2
                gap=read(mid)-value
                if gap>0:a=mid
                elif gap<0:b=mid
                else:break
            return a.union(b)
        scale=2**55
        low={n:int((scale*length_for(arb(n)/34+e)).lower().floor().unique_fmpz()) for n in keys}
        high={n:int((scale*length_for(arb(n)/34-e)).upper().ceil().unique_fmpz()) for n in keys}
        assert not has_fit(low,high)
        x={i:arb(v) for i,v in zip(orders[j],witnesses[j])}
        for i,k in product(range(9),repeat=2):
            if i!=k:
                assert abs(x[i]-x[k])>0
                assert abs(read(x[i]-x[k])-arb(labels[i][k])/34)<arb(hi)
        print("Single pole" if j==0 else "Arithmetic 21 poles","line minimax error bracket",lo,hi)


if __name__=="__main__":
    main()
