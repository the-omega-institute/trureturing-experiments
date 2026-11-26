"""Nine finite rational link configurations determine a common 48-by-48 edge response."""
from itertools import combinations_with_replacement,product
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb,arb,ctx

shell=run_path(str(Path(__file__).with_name("ch39-finite-shell-first-triplet.py")))


def orientation():
    L=shell["matrix"]({1:1,2:1,4:1})
    adjacent=[{j for j in range(24) if L[i,j]<0} for i in range(24)]
    trail,stack=[],[0]
    while stack:
        u=stack[-1]
        if adjacent[u]:
            v=min(adjacent[u])
            adjacent[u].remove(v)
            adjacent[v].remove(u)
            stack.append(v)
        else:
            trail.append(stack.pop())
    edges=list(zip(trail,trail[1:]))
    assert len(edges)==48 and trail[0]==trail[-1]
    assert all(sum(v==i for _,v in edges)==2 for i in range(24))
    return edges


def main():
    ctx.prec=160
    edges=orientation()
    code=[106*t+t*t%53 for t in range(49)]
    differences=[a-b for a in code for b in code if a!=b]
    assert len(set(differences))==49*48
    incoming=[[0]+[code[e+1] for e,(_,v) in enumerate(edges) if v==i] for i in range(24)]
    # The tenth power is never required. At time k the link is z^(k*code).
    # z is a rational unit rotation of infinite order: 2 Re(z)=6/5 is
    # rational and nonintegral, whereas a root of unity has algebraic-integer trace.
    z=acb(arb(3)/5,arb(4)/5)
    zs={i:z**i for i in set(differences)|{0}}
    largest=arb(0)
    for i,j in combinations_with_replacement(range(24),2):
        exponents=sorted({b-a for a,b in product(incoming[i],incoming[j])})
        size=len(exponents)
        assert size==(7 if i==j else 9)
        nodes=[zs[e] for e in exponents]
        for n,u in enumerate(nodes):
            dist=prod((abs(u-v) for l,v in enumerate(nodes) if n!=l),start=arb(1))
            # Sum of the absolute coefficient values in the Lagrange cardinal.
            largest=max(largest,(2**(size-1)/dist).upper())
    power=0
    while not largest<2**power:
        power+=1
    # Off-diagonals need one cardinal coefficient. At baseline the only
    # diagonal on any actual edge enters with coefficient -1; there are 15
    # remaining products from its two incident degree-four vertex stars.
    print("Largest certified row-sum bound for the Vandermonde inverses",largest.str(10),"< 2^",power)
    print("Operator recovery factor for entrywise input error",(1+62*largest).str(10))
    assert largest<2**30


if __name__=="__main__":
    main()
