"""Pure-integer certificates for the required reference count at every whole-bit tolerance."""
import json
from fractions import Fraction as Q
from pathlib import Path

directory=Path(__file__).parent


def poly(coeff,x):
    result=Q(0)
    for c in reversed(coeff):
        result=result*x+c
    return result


def main():
    levels=json.loads((directory/"ch39-five-tick-record.json").read_text())
    targets=[Q(v,2**levels["bits"]) for v in levels["levels"]]
    width=Q(1,2**levels["bits"])
    data=json.loads((directory/"ch39-kinetic-precision-data.json").read_text())
    root_den=2**data["root_bits"]
    # With m directions, this is the greatest whole bit count for which
    # the absolute tolerance 2^-bits is attainable on offsets 0 through 41.
    cutoffs=[3,9,15,22,28,35,42,50,58,66,75,83,92,102,111,121,132,143,155,168]
    for m,(model,bits) in enumerate(zip(data["models"],cutoffs,strict=True),1):
        P,D=model["P"],model["Q"]
        cuts=model["roots"]
        assert len(P)==m and len(D)==m+1 and len(cuts)==m and model["order"]==m
        assert P[0]==D[0] and D[0]>0 and all(v>0 for v in D)
        assert all(a+1<b for a,b in zip(cuts,cuts[1:])) and max(cuts)+1<0
        for k in cuts:
            left,right=Q(k,root_den),Q(k+1,root_den)
            before,after=poly(D,left),poly(D,right)
            assert before*after<0
            # P's sign throughout this interval agrees with the derivative
            # of D at its simple root, so the residue and its weight are positive.
            slope=sum(i*abs(c)*abs(left)**(i-1) for i,c in enumerate(P) if i)
            pl=poly(P,left)
            assert abs(pl)>slope/root_den and pl*after>0
        errors=[]
        for j,target in enumerate(targets):
            x=Q(j*j,1681)
            e=poly(P,x)/poly(D,x)-target
            assert max(abs(e),abs(e-width))<Q(1,2**bits)
            errors.append((e-width,e))
        signs=[]
        for i in model["alternation"]:
            lower,upper=errors[i]
            threshold=Q(1,2**(bits+1))
            assert lower>threshold or upper < -threshold
            signs.append(lower>0)
        ind=model["alternation"]
        assert len(ind)==2*m and min(ind)>0 and all(a<b for a,b in zip(ind,ind[1:]))
        assert all(a!=b for a,b in zip(signs,signs[1:]))
        # Any better normalized, pole-free rational of type (m-1,m)
        # must cross P/D between each pair AND at x=0. The numerator
        # of the difference has degree at most 2*m-1.
    print("Maximum integer precision bits by reference count 1 through 20:",cutoffs)
    print("At 32, 64, 96, 128 and 169 bits, minimum counts: 6, 10, 14, 17 and 21.")


if __name__=="__main__":
    main()
