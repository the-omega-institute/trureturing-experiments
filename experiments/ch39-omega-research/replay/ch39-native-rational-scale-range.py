"""Two unequal rational common pole scales with identically the same finite observations."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path
import json


def main():
    root=Path(__file__).resolve().parent
    book=run_path(str(root/"ch39-native-positive-integers.py"))
    word=run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    law=run_path(str(root/"ch39-two-prime-height-law.py"))
    table=json.loads((root/"ch39-five-tick-record.json").read_text())
    pre=json.loads((root/"ch39-native-positive-inverse.json").read_text())
    matrix=book["enclose_matrix"](table,word,law["rate_numerators"]())
    unit=word["DEN"]
    ident=unit*2**pre["depth"]
    den=2**table["bits"]
    rhs=table["levels"] + [0]*4
    derivative=[]
    for j in range(42):
        row=[]
        for _,h,_,_,_ in table["rows"]:
            v=Q(h*h,h*h+(j*den)**2)
            n=v*(v-1)*unit
            row.append((n//1,-((-n)//1)))
        derivative.append(row)
    derivative.extend([[(0,0)]*46]*4)
    def norm(mat, subtract_identity):
        largest=0
        for i,row in enumerate(pre["numerators"]):
            total=0
            for j in range(46):
                low=high= -ident if subtract_identity and i==j else 0
                for a,record in zip(row,mat):
                    l,h=record[j]
                    if a < 0:
                        l,h=h,l
                    low+=a*l
                    high+=a*h
                total+=max(abs(low),abs(high))
            largest=max(largest,total)
        return Q(largest,ident)
    e=Q(1,2**166)
    inv=Q(max(sum(map(abs,r)) for r in pre["numerators"]),2**pre["depth"])
    # v/[1+(1-v)d] = v + v(v-1)d + v(1-v)^2*d^2/[1+(1-v)d].
    residual=norm(matrix,True)+e*norm(derivative,False)+46*inv*e*e/(1-e)
    assert residual < Q(2492,10**6)
    b=[Q(sum(x*y for x,y in zip(r,rhs)),den*2**pre["depth"]) for r in pre["numerators"]]
    floor=min(b)-max(map(abs,b))*residual/(1-residual)
    assert floor > Q(1,2**12)
    print("All rational scales in [1-2^-166,1+2^-166] fit the same 42 and four targets exactly with weights > 2^-12.")


if __name__=="__main__":
    main()

