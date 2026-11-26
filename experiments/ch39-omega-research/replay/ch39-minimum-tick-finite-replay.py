"""Replay the upper bounds on minimum tick with integers and rational arithmetic."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path
import json


def main():
    root=Path(__file__).resolve().parent
    source=json.loads((root/"ch39-minimum-tick-finite-inputs.json").read_text())
    w=run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    rates=run_path(str(root/"ch39-two-prime-height-law.py"))["rate_numerators"]
    den=1<<source["depth"]
    for key,bits in (("2",96),("8",272)):
        rows=source["cases"][key]["rows"]
        assert sum(r[-1] for r in rows)==den and min(r[-1] for r in rows)>0
        coeff=rates(int(key))
        lo=[0]*4; hi=[0]*4
        for _,h,weight in rows:
            u,v=[-43235114615559591407865398948523*((h*a+2**311)//2**312) for a in coeff]
            for i,k in enumerate((u,v,v-u,v+u)):
                left,right=w["power"](k)[0]
                lo[i]+=weight*left
                hi[i]+=weight*right
        phase=max(max(map(abs,lo)),max(map(abs,hi)))/Q(den*w["DEN"])
        error=Q(0)
        for j,f in enumerate(source["levels"]):
            # The height can be anywhere inside its certified input cell.
            low=sum(q*h*h//(h*h+(j*den)**2) for _,h,q in rows)
            high=sum(-((-q*(h+1)**2)//((h+1)**2+(j*den)**2)) for _,h,q in rows)
            error=max(error,Q(max(abs(low-f-1),abs(high-f)),den))
        assert max(error,phase)<Q(1,2**bits)
        print("tick",key,"heights",len(rows),"time error/tolerance",float(error*2**bits),"phase/tolerance",float(phase*2**bits))


if __name__=="__main__":
    main()

