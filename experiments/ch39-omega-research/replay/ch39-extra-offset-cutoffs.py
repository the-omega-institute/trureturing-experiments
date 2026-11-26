"""Certify the first available zero index after a high-offset observation is added."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path
import json


def replay_lower(case, data, prefix):
    root=Path(__file__).resolve().parent
    cheb=run_path(str(root/"ch39-tick-separators-integers.py"))["chebyshev_upper"]
    f=1<<40; depth=case["depth"]; a=case["a"]; c=case["c"]
    levels=[Q(n,1<<288) for n in data]
    target=sum(t*(v-Q(int(t<0),2**288)) for t,v in zip(a,levels))/f
    cost=sum(map(abs,a))/Q(f)
    phase_norm=sum(map(abs,c))/Q(1<<depth)
    assert target>Q("0.00003")+cost*Q(1,2**272)+phase_norm*Q(1,2**20)
    pad=-((-phase_norm*(1<<16))//1)  # 2^112 times the 2^-96 phase error
    cap=Q(".000001")*2**112
    for h,*cosbins in prefix[:case["cut"]-1]:
        high=pad
        for j,t in zip(case["offsets"],a):
            h2=(h+int(t>=0))**2
            high-=((-t*h2*(1<<72))//(h2+(j*j<<800)))
        for t,x in zip(c,cosbins):
            lo,hi=cheb(x,5)
            high-=((-t*(hi if t>=0 else lo)*(1<<40))//(1<<depth))
        assert high<cap


def replay_upper(case, data, prefix):
    root=Path(__file__).resolve().parent
    word=run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    rates=run_path(str(root/"ch39-two-prime-height-law.py"))["rate_numerators"]()
    unit=1<<288
    counts=case["weights"]; ids=case["support"]
    assert min(counts)>0 and sum(counts)==unit and max(ids)==case["cut"]
    row=[(prefix[i-1][0]>>112,w) for i,w in zip(ids,counts)]
    for j,y in zip(case["offsets"],data):
        low=sum(w*h*h//(h*h+(j*unit)**2) for h,w in row)
        high=sum(-((-w*(h+1)**2)//((h+1)**2+(j*unit)**2)) for h,w in row)
        assert max(y+1-low,high-y)<2**16
    lo=[0]*4; hi=[0]*4
    for h,w in row:
        m,n=[-43235114615559591407865398948523*((h*r+2**311)//2**312) for r in rates]
        for j,k in enumerate((m,n,n-m,n+m)):
            a,b=word["power"](k)[0]
            lo[j]+=w*a
            hi[j]+=w*b
    assert max(map(abs,lo+hi))<word["DEN"]*2**16


def main():
    root=Path(__file__).resolve().parent
    cases=json.loads((root/"ch39-extra-offset-cutoffs.json").read_text())
    prefix=json.loads((root/"ch39-tick-separators.json").read_text())["prefix"]
    data=json.loads((root/"ch39-minimum-tick-finite-inputs.json").read_text())["levels"]
    extra=json.loads((root/"ch39-fixed-register-extra-target.json").read_text())
    run_path(str(root/"ch39-prime-clock-height-horizon.py"))["main"](5,15)
    for case in cases:
        assert case["cut"]<=len(prefix) and prefix[case["cut"]-1][0]+1<2**415
        values=[extra["floor"] if j==1682 else data[j] for j in case["offsets"]]
        replay_lower(case,values,prefix)
        replay_upper(case,values,prefix)
        print("Least highest index",case["cut"],"positive heights",len(case["support"]),"at 272-bit precision.")


if __name__=="__main__":
    main()

