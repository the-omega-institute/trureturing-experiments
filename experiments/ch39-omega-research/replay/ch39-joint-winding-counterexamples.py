"""同时允许编号和端口次数变化时，寻找比 11431 见证更便宜且能量更低的配置。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

winding=run_path(str(Path(__file__).with_name("ch39-first-arithmetic-plaquette-winding.py")))


def verify():
    ctx.prec=128
    spacing=2*arb.pi()/3
    angle=winding["angle"]
    target=abs(23226*angle(11431)-spacing)
    first_any, first_same = None,None
    for power in range(1,580):
        _,pairs=winding["near_resonances"](power,target.upper())
        for n,m in sorted(pairs):
            distance=abs(power*angle(n)-m*spacing)
            assert distance < target
            record=(power,n,m,distance)
            if first_any is None:
                first_any=record
            if first_same is None and m%3==1:
                first_same=record
    assert first_any[:3]==(515,40,2) and first_same[:3]==(579,124,1)

    # 在旧低能编号窗内，本预算只覆盖第一个正向共振。
    assert 23226*angle(10815) < 2*spacing-target
    lower=int(((spacing-target)/angle(10815)).floor().unique_fmpz())
    first_box=None
    for power in range(lower,22141):
        lo=1/(2*((spacing+target)/(2*power)).tan())
        hi=1/(2*((spacing-target)/(2*power)).tan())
        left=max(10815,winding["count"](lo.lower())+1)
        right=min(12092,winding["count"](hi.upper()))
        if left<=right:
            assert first_box is None and left==right
            first_box=(power,left,1,abs(power*angle(left)-spacing))
    assert first_box[:3]==(22140,10816,1) and first_box[-1]<target
    for record in (first_any,first_same,first_box):
        power,n,m,distance=record
        gain=2*(distance.cos()-target.cos())
        assert gain>0
        print("power, label, band",power,n,m,"distance",distance.str(11),"energy gain / g",gain.str(10),flush=True)


if __name__=="__main__":
    verify()
