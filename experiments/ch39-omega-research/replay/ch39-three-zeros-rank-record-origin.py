"""从三个实际零点用 Arb 认证六个钟相位与三个径向亏额的二进格。"""
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx


def readings(root, bits):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=384
    theta=lambda d:(1-arb(d)/1734).acos()
    omega=[theta(d)-theta(244) for d in (724,804)]
    scale=2**bits
    def cell(x):
        n=(scale*x).floor().unique_fmpz()
        assert n is not None and n<scale*x<n+1
        return int(n)
    result=[]
    for k in (-1,0,1):
        h=acb.zeta_zero(2_306_888_100_085+k).imag
        phases=[acb(0,w*h).exp() for w in omega]
        offset=shell["soft"]["point"]["study"]["offset"]
        delta=arb(23)/25*arb(2)**-20/(h*h+offset)
        result.append(tuple(cell(x) for z in phases for x in (z.real,z.imag))+(cell(delta),))
    return result


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    record=run_path(str(root/"ch39-three-zeros-integer-rank.py"))
    assert readings(root, record["BITS"])==record["RECORD"]
    print("三枚实际零点的十五个二进格全部吻合")
