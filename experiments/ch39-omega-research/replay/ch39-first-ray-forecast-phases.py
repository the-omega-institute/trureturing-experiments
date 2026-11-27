"""只读实际零点与两钟相位，复核五个预告编号的输入条件。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


FORECAST = ((1, 1632060), (-4, 2086602), (2, 2040244), (-3, 1856372), (6, 2010251))


def verify(root):
    graph = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    centers = run_path(str(root/"ch39-first-ray-light-zeros.py"))["CASES"]
    ctx.prec = 192
    theta = graph["geo"]["theta"]
    for (q,n), (charge,_,a,b) in zip(FORECAST,centers,strict=True):
        assert q == charge
        h = acb.zeta_zero(n).imag
        errors = [abs(acb(0,q*(theta(d)-theta(244))*h).exp() - acb(0,arb(target)).exp())
                  for d,target in ((724,a),(804,b))]
        assert max(e.upper() for e in errors) < arb(".002")
        print("编号、荷、两钟至盒心的弦长",n,q,*(e.str(12) for e in errors))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
