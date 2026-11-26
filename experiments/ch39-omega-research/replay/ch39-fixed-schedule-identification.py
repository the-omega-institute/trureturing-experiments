"""一个弧投影的味响应怎样鉴别时长固定、规范均值精确的时间表。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

free = run_path(str(Path(__file__).with_name("ch39-center-free-three-pulses.py")))
flavor = free["flavor"]
native = free["study"]["native"]
center = 11431


def verify():
    ctx.prec = 144
    projection = flavor["arcs"][32]
    complement = flavor["identity"]-projection
    reference = native["plane_response"](acb.zeta_zero(center).imag)
    assert acb.zeta_zero(center+1).imag > acb.zeta_zero(center).imag
    norms = []
    for first in range(1,center,256):
        count = min(256,center-first)
        zeros = acb.zeta_zeros(first,count)
        assert len(zeros)==count
        for n,z in enumerate(zeros,first):
            off_diagonal = projection*(native["plane_response"](z.imag)-reference)*complement
            norm = sum(abs(v)**2 for v in off_diagonal.entries()).sqrt()
            assert norm > 0, n
            norms.append((n,norm))
    nearest, gap = min(norms,key=lambda item:float(item[1].mid()))
    assert all(y>gap for n,y in norms if n!=nearest)
    assert nearest == 873 and gap > arb(".00010646464023357")
    print("Closest lower label",nearest,"off-diagonal norm",gap.str(17),
          "number of lower labels excluded",len(norms))


if __name__ == "__main__":
    verify()
