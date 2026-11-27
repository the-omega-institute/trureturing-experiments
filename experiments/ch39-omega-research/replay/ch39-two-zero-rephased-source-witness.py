"""固定实源与可重定相实源在两个实际制备盒上的六质量和 CP 见证。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def verify():
    orbit = run_path(str(Path(__file__).with_name("ch39-two-zero-orbit-condition.py")))
    model,up,down = orbit["calibration"]()
    adjoint = model["adj"]
    gu,gd = [adjoint(a)*a for a in (up,down)]
    c = abs(down.det()/up.det())**(arb(2)/3)
    difference = [[4*gd[i,j].real-c*gu[i,j].real for j in range(3)] for i in range(3)]
    for k,lower in enumerate((".1698",".000337",".000121")):
        assert difference[k][k]>arb(lower)
        for i in range(k+1,3):
            for j in range(k+1,3):
                difference[i][j] -= difference[i][k]*difference[k][j]/difference[k][k]

    def read(k):
        yu,yd = arb("7.06539940386823251")*up*k,down*k
        (u,U),(d,D) = [model["eig"](adjoint(y)*y) for y in (yu,yd)]
        v = adjoint(U)*D
        phase = model["mixing_angles"](U,D)[-1]
        jarlskog = (v[0,1]*v[1,1].conjugate()*v[1,2]*v[0,2].conjugate()).imag
        return [x.sqrt() for x in u+d],phase,jarlskog

    fixed = acb_mat([[1,0,-3],[0,1,0],[0,0,arb(1)/2]])
    _,delta,_ = read(fixed)
    assert arb("1.1392") < delta < arb("1.1393")

    rows = [[arb(t) for t in row.split()] for row in (
        "2.880714335494008404 0 0",
        "1.756733222323356848 .002220899690655642 0",
        "1.524463172295672026 -.001134397590249137 .001242806240043137",
    )]
    phases = [acb(0,arb(t)).exp() for t in ("0","-2.785146433349310453","2.748101856546342320")]
    rephased = acb_mat([[z*x for x in row] for row,z in zip(rows,phases,strict=True)])
    masses,phase,jarlskog = read(rephased)
    assert all(abs(x-y)<w/8 for x,y,w in zip(masses,model["ref"][:6],model["tol"][:6],strict=True))
    assert arb(".000019") < jarlskog < arb(".000020")
    assert arb("2.1732") < phase < arb("2.1734")
    print("实换基的带号 CP 相位",delta.str(12))
    print("重定相实换基的六个 Yukawa", [x.str(12) for x in masses])
    print("重定相实换基的带号 CP 相位和 J",phase.str(12),jarlskog.str(12))


if __name__ == "__main__":
    verify()
