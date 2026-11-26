"""单一伴随五维规范动能算子对两圈余差的匹配范围。"""
from pathlib import Path
from runpy import run_path

from flint import arb

light = run_path(str(Path(__file__).with_name("ch39-ten-light-two-loop-residual.py")))


def verify():
    read,mz,r,vector,projected = light["verify"]()
    def window(low, high, steps=64):
        return read((arb(low)/mz).log(),(arb(high)/mz).log(),steps,arb(3))[3]
    left, right = "4e15","9.5e15"
    assert window(left,left)[1] > 0 and window(right,right)[1] < 0
    # Exclude other scales by covering both exterior intervals.
    for a,b,sign in (("1e13",left,1),(right,"1e17",-1)):
        l,h = (arb(a)/mz).log(),(arb(b)/mz).log()
        for j in range(12):
            boxes = read(l+(h-l)*j/12,l+(h-l)*(j+1)/12,64,arb(3))[3]
            assert sign*boxes[1] > 0
    P,_,a5 = window(left,right)
    kinetic = P/projected
    scaled = kinetic*arb(60).sqrt()/(4*arb.pi())
    assert scaled < -arb(".047")
    assert a5 > 0
    # With canonical adjoint scalar kinetic term the off-diagonal generator
    # sees the average eigenvalue of its weak and colour endpoints.
    ax = a5+kinetic/2
    mu = arb(left).union(arb(right))
    mx = min(v.lower() for v in vector).union(min(v.upper() for v in vector)).sqrt()*mu/r
    v = (60*ax/(4*arb.pi())).sqrt()*mx/5
    assert ax > 0 and v > 0
    reduced_planck = arb("2.435e18")
    coefficient = scaled*reduced_planck/v
    assert arb("45") < a5 < arb("56")
    assert -arb(".13") < scaled < -arb(".047")
    assert -arb(150) < coefficient < -arb("28.7")
    print("Matching scale bounds, common inverse and c*v/M interval",left,right,a5.str(9),scaled.str(9))
    print("Single-adjoint mass calibration: vev GeV and reduced Planck coefficient",
          v.str(9),coefficient.str(9))
    print("Minimum |c_5| at the reduced Planck scale",(-coefficient).lower().str(10))


if __name__=="__main__":
    verify()
