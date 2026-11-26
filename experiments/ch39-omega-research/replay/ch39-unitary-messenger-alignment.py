"""缩放酉矩阵中介的共同酉方向与中心复场的精确能量比较。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

scalar = run_path(str(Path(__file__).with_name("ch39-free-complex-messengers.py")))


def verify():
    scalar["verify"]()
    ctx.prec = 160
    radius = arb(".9999").union(arb("1.0001"))
    h = radius**2
    energy = arb("3e-18")*(arb(25000)/42)**4
    budget = arb("4e-19")
    max_q = 60000
    closure_error = (budget/(2*energy)).sqrt()
    sin_chi = (budget/(energy*h.lower()**2)).sqrt()
    assert sin_chi.asin() < arb("2e-6")
    rho = (1-sin_chi**2).sqrt().union(arb(1))
    cosine = ((radius+arb(0,closure_error.upper()))**2-2*radius**2)/(2*h*rho)
    beta = cosine.acos()
    assert 2 < beta and beta < arb("2.2")
    # 链端项迫使 cos(delta)>0.99，并使费米子两分支都处在严格凸角区间。
    cos_delta = (1+(h*rho)**2-max_q*budget/energy)/(2*h*rho)
    assert cos_delta > arb(".99") and arb(".99").acos()+sin_chi.asin()<arb.pi()/2
    print("Relative unitary eigenphase half split below",sin_chi.asin().upper().str(10),
          "before comparison with the scalar minimum.")


if __name__ == "__main__":
    verify()
