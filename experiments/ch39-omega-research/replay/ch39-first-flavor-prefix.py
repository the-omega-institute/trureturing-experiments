"""四项味读数的首个可行零点前缀及其显式双态混合。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    ctx.prec = 192
    flavor = run_path(str(Path(__file__).with_name("ch39-fixed-local-flavor.py")))
    coef = [arb(k)/2048 for k in flavor["codes"]]
    half = arb(".000055")
    upper_s = (arb(".22501")+half)**2
    lower_t = (arb(".04183")-half)**2
    selected = {}
    scores = []
    for n,zero in enumerate(acb.zeta_zeros(1,50),1):
        w = flavor["mixing"](flavor["native"]["plane_response"](zero.imag),coef)
        a,b,u,d = [abs(w[i,k])**2 for i,k in ((0,1),(1,2),(0,2),(1,0))]
        j = (w[0,0]*w[1,0].conjugate()*w[1,1]*w[0,1].conjugate()).imag
        if n<50:
            scores.append(42*(upper_s*(1-u)-a)+1969*(b-lower_t*(1-u)))
        if n in (6,50):
            selected[n] = [a,b,u,d,j]
    assert max(x.upper() for x in scores)<arb("-.016")
    weight = arb(19332)/10**6
    a,b,u,d,j = [weight*x+(1-weight)*y for x,y in zip(selected[6],selected[50])]
    ss,tt = a/(1-u),b/(1-u)
    real = (1-u)*(d-ss*(1-tt)-(1-ss)*tt*u)/2
    values = [ss.sqrt(),tt.sqrt(),u.sqrt(),acb(real,j).arg()]
    radii = [half]*3+[arb(".006")]
    assert all(abs(x-y)<h for x,y,h in zip(values,flavor["shell"]["centers"],radii))
    print("First feasible prefix 50; (6,50) weights",weight.str(7),(1-weight).str(7),
          "readouts",*(x.str(12) for x in values),"lower-prefix certificate",max(scores).str(9))


if __name__ == "__main__":
    verify()
