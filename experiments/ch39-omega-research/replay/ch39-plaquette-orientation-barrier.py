"""三路径取向的连续势垒与有限链零取向切面的精确最低能。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
from sympy import cancel, diff, symbols

root = Path(__file__).parent
q = arb


def verify():
    # t=2cos(phi) 参数化在两个等长路径之间交换取向的连续变形。
    l,k,t,r,s = symbols("l k t r s", positive=True)
    c = k/(2*l+k*(1+t**2/2))
    small = 1-t*c*(t-1)/2
    large = 1+c*(t-1)
    energy = 2*l*(2*(r-1)**2+(s-1)**2)+k*(t*r-s)**2
    solution = {r:small,s:large}
    profile = 2*l*c*(t-1)**2
    assert all(cancel(diff(energy,v).subs(solution))==0 for v in (r,s))
    assert cancel(energy.subs(solution)-profile)==0
    assert cancel(diff(profile,t)-2*l*k*(t-1)*(4*l+k*(2+t))/(2*l+k*(1+t**2/2))**2)==0
    barrier = 2*l*k/(2*l+3*k)
    assert cancel(profile.subs(t,2)-barrier)==0
    low_radius = l/2+l*k/(4*(l+k))
    assert cancel(low_radius-barrier-l*(2*l+k)**2/(4*(l+k)*(2*l+3*k)))==0

    d,j,a = symbols("d j a", positive=True)
    defect = j/(d+3*j)
    # 实轴上两短一长时的有限计次约束极小值。
    quadratic = d*(2*(r-a)**2+(s-a)**2)+j*(2*r-s)**2
    point = {r:a*(1-defect),s:a*(1+defect)}
    assert all(cancel(diff(quadratic,v).subs(point))==0 for v in (r,s))
    assert cancel(quadratic.subs(point)-d*j*a**2/(d+3*j))==0
    ctx.prec = 128
    assert 3*q(".15")+q(2)*4/3 < 4
    assert (8/(8+q(".15")))*2/3 > q(1)/2
    theta = run_path(str(root/"ch39-first-arithmetic-plaquette-winding.py"))["angle"](11431)
    radial, closure = q(5), q(2)
    bare = 2*radial*closure/(2*radial+3*closure)
    for p in (23226,46453):
        cosine = (p*theta).cos()
        u, j = (2+cosine)/p, closure-cosine/p
        d = 2*radial+u
        a = 2*radial/d
        cut = 3*a*u+d*j*a**2/(d+3*j)
        above = bare+q(136)/(9*p)
        assert bare < cut < above
        assert q(136)/(9*p)/(q(8)/7) < q(".00057")
        print("链长",p,"取向切面最低能",cut.str(13),"连续路径最高能上界",above.str(13))


if __name__ == "__main__":
    verify()
