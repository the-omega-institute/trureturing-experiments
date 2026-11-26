"""Lower comparison for a running scalar partial-wave expectation in the X sector."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

TEN_ONE, TEN_SELF, TEN_CROSS = run_path(str(Path(__file__).with_name("ch39-two-loop-product-beta.py")))["target_inverse_coefficients"]()

from flint import arb, ctx

ROOT = Path(__file__).parent
MATCH = run_path(str(ROOT/"ch39-quartic-soft-adjoint-match.py"))
a = MATCH["a"]


def gauge_source():
    dim = 100
    ca,cb = F(18,5), F(99,20)
    def same(c, adj, count):
        return F(6,dim+1)*(dim*c*c/count+c*c-adj*c/4)
    aa,bb,ab = same(ca,5,24),same(cb,10,99),F(12,dim+1)*ca*cb
    return a(aa),a(bb),a(ab)


def crossing(states, species, initial):
    state = states[0]
    for point in states[1:]:
        state = [v.union(w) for v,w in zip(state,point)]
    g,h,x,_,_,y,*_ = state
    aa,bb,ab = gauge_source()
    k = 1/(16*arb.pi()**2)
    c5,c10,ybox = a("18/5"),a("99/20"),a("120/101")
    def rhs(g,h,x,y):
        z,w = 1/g,1/h
        casimir = a("54/5")*z+a("297/20")*w
        return [2*k*(a("10/3")-(a("6796/15")*z+594*w-45*(x+y))*k),
                2*k*(TEN_ONE+(TEN_SELF*w-TEN_CROSS*z+15*(x+y))*k),
                2*k*x*(13*x+10*y-casimir),2*k*y*(10*x+13*y-casimir)]
    step = a("1/4096")
    radius = (-step).union(step)
    lower = a(initial).lower()
    cap = 8*arb.pi()/202
    for j in range(12288):
        tube = [v+radius for v in (g,h,x,y)]
        slope = rhs(*tube)
        assert all(abs(ds)<1 for ds in slope)
        z,w = 1/tube[0],1/tube[1]
        strength = tube[2 if species == 'X' else 3]
        radial = lower+radius
        derivative = k*(416*(radial*radial)+(-12*(c5*z+c10*w)+12*strength)*radial
                        +aa*z*z+bb*w*w+ab*z*w-ybox*(strength*strength))
        speed = derivative.lower()
        assert abs(speed)<1
        lower = (lower+step*speed).lower()
        if lower>cap:
            t = (j+1)*step
            assert t < (a(".66") if species=="X" else a(3))
            print(species,"contact crossing by t",t.str(10),"scale ratio",t.exp().str(9))
            return t
        g,h,x,y = [v+step*ds for v,ds in zip((g,h,x,y),slope)]
    raise AssertionError("No certified crossing in the tested interval")


def verify():
    ctx.prec = 140
    points = list(MATCH["matched_points"]())
    states = [state for _,_,_,_,state in points]
    t = crossing(states, "X", "242/2424")
    assert all(v*t.exp() < a("1.15e15") for _,_,_,v,_ in points)
    t_y = crossing(states, "Y", 0)
    assert all(v*t_y.exp() < a("1.2e16") for _,_,_,v,_ in points)


if __name__ == "__main__":
    verify()
