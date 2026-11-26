"""邻号线性段结束后，11406 进入最优混合直至相位触及下边。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_poly, ctx

root = Path(__file__).parent
box = run_path(str(root/"ch39-zero-port-flavor-box.py"))
study, sharp, joint = box["study"], box["sharp"], box["joint"]


def verify():
    ctx.prec = 240
    cost, strength, tangent = joint["port_cost"](include_tangent=True)
    scale = -strength/tangent
    ids = (11430,11431,11432,11406)
    heights = [acb.zeta_zero(n).imag for n in ids]
    xs = [study["center"]/(h*h+study["offset"])-1 for h in heights]
    xs[1] = arb(0)
    costs = [cost(x)[0] for x in xs]
    costs[1] = arb(0)
    inputs = [study["native"]["plane_response"](h) for h in heights]
    radius = arb(2)**-28
    theta = arb(0).union(1)
    centers = box["centers"]
    obs, derivatives = zip(*(box["tight_data"](t,centers,[radius]*4) for t in inputs))
    T = arb(".041775")**2
    def second(o):
        return o[1]-T*(1-o[2])
    margins = list(map(second,obs))

    def geometry(t):
        v = [x+t*scale*c for x,c in zip(xs,costs)]
        w = -v[2]/(v[0]-v[2])
        k = (v[3]-v[2])/(v[0]-v[2])
        return v,w,k
    v,w,k = geometry(theta)
    def difference(y, k):
        return y[3]-k*y[0]-(1-k)*y[2]
    dc, dm = difference(costs,k), difference(margins,k)
    multiplier = dc/dm
    assert dc>0 and dm>0
    beta=(costs[0]-costs[2]-multiplier*(margins[0]-margins[2]))/(v[0]-v[2])
    intercept=costs[2]-beta*v[2]-multiplier*margins[2]
    assert 1-theta*scale*beta>0 and 1+(1-theta)*scale*beta>0
    # 在三点外部支撑上，给定外部权重 r 后的代价是 r 的仿射函数。
    examined = arb(".0278").union(arb(".22"))
    def numerator(r):
        return -(1-r)*margins[1]-r*(w*margins[0]+(1-w)*margins[2])
    corners = []
    for j in range(4):
        derivative = [p[j][1]+T*p[j][2] for p in derivatives]
        d_pair = w*derivative[0]+(1-w)*derivative[2]
        # 对整个 r 区间检查仿射表达式的两个端点。
        ends = []
        for r in (examined.lower(),examined.upper()):
            ends.append(-((1-r)*derivative[1]+r*d_pair)*dm-numerator(r)*difference(derivative,k))
        sign = ends[0].union(ends[1])
        assert 0 not in sign
        corners.append(-1 if sign>0 else 1)
    print("Cost minimizing taste signs",corners,flush=True)
    assert corners == [1,-1,-1,-1]
    corner = [c+s*radius for c,s in zip(centers,corners)]
    best = [box["flavor_data"](t,corner)[0] for t in inputs]
    margin = list(map(second,best))

    def branch(t):
        _, weight, projection = geometry(t)
        pair = weight*margin[0]+(1-weight)*margin[2]
        surplus = difference(margin,projection)
        start = margin[1]/(margin[1]-pair)
        r = arb_poly([0,1])
        new = ((margin[1]-pair)*r-margin[1])*(1/surplus)
        weights = [weight*r-projection*new,1-r,(1-weight)*r-(1-projection)*new,new]
        a,b,u,h,j = [sum(y*o[k] for y,o in zip(weights,best)) for k in range(5)]
        real = (1-u)**2*(u-h)-u*(1-u)*(a+b)+(1+u)*a*b
        imag = 2*(1-u)*j
        boundary = imag-arb("1.141").tan()*real
        upper = imag-arb("1.153").tan()*real
        assert boundary.derivative()(examined)<0 and boundary(arb(".22"))<0
        assert boundary(arb(".2"))>0 and upper(examined)<0
        assert real(examined)>0 and imag(examined)>0
        assert (a-arb(".224955")**2*(1-u))(start)>0
        assert (a-arb(".224955")**2*(1-u)).derivative()(examined)>0
        assert (a-arb(".225065")**2*(1-u))(examined)<0
        assert sharp["lo"]<u(examined)<sharp["hi"]
        assert all(y(start)>0 and y.derivative()(examined)>0 for y in (weights[0],weights[2]))
        assert new.derivative()(examined)>0
        left,right=arb(".2"),arb(".22")
        for _ in range(76):
            mid=(left+right)/2
            val=boundary(mid)
            if 0 in val:
                break
            if val>0:
                left=mid
            else:
                right=mid
        end=left.union(right)
        assert boundary(left)>0>boundary(right)
        return start,end,weights
    for t in (theta,arb(0),arb(1)/2,arb(1)):
        first,last,weights=branch(t)
        assert examined.lower()<first<last<examined.upper()
        print("Budget direction", t.str(3),"branch from",first.str(18),"to",last.str(18),flush=True)

    # 定标正半径控制的全谱支撑平面。尾部使用原径向轴给出的 |V_cb|² 上界。
    for n in (10699,12251):
        h = acb.zeta_zero(n).imag
        x = study["center"]/(h*h+study["offset"])-1
        c,dx=cost(x)
        derivative=dx/(1+theta*scale*dx)
        assert (derivative>beta if n==10699 else derivative<beta)
        assert c-beta*(x+theta*scale*c)>intercept+multiplier*arb(104)/867
    coeff_box = [c+arb(0,radius) for c in centers]
    floor = arb(1)
    for start in range(10700,12251,256):
        count = min(256,12251-start)
        zeros=acb.zeta_zeros(start,count)
        assert len(zeros)==count
        for n,z in enumerate(zeros,start):
            if n in ids:
                continue
            x=study["center"]/(z.imag**2+study["offset"])-1
            c=cost(x)[0]
            t=study["native"]["plane_response"](z.imag)
            m=second(box["observations"](box["flavor"]["mixing"](t,coeff_box)))
            slack=c-beta*(x+theta*scale*c)-intercept-multiplier*m
            assert slack>0,(n,str(slack))
            floor=min(floor,slack.lower())
    print("Full-spectrum minimum finite strict slack",floor.str(12))


if __name__ == "__main__":
    verify()
