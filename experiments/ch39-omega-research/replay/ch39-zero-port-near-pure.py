"""三邻号近纯混合达到规范与端口最小误差的范围。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
box = run_path(str(root/"ch39-zero-port-flavor-box.py"))
study, sharp, joint = box["study"], box["sharp"], box["joint"]


def verify():
    ctx.prec = 240
    cost, strength, tangent = joint["port_cost"](include_tangent=True)
    scale = -strength/tangent
    heights = [acb.zeta_zero(n).imag for n in (11430,11431,11432)]
    xs = [study["center"]/(h*h+study["offset"])-1 for h in heights]
    xs[1] = arb(0)
    targets = [cost(xs[0])[0], arb(0), cost(xs[2])[0]]
    def left_weight(theta):
        l, r = [xs[i]+theta*scale*targets[i] for i in (0,2)]
        return -r/(l-r)
    interval = arb(0).union(1)
    shifted = [x+interval*scale*c for x,c in zip(xs,targets)]
    beta = (targets[0]-targets[2])/(shifted[0]-shifted[2])
    assert 1-interval*scale*beta > 0 and 1+(1-interval)*scale*beta > 0
    weight = left_weight(interval)
    assert 0 < weight < 1
    inputs = [study["native"]["plane_response"](h) for h in heights]
    delta = arb(2)**-28
    centers = box["centers"]
    details = [box["tight_data"](t,centers,[delta]*4) for t in inputs]
    obs, derivatives = zip(*details)
    threshold = arb(".041775")**2
    def second(o):
        return o[1]-threshold*(1-o[2])
    base = second(obs[1])
    end = weight*second(obs[0])+(1-weight)*second(obs[2])
    limit = base/(base-end)
    assert 0 < limit < 1
    signs = []
    for j in range(4):
        def dm(v): return v[1]+threshold*v[2]
        dc = dm(derivatives[1][j])
        dp = weight*dm(derivatives[0][j])+(1-weight)*dm(derivatives[2][j])
        gradient = (base*dp-end*dc)/(base-end)**2
        assert 0 not in gradient
        signs.append(1 if gradient>0 else -1)
    assert signs == [1, -1, -1, -1]
    best = [c+sign*delta for c, sign in zip(centers,signs)]
    obs = [box["flavor_data"](t,best)[0] for t in inputs]
    pair = [weight*a+(1-weight)*b for a,b in zip(obs[0],obs[2])]
    base, end = second(obs[1]), second(pair)
    fraction = base/(base-end)
    for i in range(32):
        r = fraction*arb(i).union(i+1)/32
        a,b,u,v,j = [c+r*(e-c) for c,e in zip(obs[1],pair)]
        assert arb(".224955")**2 < a/(1-u) < arb(".225065")**2
        assert b/(1-u) < arb(".041885")**2
        assert sharp["lo"] < u < sharp["hi"]
        real = ((1-u)*(u-v)-u*(a+b)+(1+u)/(1-u)*a*b)/2
        phase = acb(real,j).arg()
        assert arb("1.141") < phase < arb("1.153")
    for theta in (arb(0),arb(1)/2,arb(1)):
        w = left_weight(theta)
        end = w*second(obs[0])+(1-w)*second(obs[2])
        r = base/(base-end)
        mean_cost=w*targets[0]+(1-w)*targets[2]
        print("Near-pure outside limit",theta.str(3), r.str(22), "gauge and port margins per outside mass",
              (theta*scale*mean_cost).str(15),((1-theta)*strength*mean_cost).str(15))
    print("Flavor corner maximizing the interval", signs)


if __name__ == "__main__":
    verify()
