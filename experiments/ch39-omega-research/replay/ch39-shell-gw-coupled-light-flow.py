"""Validated two-loop gauge and one-loop single-doublet matrix flow with a moving end scale."""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb

coeff = run_path(str(Path(__file__).with_name("ch39-shell-gw-filtered-two-loop-crossing.py")))
flow = coeff["flow"]
q = coeff["q"]
loop, eye = flow["loop"],flow["eye"]


def jets(matrices,gs,order):
    ys, gauges, mult, traces = [matrices],[gs],[],[]
    square, factors = [],[]
    for k in range(order):
        H = [sum((ys[j][i]*flow["adj"](ys[k-j][i]) for j in range(k+1)),acb_mat(3,3))
             for i in range(3)]
        tr = [sum(m[i,i].real for i in range(3)) for m in H]
        traces.append(tr)
        trace = 3*(tr[0]+tr[1])+tr[2]
        t = [trace-sum(c*a for c,a in zip(row,gauges[k])) for row in flow["C"]]
        diff = q(3,2)*(H[0]-H[1])
        mult.append([t[0]*eye+diff,t[1]*eye-diff,t[2]*eye+q(3,2)*H[2]])
        ys.append([sum((mult[j][i]*ys[k-j][i] for j in range(k+1)),acb_mat(3,3))/(k+1)
                   for i in range(3)])
        square.append([sum(gauges[j][i]*gauges[k-j][i] for j in range(k+1)) for i in range(3)])
        factors.append([(b if k==0 else 0)+(sum(c*a for c,a in zip(B,gauges[k]))
                           -sum(c*a for c,a in zip(D,tr)))/loop
                        for b,B,D in zip(flow["beta"],coeff["base"],coeff["yukawa"])])
        gauges.append([2*sum(square[j][i]*factors[k-j][i] for j in range(k+1))/(k+1)
                       for i in range(3)])
    return ys,gauges


def advance(data,g,h):
    # The last interval can end at any point specified by the mass parameter.
    pad = arb(0,".025")
    box = acb_mat([[acb(pad,pad) for j in range(3)] for i in range(3)])
    ys,gs = jets([y+box for y in data],[v+pad for v in g],1)
    assert all(abs(h*z)<q(".0249") for M in ys[1] for z in M.entries())
    assert all(abs(h*v)<q(".0249") for v in gs[1])
    span = q(0).union(h)
    tube = [y+v*span for y,v in zip(data,ys[1])]
    gtube = [y+v*span for y,v in zip(g,gs[1])]
    at = jets(data,g,4)
    entire = jets(tube,gtube,5)
    out = []
    for series,remainder in zip(at,entire):
        vals = remainder[5]
        for term in reversed(series):
            vals = [b+h*a for a,b in zip(vals,term)]
        out.append(vals)
    return *out, tube, gtube


def light(end):
    h = q(1,1024)
    N = int(float((end/h).lower()))
    x = N*h
    assert 0 <= end-x < 2*h
    data = flow["inputs"]()
    gs = flow["start"]
    for i in range(N):
        data,gs,_,_ = advance(data,gs,h)
    reference = flow["gauge"](x,(flow["start"],flow["beta"]))
    z = [4*arb.pi()*(1/g-1/v) for g,v in zip(gs,reference)]
    data,gs,tube,gtube = advance(data,gs,end-x)
    traces = [sum(abs(c).abs_upper()**2 for c in M.entries()).union(0) for M in tube]
    for j,(B,D) in enumerate(zip(coeff["base"],coeff["yukawa"])):
        correction = sum(d*t for d,t in zip(D,traces))-sum(b*v for b,v in zip(B,gtube))
        z[j] += (end-x)*loop*correction/(32*arb.pi()**3)
    return data,z


if __name__=="__main__":
    for n in (43,195):
        _,m = flow["masses"](n)
        s = q(".30").union(q(".39"))
        end = ((125**2+(m[0]**2-125**2)*(2*s).exp()).sqrt()/1000).log()/loop
        data,z = light(end)
        top = abs(data[0][2,2])
        assert q(".52") < top < q(".55") and all(abs(v)<2 for v in z)
        print("zero",n,"top entry across moved threshold",top.str(9),"two-loop gauge corrections",
              [v.str(6) for v in z])
