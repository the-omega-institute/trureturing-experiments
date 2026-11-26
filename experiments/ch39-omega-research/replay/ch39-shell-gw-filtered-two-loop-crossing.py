"""Linearized two-loop logarithmic matching of the locally filtered shell model."""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat

folder = Path(__file__).resolve().parent
flow = run_path(str(folder/"ch39-shell-gw-first-heavy-flow.py"))
q = lambda n,d=1: arb(n)/d
base = [[q(199,50),q(27,10),q(44,5)], [q(9,10),q(35,6),12],[q(11,10),q(9,2),-26]]
extra = [[q(9,50),q(9,10),0],[q(3,10),q(13,6),0],[0,0,0]]
yukawa = [[q(17,10),q(1,2),q(3,2)],[q(3,2),q(3,2),q(1,2)],[2,2,0]]
db = [q(1,10),q(1,6),0]


def gauge_integrals(Q,m):
    G = list(flow["start"])
    result = [q(0) for _ in range(3)]
    events = [q(1000),*m,Q]
    for k in range(len(events)-1):
        dt = (events[k+1]/events[k]).log()
        vals = []
        for j,b in enumerate(flow["beta"]):
            slope = 2*(b+k*db[j])/flow["loop"]
            if j==1 and k==19:
                vals.append(dt*G[j])
            else:
                vals.append(-(1-slope*G[j]*dt).log()/slope)
            G[j] = G[j]/(1-slope*G[j]*dt)
        for i in range(3):
            result[i] += sum((base[i][j]+k*extra[i][j])*v for j,v in enumerate(vals))
    return result


def low_traces(end):
    unit = q(1,1024)
    N = int(float((end/unit).mid()))
    assert N*unit < end < (N+1)*unit
    data = flow["inputs"]()
    total = [q(0) for _ in range(3)]
    params = (flow["start"],flow["beta"])
    pad = arb(0,".02")
    box = flow["acb_mat"]([[flow["acb"](pad,pad) for j in range(3)] for i in range(3)])
    for k in range(N+1):
        x = k*unit
        h = unit if k<N else end-x
        span = q(0).union(h)
        slopes = flow["jets"]([v+box for v in data],x+span,1,params)[1]
        assert all(abs(h*z)<q(".0199") for M in slopes for z in M.entries())
        for i,A in enumerate(data):
            band = A+span*slopes[i]
            lo = sum(abs(z).abs_lower()**2 for z in band.entries())
            hi = sum(abs(z).abs_upper()**2 for z in band.entries())
            total[i] += h*lo.union(hi)*flow["loop"]
        data = flow["advance"](data,x,h,params)
    return total, sum(flow["fit"]["hs"](M,M).real for M in data)


def check(n):
    Q,m = flow["masses"](n)
    gauge = gauge_integrals(Q,m)
    quarks,light = low_traces((m[0]/1000).log()/flow["loop"])
    F,D = flow["fit"]["readouts"](n)
    early = (m[2]/m[0]).log()/flow["loop"]
    late = (Q/m[2]).log()/flow["loop"]
    # Integrate the scalar majorant up to Q; the two first extra columns vanish.
    a = light/(1-10*light*early)
    heavy_sum = q("7.9")**2*(q(81,977))**2*sum(flow["fit"]["hs"](M,M).real for M in (F,D))
    z = a+heavy_sum
    assert 10*early*light < 1 and 10*late*z < 1
    high = -flow["loop"]/10*((1-10*early*light)*(1-10*late*z)).log()
    base_delta = [(sum(c*y for c,y in zip(row,quarks))-g)/(32*arb.pi()**3)
                  for row,g in zip(yukawa,gauge)]
    # Keep the correlated high-scale simplex instead of three box intervals.
    derivative = sum(1-(125/mass)**2 for mass in m)
    slopesQ = [b+53*d for b,d in zip(flow["beta"],db)]
    # alpha_i^-1 difference changes by (-slopesQ * dlogQ + db * derivative * dlogM)/2pi.
    jac = arb_mat([[-(slopesQ[i]-slopesQ[2]), (db[i]-db[2])*derivative] for i in (0,1)])
    def corrections(y):
        residual = [v+sum(c*f for c,f in zip(row,y))/(32*arb.pi()**3) for row,v in zip(yukawa,base_delta)]
        tau,scale = jac.solve(-2*arb.pi()*arb_mat([[residual[i]-residual[2]] for i in (0,1)])).entries()
        alpha = residual[2]-slopesQ[2]*tau/(2*arb.pi())
        return tau,scale,alpha
    values = []
    for vertex in (-1,0,1,2):
        values.append(corrections([high if j==vertex else q(0) for j in range(3)]))
    def hull(vals):
        x = vals[0]
        for v in vals[1:]:
            x = x.union(v)
        return x
    tau,scale,alpha = [hull(z) for z in zip(*values)]
    assert q("-.111") < tau < q("-.095") and q(".331") < scale < q(".348")
    assert q(".221") < alpha < q(".240")
    corrected = Q*tau.exp()
    assert q("4.47e13") < corrected < q("4.56e13")
    print("zero",n,"two-loop delta log Q, log M, inverse alpha",
          *[(float(v.lower()),float(v.upper())) for v in (tau,scale,alpha)],
          "Q in GeV from exponentiated corrected log", (Q*tau.exp()).str(7))


if __name__=="__main__":
    for n in (43,195):
        check(n)
