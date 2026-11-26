"""用两个动量、三项标量读数排除有隙的混合正矩阵谱。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

model = run_path(str(Path(__file__).with_name("ch39-two-dirac-current-witness.py")))
shape, mat, poles, q = (model[key] for key in ("shape","mat","poles","q"))
family = arb_mat([[0],[0],[1]])
probe = arb_mat([[28],[-9],[18]])


def form(left, matrix, right):
    return (left.transpose()*matrix*right)[0,0]


def verify():
    time = poles["study"]["nominal"][2]
    target = poles["study"]["comp"]["model"]["target"](11430,time)
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        (scale,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
        neutral = poles["reference_readings"](logs)[3]
        e, b, d = neutral**2/scale,c/rho,eta/rho
        h = mat["flow"]["running"]["hypercharge_window"](time+neutral.log()).mid()**2
        k = 1/(poles["offset"]+target/(15*scale))
        origin = mat["massive_exchange"](e,z,b,d)
        bare1, bare0 = (origin*coupling+mat["hypercharge"]*h for coupling in (k,1/rho))
        data = [shape["probe"](s,e,z,b,d,masses[0]) for s in (q(".5"),q(1))]
        eye = arb_mat([[1,0,0],[0,1,0],[0,0,1]])
        def dressed(coupling):
            result = []
            for p,a,_ in data:
                bare = p*coupling+mat["hypercharge"]*h
                result.append((eye-bare*a).inv()*bare)
            return result
        graphs, direct = dressed(k), dressed(1/rho)
        x1,x2 = (r[2] for r in data)
        gap = 4*masses[0]**2/e
        cap = (x1-1)/(x2-1)*(gap+x2)/(gap+x1)
        first,second = graphs
        limits = []
        for log_x,sign in ((q(".5"),-1),(q(1),1)):
            _,_,dest,delta,*_ = model["packet"](log_x,e,z,b,d,masses[0],1/rho,k,h)
            limits.append(delta+sign*q("2e-8")*dest)
        low,high = limits
        metric = second-first*cap
        assert low > 0 and high > 0 and mat["positive"](metric)
        inv_metric = metric.inv()
        u,v,w = first*family,second*family,second*probe
        axes = arb_mat([[v[i,0],w[i,0]] for i in range(3)])
        coefficients = ((".80187213","-.01114138"),(".80190694","-.01114168"),
                        (".80205761","-.01114293"))
        weights = arb_mat([[q(s)] for s in coefficients[case]])
        alpha,beta = weights.entries()
        residual = u-axes*weights
        square = form(residual,inv_metric,residual)
        noise = max(arb(0), square.upper()).sqrt()
        slack = (low/cap).sqrt()-abs(alpha)*high.sqrt()-noise
        assert slack > 0 and not beta.contains(0)
        increment = (slack/abs(beta))**2
        ref = form(probe,direct[1],probe)/form(probe,bare0,probe)
        smallest = (form(probe,second,probe)+increment)/form(probe,bare1,probe)
        difference = (smallest/ref-1)*10**6
        print("case",case,"rational mixed-current excess (ppm)",difference.lower().str(11),
              "witness coefficients",alpha.str(9),beta.str(9),flush=True)
        assert difference > (87,95,131)[case]


if __name__=="__main__":
    verify()
