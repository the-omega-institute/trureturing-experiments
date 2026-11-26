"""在 1% 离散度与半数占据两种约束下，认证全谱混合的唯一取等解。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat

sharp = run_path(str(Path(__file__).with_name("ch39-sharp-gauge-flavor-variance.py")))
study = sharp["study"]
lo, hi = sharp["lo"], sharp["hi"]
cases = (
    ("dispersion", (11159,11431,11432,11559,11831), arb(".01")**2),
    ("share", (10762,11159,11431,11432,11831), arb(1)/2),
)


def entries(case, n, x):
    if case == "dispersion":
        return x*x, arb(int(n==11431))
    return arb(int(n==11431)), x*x


def candidate(case, labels, budget):
    rows = [sharp["row"](acb.zeta_zero(n)) for n in labels]
    gaps = [x for x,_ in rows]
    constrained, costs = zip(*(entries(case,n,x) for n,x in zip(labels,gaps)))
    probs, reference = sharp["attaining_distribution"](rows,[list(constrained)],[budget])
    value = sum(w*c for w,c in zip(probs,costs))
    print(case, "support", labels, "weights", *(w.str(16) for w in probs),
          "minimum", value.str(18), flush=True)
    if case=="share":
        print("minimum standard deviation at most half weight", value.sqrt().str(18), flush=True)
    return rows, constrained, costs, reference


def endpoints(case, rows, constrained, costs, reference):
    dets = []
    inequalities = []
    for U in (lo,hi):
        tests = sharp["phase_tests"](U,reference)
        m = arb_mat([[1,x,y,*tests(o)] for (x,o),y in zip(rows,constrained)])
        c = m.solve(arb_mat([[y] for y in costs]))
        dets.append(m.det())
        assert c[2,0]<0 and c[3,0]>0 and c[4,0]>0
        scale = -c[2,0] if case=="dispersion" else arb(1)
        def score(obs, tests=tests, c=c, scale=scale):
            a,p = tests(obs)
            return (c[3,0]*a+c[4,0]*p)/scale
        ceiling = sharp["compare_torus"](c[4,0]/scale,score)
        for n in (8499,16001):
            h = acb.zeta_zero(n).imag
            x = study["center"]/(h*h+study["offset"])-1
            slope = c[1,0]/scale
            assert (2*x>slope if n==8499 else 2*x<slope)
            assert x*x-slope*x > c[0,0]/scale+ceiling
        inequalities.append((c,scale,score))
    assert dets[0]*dets[1] > 0
    return inequalities


def neighboring_mixture():
    rows = [sharp["row"](acb.zeta_zero(n)) for n in (11430,11431,11432)]
    left, right = rows[0][0], rows[2][0]
    assert left>0>right
    share = -right/(left-right)
    center = rows[1][1]
    pair = [share*a+(1-share)*b for a,b in zip(rows[0][1],rows[2][1])]
    second_lower = arb(".041775")**2
    def second_margin(obs):
        return obs[1]-second_lower*(1-obs[2])
    start,end = second_margin(center),second_margin(pair)
    assert start>0>end
    cutoff = start/(start-end)
    for fraction in (arb(0),cutoff):
        a,b,u,v,j = [x+fraction*(y-x) for x,y in zip(center,pair)]
        assert arb(".224955")**2 < a/(1-u) < arb(".225065")**2
        assert lo < u < hi and b/(1-u) < arb(".041885")**2
    for i in range(16):
        t = cutoff*arb(i).union(i+1)/16
        a,b,u,v,j = [x+t*(y-x) for x,y in zip(center,pair)]
        real = ((1-u)*(u-v)-u*(a+b)+(1+u)/(1-u)*a*b)/2
        phase = acb(real,j).arg()
        assert arb("1.141") < phase < arb("1.153")
    print("Nearest pair, outside weight through",cutoff.str(16),"variance per outside weight",
          (-left*right).str(16),"99 percent relative deviation",(-left*right/100).sqrt().str(16))


def verify():
    neighboring_mixture()
    results = []
    for case,labels,budget in cases:
        results.append((case,labels,endpoints(case,*candidate(case,labels,budget))))
    minimum_slack = {}
    for start in range(8500,16001,256):
        count = min(256,16001-start)
        zeros = acb.zeta_zeros(start,count)
        assert len(zeros) == count
        for n,z in enumerate(zeros,start):
            x,obs = sharp["row"](z)
            for case,labels,duals in results:
                if n in labels:
                    continue
                constraint,cost = entries(case,n,x)
                for c,scale,score in duals:
                    slack = (cost-c[0,0]-c[1,0]*x-c[2,0]*constraint-scale*score(obs)).lower()
                    assert slack > 0
                    if case not in minimum_slack or slack < minimum_slack[case]:
                        minimum_slack[case] = slack
    for case,slack in minimum_slack.items():
        print(case,"smallest finite dual slack",slack.str(12))


if __name__ == "__main__":
    study["torus_bound"]()
    verify()
