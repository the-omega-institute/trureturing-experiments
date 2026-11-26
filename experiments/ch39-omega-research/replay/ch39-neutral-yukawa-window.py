"""Bound the one-loop neutral gauge contribution to the Yukawa flow below M3."""
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
study = run_path(str(root/"ch39-top-threshold-resolution.py"))
q, retuned, empirical = study["q"], study["retuned"], study["empirical"]
comp, adjoint = study["comp"], study["adjoint"]


def charges(logs):
    pole, offsets, (alpha,beta,family), _ = retuned["current"]["light_current"](logs)
    return pole,offsets,*charge_bounds(alpha,beta,family)


def charge_bounds(alpha,beta,family):
    couplings = {name:[(d*alpha+r*beta+sign*x*family)**2 for x in (1,-2)]
                 for name,(d,r,sign) in retuned["current"]["weyl_charges"].items()}
    pairs = [(couplings[a],couplings[b]) for a,b in (("Q","uc"),("Q","dc"),("L","ec"))]
    low = [min(left)+min(right) for left,right in pairs]
    high = [max(left)+max(right) for left,right in pairs]
    trace = [sum(2*s[0]+s[1] for s in pair) for pair in pairs]
    return low,high,trace


@lru_cache(None)
def before(index,upper):
    bands = study["tables"]["table_bands"]()
    word = empirical["endpoint"]
    state = [q(Fraction(v,word[0])**2) for v in word[1:]]
    state += [q(row[index if j==2 else (-1 if upper else 0)][0]) for j,row in enumerate(bands)]
    for _ in range(160):
        state = comp["advance"](state,q(23)/160,upper,flow=adjoint["flow"](False,False))
    return state


def endpoint(index,upper,logs,parameters,time,*,check_speeds=False):
    pole, offsets, low, high, _ = parameters
    h3,h8,cal = map(q,logs)
    onset = cal+pole.log()/2
    state = before(index,upper)
    sm = adjoint["flow"](False,False)
    segments = ((q(23),onset,10,sm),
                (onset,h3,8,{**sm,"NEUTRAL_YUKAWA":low if upper else high}),
                (h3,h8,30,adjoint["flow"](True,False)),
                (h8,time,8,adjoint["flow"](True,True)))
    for a,b,n,flow in segments:
        assert b>a
        for _ in range(n):
            state,tube = comp["advance"](state,(b-a)/n,upper,flow=flow,need_tube=True)
            if check_speeds:
                assert all(abs(v)<2 for v in tube)
    target = comp["model"]["target"](11431,time)
    return [1/(r*v)-d-target for r,v,d in zip(adjoint["ratios"],state[:3],offsets,strict=True)]


def verify():
    logs = retuned["logs"]
    h3,_,cal = map(q,logs)
    params = charges(logs)
    pole,_,low,high,trace = params
    tau = 3*(h3-cal-pole.log()/2)/comp["loop"]
    for value,reference in zip(trace,(".005723122086",".007429651621",".007257433705"),strict=True):
        assert abs(tau*value-q(reference)) < q("1e-12")
    print("direct negative log determinants",*(v.str(12) for v in (tau*x for x in trace)))
    print("top singular-value loss limits",(tau*low[0]).str(12),(tau*high[0]).str(12))
    verify_classification(logs,params,q(".00000007"),q("3.41e-7"))


def verify_classification(logs,params,trial_shift,allowance,maximum_log_rates=None):
    h3,_,cal = map(q,logs)
    pole,_,low,high,_ = params
    for index,offset in enumerate((".00000266",".00000212")):
        time = cal+q(offset)
        lo,hi = [endpoint(index,upper,logs,params,time) for upper in (False,True)]
        assert hi[1] < -arb(2)**-21 and lo[1]-hi[0] > arb(2)**-20
        for upper in (False,True):
            endpoint(index,upper,logs,params,q(33),check_speeds=True)
        print("excluded group",index,(-hi[1]).lower().str(10),(lo[1]-hi[0]).lower().str(10))
    word = empirical["endpoint"]
    tables = study["tables"]["table_bands"]()
    inputs = []
    for upper in (False,True):
        inputs.append([q(abs(Fraction(v,word[0]))) for v in word[1:]]
                      +[q(row[2 if i==2 else (-1 if upper else 0)][0]).sqrt()
                        for i,row in enumerate(tables)])
    onset = cal+pole.log()/2
    trial = cal+trial_shift
    sm = adjoint["flow"](False,False)
    h8 = q(logs[1])
    unit_box = arb(0).union(arb(2))
    for flow in (sm,adjoint["flow"](True,False),adjoint["flow"](True,True)):
        speeds = [-2*(beta+sum(b*unit_box-d*3*unit_box for b,d in zip(bs,ds))/comp["loop"])
                  /(ratio*comp["loop"]) for beta,bs,ds,ratio in
                  zip(flow["BETA"],flow["B"],flow["D"],adjoint["ratios"],strict=True)]
        weak_rate = params[1][1]/(h3-onset) if maximum_log_rates is None else maximum_log_rates[1]
        assert speeds[0]<0 and speeds[1]>weak_rate>0
    segments = ((q(0),onset,160,sm),(onset,h3,8,{**sm,"NEUTRAL_YUKAWA_BOUNDS":(high,low)}),
                (h3,h8,30,adjoint["flow"](True,False)),(h8,trial,8,adjoint["flow"](True,True)),
                (trial,q(33),64,adjoint["flow"](True,True)))
    empirical["triplet"]["verify_ground"](inputs[0],segments,3,trial,allowance,q(".136"),
        upper_values=inputs[1],logarithm=(onset,h3,params[1]),stationary_readout=True,
        maximum_log_rates=maximum_log_rates)


if __name__ == "__main__":
    verify()
