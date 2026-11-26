"""将公共中心商的两种族流反例代入 1 TeV 输入盒的三支规范跑动。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
source = run_path(str(root/"ch39-published-split-adjoint.py"))
data, comp, q = (source[key] for key in ("data","comp","q"))
ratios = source["ratios"]


def fermion(flow, multiplicity, dimension, weak_casimir, colour_casimir, hypercharge):
    casimirs = [q(3)*q(hypercharge)**2/5,q(weak_casimir),q(colour_casimir)]
    # dim includes both the colour and weak degeneracies.
    indices = [multiplicity*dimension*c/n for c,n in zip(casimirs,(1,3,8),strict=True)]
    beta = [v+4*index/3 for v,index in zip(flow["BETA"],indices,strict=True)]
    matrix = [[row[j]+index*(4*casimirs[j]+q(20)*adjoint*(i==j)/3)
               for j in range(3)] for i,(row,index,adjoint) in enumerate(
                   zip(flow["B"],indices,(0,2,3),strict=True))]
    return {**flow,"BETA":beta,"B":matrix}


def field_flow(total, enabled):
    triplet, octet, a, b = enabled
    flow = source["flow"](triplet,octet)
    if a and total == 46:
        flow = fermion(flow,2,2,"3/4",0,"-1/2")
    if a and total == 60:
        flow = fermion(flow,1,1,0,0,1)
    if b:
        flow = fermion(flow,1 if total==46 else 2,6,0,"10/3","4/3")
        if total==46:
            flow = fermion(flow,1,3,0,"4/3","4/3")
    return flow


def thresholds(total):
    model = run_path(str(root/"ch39-quotient-dynkin-budget.py"))
    poles = model["poles"]
    bounds = []
    time = poles["study"]["nominal"][2]
    for i,logs in enumerate(poles["retune"]["log_sets"]):
        (scale,_,_,_,_,masses),_,_ = poles["retune"]["predictions"](logs,time)
        neutral = poles["reference_readings"](logs)[3]
        e = neutral**2/scale
        if total==46:
            r1 = q(model["mass46"][i])
            ma,mb = (e/r1).sqrt(),masses[0]
        else:
            r1,r2 = [q(v)+arb(0,q(".000001").upper()) for v in model["mass60"][i]]
            ma,mb = (e/r1).sqrt(),(e/r2).sqrt()
        cal = q(logs[2])
        assert q("26.81") < cal < q("26.82") < cal+(1000*masses[0]).log()
        bounds.append([cal+x.log() for x in (*masses,ma,mb)])
    merged = [a.union(b).union(c) for a,b,c in zip(*bounds,strict=True)]
    # Round switching times out to a 0.001 mesh. Each bound integrates
    # with a definite field content and sandwiches all three spectra.
    return [(Fraction(int((v*1000).lower().floor().unique_fmpz()),1000),
             Fraction(int((v*1000).upper().ceil().unique_fmpz()),1000)) for v in merged]


def verify_flow(total):
    bands = thresholds(total)
    times = sorted({Fraction(0),Fraction("26.81"),Fraction("26.82"),
                    *(value for band in bands for value in band)})
    states = [data["corner"](1,upper) for upper in (False,True)]
    minima, matching_minima = [], []
    for left,right in zip(times,times[1:]):
        mid = (left+right)/2
        flows = [field_flow(total,[mid >= band[int(upper)] for band in bands]) for upper in (False,True)]
        steps = 1+int(8*(right-left))
        dt = q((right-left)/steps)
        for _ in range(steps):
            tubes = []
            for i in range(2):
                states[i],tube = comp["advance"](states[i],dt,bool(i),need_tube=True,flow=flows[i])
                tubes.append(tube)
            # The lower corner gives larger gauge couplings.
            readings = [(1/(r*tubes[0][j])).union(1/(r*tubes[1][j]))
                        for j,r in enumerate(ratios)]
            assert all(value > 0 for value in readings)
            residual = (max(y.lower() for y in readings)-min(y.upper() for y in readings))/2
            minima.append(residual.lower())
            if left >= Fraction("26.81"):
                matching_minima.append(residual.lower())
    print("index",total,"threshold intervals",[(float(a),float(b)) for a,b in bands],
          "full-range residual lower",min(minima).str(10),
          "matching window lower",min(matching_minima).str(10),flush=True)
    global_limit, local_limit = ((".0470",".2071") if total==46 else (".0661",".2943"))
    assert min(minima) > q(global_limit) and min(matching_minima) > q(local_limit)


def all_minimal_spectra():
    # Compare with the smallest SM gauge trajectory and largest Yukawa
    # trajectory in the empirical box; all added matter is Yukawa-inert.
    state = data["corner"](1,True)
    flow = source["flow"](False,False)
    assert all(v >= 0 for v in flow["B"][0]+flow["D"][0])
    for _ in range(200):
        state = comp["advance"](state,q("26.81")/200,True,flow=flow)
    _, tube = comp["advance"](state,q(".01"),True,flow=flow,need_tube=True)
    sm_upper = (1/(ratios[0]*tube[0])).upper()
    target_lower = comp["model"]["target"](11430,q("26.82")).lower()
    decrement = q(36)/(10*arb.pi()**2*ratios[0])
    deficit = target_lower-sm_upper+decrement
    print("all index 46: SM upper",sm_upper.str(11),"hypercharge decrement lower",
          decrement.lower().str(11),"arithmetic target lower", target_lower.str(11),
          "fixed-scale deficit lower",deficit.lower().str(11),flush=True)
    assert deficit > q(".3527")


if __name__ == "__main__":
    all_minimal_spectra()
    for total in (46,60):
        verify_flow(total)
