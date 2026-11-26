"""Interval crossing of two-loop gauge trajectories with one split mirror copy."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
from flint import arb, ctx

sm = run_path(str(Path(__file__).with_name("ch39-two-loop-shell-separation.py")))
q, loop = sm["q"], sm["LOOP"]
db = list(map(q, (F(4,15),F(4,3),F(2,3))))
dB = sm["matrix"](["7/75 3/5 4/15", "1/5 49/3 4", "1/30 3/2 38/3"])


def threshold_crossing(events, subdivisions=256):
    start = [1/q(v)**2 for v in (".467453",".63811","1.0583")]
    yuk = sm["hull_inputs"](".00000615 .00311 .8616 .0000135 .000268 .01401 .0000028227 .00059465 .010101",F(1,100))
    y = [2*v.log() for v in yuk]
    correction = [arb(0) for _ in range(3)]
    answer = None
    step = q(F(1,subdivisions))
    for j in range(28*subdivisions):
        t = F(j,subdivisions)
        t1 = t+F(1,subdivisions)
        time = q(t).union(q(t1))
        durations = [q(max(0,t-high)).union(q(max(0,t1-low))) for low,high,_,_ in events]
        active = [q(int(t >= high)).union(q(int(t1 > low))) for low,high,_,_ in events]
        inv = [start[i]+correction[i]+arb(0,".0001")
               -2*(sm["BETA"][i]*time+sum(u*event[2][i] for u,event in zip(durations,events)))/loop
               for i in range(3)]
        assert all(a > q("1/2") for a in inv)
        gs = [1/a for a in inv]
        ys = [(v+arb(0,".004")).exp() for v in y]
        tr = [sum(ys[i:i+3]) for i in (0,3,6)]
        dx = [-2*(sum((sm["B"][i][k]+sum(u*event[3][i][k] for u,event in zip(active,events)))*gs[k] for k in range(3))
                 -sum(sm["D"][i][k]*tr[k] for k in range(3)))/loop**2 for i in range(3)]
        slope = [dx[i]-2*(sm["BETA"][i]+sum(u*event[2][i] for u,event in zip(active,events)))/loop for i in range(3)]
        assert slope[0]-slope[2] < q("-.10")
        dy = []
        for f,row in enumerate(sm["C"]):
            other = arb(0).union(tr[1-f]) if f<2 else arb(0)
            gauge = sum(c*v for c,v in zip(row,gs))
            dy.extend((3*(v-other)+2*(3*tr[0]+3*tr[1]+tr[2]-gauge))/loop for v in ys[3*f:3*f+3])
        assert all(abs(step*v) < q(".0001") for v in dx)
        assert all(abs(step*v) < q(".004") for v in dy)
        difference = inv[0]-inv[2]
        if j == 0:
            assert difference > 0
        if t1 == 28:
            assert difference < 0
        if difference.contains(0):
            current = (time,4*arb.pi()*(inv[1]-inv[2]),4*arb.pi()*inv[2])
            answer = current if answer is None else tuple(x.union(z) for x,z in zip(answer,current))
        correction = [x+step*v for x,v in zip(correction,dx)]
        y = [x+step*v for x,v in zip(y,dy)]
    assert answer is not None
    return answer


def crossing(low, high, subdivisions=256):
    return threshold_crossing([(low,high,db,dB)], subdivisions)


def verify():
    ctx.prec = 160
    lo,hi = F(2), F(9,4)
    for a,b in ((lo,lo),(hi,hi),(lo,hi)):
        t,difference,alpha = crossing(a,b)
        if a == b == lo:
            assert q("-.033") < difference < q("-.011")
        elif a == b == hi:
            assert q(".0032") < difference < q(".022")
        else:
            assert q("4.676e14") < 1000*t.exp() < q("4.827e14")
            assert q("38.69") < alpha < q("38.77")
        print(a,b,"t",t.lower().str(10),t.upper().str(10),
              "difference",difference.lower().str(10),difference.upper().str(10),
              "alpha",alpha.lower().str(10),alpha.upper().str(10))


if __name__ == "__main__":
    verify()
