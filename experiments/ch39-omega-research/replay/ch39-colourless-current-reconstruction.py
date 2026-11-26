"""由小表示中的色单态电中性权重建判别电流，核验有限精度传播。"""
from pathlib import Path
from runpy import run_path

from flint import arb_mat

protocol = run_path(str(Path(__file__).with_name("ch39-rational-momentum-calibration.py")))
cal, mat, poles, q = (protocol[key] for key in ("cal","mat","poles","q"))
colour = arb_mat([[3],[0],[0]])
right = arb_mat([[0],[1],[0]])
higgs = arb_mat([[3],[-1],[2]])
family = arb_mat([[0],[0],[1]])


def condition(exchange, first, second, count):
    def value(source):
        return cal["form"](source,exchange,source)
    plus = (1+count)*value(first)+(count**2+count)*value(second)
    minus = count*value(first-second)
    direct = value(first+second*count)
    assert direct > 0 and (plus-minus).overlaps(direct)
    return (plus+minus)/direct


def check_neutral_weights():
    # (D, R, X), SU(4) antifundamental count, SU(2)_L and SU(2)_R labels.
    fields = (((3,0,0),1,1,0), ((0,1,0),0,1,1), ((3,-1,0),1,0,1),
              ((3,-1,2),1,0,1), ((0,-1,2),0,1,1), ((0,0,1),0,0,0))
    for (d,r,_), anti, left, right_spin in fields:
        assert d==3*anti and r in range(-right_spin,right_spin+1,2)
        assert (anti+left+right_spin)%2==0
        twice_t3 = -(d//3+r)
        assert twice_t3 in range(-left,left+1,2)
        assert (4**anti)*(left+1)*(right_spin+1)<=8


def verify():
    check_neutral_weights()
    time = poles["study"]["nominal"][2]
    for i, logs in enumerate(poles["retune"]["log_sets"]):
        (length,z,rho,eta,c,m),_,_ = poles["retune"]["predictions"](logs,time)
        e = protocol["anchor"]/length
        h = mat["flow"]["running"]["hypercharge_window"](time+protocol["anchor"].log()/2).mid()**2
        k = 1/rho
        origin = mat["massive_exchange"](e,z,c/rho,eta/rho)*k+mat["hypercharge"]*h
        p,a,_ = cal["shape"]["probe"](protocol["logs_x"][1],e,z,c/rho,eta/rho,m[0])
        tree = p*k+mat["hypercharge"]*h
        second = (cal["eye"]-tree*a).inv()*tree
        hyper_error = condition(origin,colour,right,9)
        assert q(".5")*hyper_error < 1
        print("case",i,"hypercharge amplification",hyper_error.str(9),flush=True)
        for n in (9,10):
            start = condition(origin,colour,higgs,3*n)
            end = condition(second,colour,higgs,3*n)
            # The nominal signals are each normalized by the family exchange
            # at the same fixed subtraction point.
            room = q("2e-5")/(end+(1+q("2e-5"))*start)
            assert room > q("5e-7")
            print("n",n,"amplification",start.str(8),end.str(8),
                  "common mixed precision (ppm)",(room*10**6).str(8),flush=True)
        print("hypercharge common precision (ppm)",(1/hyper_error).str(8),flush=True)


if __name__ == "__main__":
    verify()
