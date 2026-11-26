"""分类四项 CKM 窗口的哪些子集会迫使零点 11431 出现。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

control = run_path(str(Path(__file__).with_name("ch39-arbitrary-zero-control.py")))
free = control["free"]
band_data = (
    # L 下上界，ub 相对模长的两端，各端的规范、端口、中心及范数乘子
    ("0", ".003", (".86", "-2.489", ".793", "4.10", "0"),
                   ("1.14", "-2.374", ".713", "4.08", "0")),
    (".003", ".0045", (".79", "-2.115", ".895", "3.06", "438"),
                      ("1.21", "-2.071", ".918", "3.16", "433")),
    (".0045", ".0065", (".69", "-.758", "3.790", "0", "2617"),
                       ("1.31", "-.681", "3.905", "0", "2693")),
    (".0065", ".00962", (".55", "1.001", "5.996", "0", "3716"),
                        ("1.45", "1.329", "6.444", "0", "3964")),
)


def verify_exclusions(data):
    values, projectors = data[-2:]
    tests = control["base_tests"](projectors)
    bonus, center_price, scores = {}, {}, []
    ua, up = arb(".000055"), arb(".006")
    tilts = (arb(".011"), arb(".002"))
    def add_score(name, weight, angle, tilt_sign, phase_ratio, low, high, multipliers):
        g, h, k, b = map(arb, multipliers)
        phase_coeffs = (phase_ratio, phase_ratio, -1, -phase_ratio)
        phase = sum((-acb(0, s)*m for s,m in zip(phase_coeffs,projectors)),acb_mat(2,2))
        matrix = (weight*abs(values[angle])/ua*(1-acb(0,tilt_sign*tilts[angle]))*projectors[angle]
                  + (1-weight)*phase/up)
        bonus[name], center_price[name] = b, k
        tests.append((name, matrix, g, h))
        scores.append((name, matrix, weight, angle, phase_ratio, arb(low), arb(high)))

    for sign in (-1, 1):
        for low, high, mult in (
            ("0",".003",("-3.002","1.443","1.932","0")),
            (".003",".005",("-1.744",".947","0","494")),
            (".005",".00962",("-.305","2.055","0","1531")),
        ):
            add_score(f"second {low} {sign}",arb(1),1,sign,arb(0),low,high,mult)
        for low, high, *ends in band_data:
            for rho, *mult in ends:
                add_score(f"first-phase {low} {rho} {sign}",arb(".34"),0,sign,arb(rho),low,high,mult)
    share = arb(".06")
    bounds = control["integrated_bounds"](data, tests, bonus, center_price, 1024, share)
    norm = bounds["norm"]
    assert norm < arb(".00962")
    axes = (1,5,29)
    lens = [((1-arb(axes[i])**2/867)*(1-arb(axes[j])**2/867)).sqrt()/abs(v)
            for (i,j),v in zip(control["entries"],values)]
    assert abs(values[2])*(1+norm*lens[2]) < arb(".006")
    modulus_floor = [arb(".224955")*(1-arb(".006")**2).sqrt(),
                     arb(".041775")*(1-arb(".006")**2).sqrt()]
    imag = {i:bounds[name+"1"].max(bounds[name+"-1"]) + norm**2*lens[i]/2
            for i,name in ((0,"us"),(1,"cb"),(3,"cs"))}
    assert all(imag[i]/(2*(1-norm*lens[i])) < tilts[i] for i in (0,1))
    for low, high, left, right in band_data:
        span = arb(high)*lens[2]
        assert 0 < arb(left[0]) < 1-span and 1+span < arb(right[0])
    q0 = (values[0]*values[1]*values[2].conjugate()*values[3].conjugate()).arg()
    assert 0 < q0 < arb.pi()/2 and sum((norm*f).asin() for f in lens) < 1
    worst = [arb(-100),arb(-100)]
    for name, matrix, weight, angle, rho, low, high in scores:
        upper = (weight*(abs(values[angle])-modulus_floor[angle])/ua + bounds[name]
                 -bonus[name]*low + control["nuclear"](matrix)*high**2/2)
        if angle == 0:
            errors = {}
            turns = {}
            for i in (0,1,3):
                radius = high*lens[i]
                a,b = 1-radius,1+radius
                if i == 0:
                    a,b = modulus_floor[0]/abs(values[0]),arb(".225065")/abs(values[0])
                turns[i] = (imag[i].min(radius)/a).asin()
                errors[i] = abs(1-a).max(abs(1-b))*turns[i] + b*turns[i]**3/6
            ref = q0-arb("1.141")
            max_turn = abs(ref) + sum(turns.values())
            assert max_turn < arb.pi()/2
            upper += (1-weight)/up*rho*(ref + sum(errors.values()) + max_turn**3/6)
        worst[angle] = worst[angle].max(upper.upper())
        assert upper < 0, (name,upper.str(8))
    print("Certified center share exceeds 6% for first angle plus phase, and for second alone")
    print("Score upper bounds", *(v.upper().str(12) for v in worst))


def verify_examples(data):
    point, scale, eps = data[:3]
    labels = (11415,11430,11432,11405)
    points = [point(n) for n in labels]
    blocks = [free["study"]["native"]["plane_response"](h) for h,_,_ in points]
    coeff = [v+arb(0,arb(2)**-28) for v in free["box"]["centers"]]
    for selected, counts in (
        ((0,2), (2644009509000,208801975260251,779542857116749,9011158114000)),
        ((2,3), (17285861547000,159524003534400,819576703449600,3613431469000)),
    ):
        assert sum(counts) == 10**15
        times = [arb(k)/10**15 for k in counts]
        gauge = sum(t*x for t,(_,x,_) in zip(times,points))
        port = gauge+scale*sum(t*c for t,(_,_,c) in zip(times,points))
        assert abs(gauge)<eps and abs(port)<eps
        steps = {n:free["ordered"]["evolution"](block,t,coeff) for n,block,t in zip(labels,blocks,times)}
        obs = free["ordered"]["readings"](steps,labels[::-1])
        widths = [arb(".000055")]*3+[arb(".006")]
        residuals = [(v-m)/w for v,m,w in zip(obs,free["flavor"]["shell"]["centers"],widths)]
        assert all(abs(residuals[i])<1 for i in selected)
        print("Center-free example",selected,"normalized residuals",*[v.str(6) for v in residuals])


if __name__ == "__main__":
    data = control["background"]()
    verify_examples(data)
    verify_exclusions(data)
