"""原生十二生成元一圈响应中，五种 Weyl 多重态所需的最小谱分裂。"""
from fractions import Fraction as F

from flint import arb, ctx


def dot(weights, values):
    return sum(a*b for a, b in zip(weights, values))


def ball(value):
    value = F(value)
    return arb(value.numerator)/value.denominator


def verify():
    ctx.prec = 192
    # b is defined by d(alpha^-1)/d ln(mu) = -b/(2 pi).
    gauge = (0, F(-22, 3), -11)
    higgs = (F(1, 10), F(1, 6), 0)
    weyl = [(F(1, 15), 1, F(2, 3)), (F(8, 15), 0, F(1, 3)),
            (F(2, 15), 0, F(1, 3)), (F(1, 5), F(1, 3), 0), (F(2, 5), 0, 0)]
    beta = tuple(g + h + 3*sum(f) for g, h, f in zip(gauge, higgs, zip(*weyl)))
    assert beta == (F(41, 10), F(-19, 6), -7)
    c, j = (1, -3, 2), (5, -3, -2)
    assert dot(c, gauge) == dot(j, higgs) == 0
    assert dot(c, higgs) == F(-2, 5) and dot(j, gauge) == 44
    c_weights = [dot(c, row) for row in weyl]
    j_weights = [dot(j, row) for row in weyl]
    assert c_weights == [F(k, 5) for k in (-8, 6, 4, -4, 2)]
    assert j_weights == [-4, 2, 0, 0, 2]
    for weights, maximum in ((c_weights, F(12, 5)), (j_weights, 4)):
        assert sum(weights) == 0 and -sum(min(w, 0) for w in weights) == maximum
    mz, em, weak = arb("91.1876"), arb("127.930"), arb(".23129")
    initial = (3*(1-weak)*em/5, weak*em, 1/arb(".1180"))
    c0, j0 = dot(c, initial), dot(j, initial)
    lower, upper = [(arb(v)/mz).log() for v in ("1e13", "1e17")]
    tau = 2*arb.pi()
    required_higgs = -5*arb.pi()*c0
    assert arb("201.26604") < required_higgs < arb("201.26605")
    assert c0 + upper/(5*arb.pi()) < arb("-10.60831")

    def check_crossing(scale_log, gauge_log, spread, sm_baseline):
        # All three generations use Q,L at spread and u,d,e at zero.
        assert 0 < spread < scale_log and 0 <= gauge_log <= scale_log
        exposed = (spread, 0, 0, spread, 0)
        result = []
        for i in range(3):
            contribution = (sm_baseline*ball(beta[i])*scale_log + ball(gauge[i])*gauge_log
                            + ball(higgs[i])*scale_log
                            + 3*sum(ball(row[i])*f for row, f in zip(weyl, exposed)))
            result.append(initial[i] - contribution/tau)
        assert all((v - result[0]).contains(0) for v in result[1:]) and result[0] > 0
        return result[0]

    # For total responses in [0,L], C=0 requires D >=
    # (-5 tau C0 - 2 L)/36, with D the spread of generation means.
    # This and its ratio to L decrease with L. The upper endpoint realizes both.
    spread = (-5*tau*c0 - 2*upper)/36
    gauge_log = (tau*j0 + 12*spread)/44
    common = check_crossing(upper, gauge_log, spread, 0)
    assert arb("9.25750096") < spread < arb("9.25750097")
    assert arb(".26731811") < spread/upper < arb(".26731812")
    print("Total responses: minimum log spread, normalized spread, common inverse",
          spread.str(15), (spread/upper).str(15), common.str(15))
    print("Required equal-generation Higgs log; optimal gauge log",
          required_higgs.str(15), gauge_log.str(15))

    # A separate model adds one such set of responses to the full SM.
    # C gives (-5 tau C0 - 4 L)/36. J and nonnegative extra gauge exposure
    # give (44 L - tau J0)/12. Their maximum is minimized at their intersection.
    turn = (3*tau*j0 - 5*tau*c0)/136
    assert lower < turn < upper
    extra = (-5*tau*c0 - 4*turn)/36
    assert (extra - (44*turn - tau*j0)/12).contains(0)
    common = check_crossing(turn, arb(0), extra, 1)
    assert arb("7.93725573") < extra < arb("7.93725574")
    assert arb(".27184504") < extra/turn < arb(".27184505")
    print("Additional responses: crossing scale, minimum log spread, normalized spread, common inverse",
          (mz*turn.exp()).str(15), extra.str(15), (extra/turn).str(15), common.str(15))


if __name__ == "__main__":
    verify()
