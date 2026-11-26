"""等权半周配对在连续通过同一双钟干涉端口后的 CP 亏额。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_series, ctx

cp = run_path(str(Path(__file__).with_name("ch39-cp-ensemble-test.py")))


def read_pair(t, power, radius, spectral, polys):
    # 从 c 的二阶 Taylor 系数换元为 t 的系数，再作成功概率加权。
    total = arb_series([0],prec=3)
    up = arb_series([0],prec=3)
    current = arb_series([0],prec=3)
    for sign in (-1, 1):
        c = sign*radius*t
        source = cp["two"]["probabilities"](c,spectral)[0]
        cp_value = cp["cp_value"](c,polys,jet=True)
        weight = arb_series([1+sign*t,sign],prec=3)**power
        powers = [(sign*radius)**i for i in range(3)]
        up += weight*arb_series([source[i]*powers[i] for i in range(3)],prec=3)
        current += weight*arb_series([cp_value[i]*powers[i] for i in range(3)],prec=3)
        total += weight
    return up/total, current/total


def verify():
    ctx.prec, ctx.cap = 384, 3
    flavor = cp["flavor"]
    axes = [flavor["singlet_at"]({724:v,804:acb(1)})[1] for v in (acb(1),acb(0,1))]
    radius = sum(v*v for v in axes).sqrt()
    assert (axes[0]*(7*arb.pi()/64).cos()-abs(axes[1]))/radius > arb(".94")
    spectral = cp["two"]["spectral_data"]()
    polys = cp["cp_data"](spectral,flavor["reference_matrix"]())
    ref = flavor["singlet"](11661)[1]
    up0 = cp["two"]["probabilities"](ref,spectral)[0][0]
    j0 = cp["cp_value"](ref,polys)
    floors = (".25773602", ".00248489", ".0000141024", ".0000000792649")
    for power, floor in enumerate(floors, 1):
        def readings(t):
            return read_pair(t,power,radius,spectral,polys)
        pending = [(arb(".94"),arb(1))]
        while pending:
            a,b = pending.pop()
            u,j = readings(a.union(b))
            if u[1] < 0 and j[2]*u[1]-j[1]*u[2] > 0:
                continue
            assert b-a > arb("1e-9")
            middle = (a+b)/2
            pending.extend(((a,middle),(middle,b)))
        lo,hi = arb(".94"),arb(1)
        assert readings(lo)[0][0] > up0 > readings(hi)[0][0]
        for _ in range(100):
            mid = (lo+hi)/2
            value = readings(mid)[0][0]
            if value > up0:
                lo = mid
            elif value < up0:
                hi = mid
            else:
                break
        t = lo.union(hi)
        loss = 1-readings(t)[1][0]/j0
        assert loss > arb(floor)
        print("连续通过次数、命中 U₀ 的单对相位、平均 J 的相对亏额下界",
              power,t.str(14),loss.str(14))


if __name__ == "__main__":
    verify()
