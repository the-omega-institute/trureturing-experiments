"""五项同精度时钟的最小最大误差和第一次共同二进排除。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


model = run_path(str(Path(__file__).with_name("ch39-global-fifth-time.py")))
shared = model["shared"]


def stationary_upper(readings, node, data):
    invariants, reference, lambdas, slopes = data
    moments = [v/200**j for j,v in enumerate(shared["clock"]["clock"]["cleared_moments"](readings,4))]
    middle = arb(node)
    radius = arb("1e-10")
    rule = model["stationary_rule"](model["endpoints"],moments,lambdas,slopes,
                                    middle-radius,middle+radius,-1)
    assert all(w > 0 for w in rule[1])
    assert 0 < rule[0][1] < rule[0][2]
    for j in range(1,5):
        assert model["integrate"](rule,lambda x:model["time_kernel"](x,j)).overlaps(readings[j])
    model["certify_maximum"](rule,invariants,reference)
    coefficients = model["certify_dual"](rule,lambdas,slopes,1)
    assert all((-1)**j*coefficients[j] > 0 for j in range(1,5))
    return model["integrate"](rule,lambda x:model["time_kernel"](x,5))


def verify():
    ctx.prec, ctx.cap = 768, 3
    data = shared["flavor_data"]()
    times = shared["clock"]["clock"]["times"]()
    for digits,side in (("2.55974507429e-10",-1), ("2.55974507430e-10",1)):
        error = arb(digits)
        readings = [1,*[times[j]+(-1)**j*error for j in range(1,5)]]
        upper = stationary_upper(readings,".00118060765594",data)
        assert side*(upper+error-times[5]) > 0 and upper < times[5]
    for bits,fraction,node in ((30,".396",".00118423436133"),(31,"1",".00112970122029")):
        low = [(v*2**bits).floor()/2**bits for v in times[:6]]
        corner = [low[j]+(arb(2)**-bits if j%2==0 else 0) for j in range(5)]
        readings = [1,*[times[j]+arb(fraction)*(corner[j]-times[j]) for j in range(1,5)]]
        upper = stationary_upper(readings,node,data)
        if bits == 30:
            assert all(low[j] < readings[j] < low[j]+arb(2)**-bits for j in range(1,5))
            assert low[5] < upper < low[5]+arb(2)**-bits
        else:
            assert (low[5]-upper)*2**bits > arb("45.57")
        print("共同位数和第五项相对该格左端的格数", bits, ((upper-low[5])*2**bits).str(12))
    print("最大绝对误差下确界包围", "2.55974507429e-10", "2.55974507430e-10")


if __name__ == "__main__":
    verify()
