"""十二次的全类排除，含两个味部门及中性第11431枚零点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


PHASES = {
    -4: (
        (("0", "0"), ("0", "-1e-8"), ("1e-8", "1e-8")),
        ((".3590509", "-.412163"),),
    ),
    2: (
        (("0", "0"), ("1e-8", "0"), ("1e-8", "1e-8"), ("-1e-6", "-1e-6")),
        (("-.2266231", "-1.4679559"), ("-1.4367452", ".0474222"), (".3766593", ".9313407")),
    ),
}


def verify():
    replay = Path(__file__).resolve().parent
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    selector = run_path(str(replay/"ch39-light-spectrum-zero-selector.py"))
    flavor = run_path(str(replay/"ch39-global-light-action-degree.py"))
    ctx.prec = 384
    degree = 12
    vectors = flavor["flavor_vectors"](shell,pump,selector,PHASES,degree)
    block = list(pump["rotation_bases"](shell))[1]
    adj = pump["adjoint"]

    def energy(n):
        matrix = shell["laplacian"](0,acb.zeta_zero(n).imag)
        return flavor["lowest"](adj(block)*matrix*block)

    lower,middle,upper = [energy(n) for n in (11430,11431,11432)]
    assert lower < middle < upper
    a,b,c = [selector["chebyshev_moments"]([e,e],degree) for e in (lower,middle,upper)]
    prev = [(x-y)/(middle-lower) for x,y in zip(a,b,strict=True)]
    following = [(z-y)/(upper-middle) for z,y in zip(c,b,strict=True)]
    # 两个割线按正比例合成。组合取严格符号必有一个实际竞争者取同号。
    combined = [(x+y)/(upper-lower) for x,y in zip(prev,following,strict=True)]
    for v in (following,combined):
        length = sum(x*x for x in v).sqrt()
        assert length > 0
        vectors.append([x/length for x in v])
    radius = selector["interior_radius"](vectors)
    assert radius > arb("1.19e-11")
    print("十二维归一位移单形的原点面距下界",radius.str(15))


if __name__ == "__main__":
    verify()
