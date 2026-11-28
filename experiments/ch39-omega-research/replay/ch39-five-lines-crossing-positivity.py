"""五个单位置例外内，牛顿目标与指定能标交汇强制全跑动区正性。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def verify():
    edges,rows,heat,_ = load("two-copy-gram-cutoff-two-exceptions")["model_inputs"]()
    cmin = arb(".001449")
    pi = arb.pi()
    slope = sum(a for _,a,_ in rows)/(2*pi)
    loss = sum((-w*a).max(0)*cmin*b*(1-heat(2,z/(cmin*b)))
               for _,b,w in edges for _,a,z in rows)/(2*pi)
    assert slope < 0
    for endpoint in (0,1):
        terms = [slope*line[2]*line[endpoint] for line in edges]
        improvements = [(slope*w*a).max(slope*w*b)-v
                        for (a,b,w),v in zip(edges,terms,strict=True)]
        top = sorted((v.upper() for v in improvements),reverse=True)[:5]
        bound = cmin*(sum(terms)+sum(top)).max(0)+loss
        assert bound < arb(".999446")
    model,_,_,initial,_ = load("chiral-minimal-neutral-gravity")["chiral"]["source"]["mean"]["inputs"]()
    fields = load("signed-gravity-unification")["gauge_weights"](model)
    L = (arb(10**16)/model["mz"]).log()
    bare = [a-b*L/(2*pi) for a,b in zip(initial,model["beta"],strict=True)]
    b1,b3 = [-sum(b[j] for b in fields)/(4*pi) for j in (0,2)]
    theta = b3/(b3-b1)
    q = [(theta*b[0]+(1-theta)*b[2])/(4*pi) for b in fields]
    const,_,_,heavy = load("upward-crossing-cutoff-inverse-moment")["inverse_moment_terms"](
        ROOT,model,initial,rows)
    z = [m for _,m in heavy]
    h = [-d/(4*pi) for d,_ in heavy]
    D = theta*bare[0]+(1-theta)*bare[2]+sum(v*m.log() for v,m in zip(q,z,strict=True))
    lam = -sum(v*m for v,m in zip(q,z,strict=True))/sum(v*m for v,m in zip(h,z,strict=True))
    p = [-a-lam*b for a,b in zip(q,h,strict=True)]
    offset = D+lam*(L-const)
    # 由 lambda 定义保证 sum(p*z)=0。先在原点附近用三阶余项给出正导数。
    small = -sum(v*m*m for v,m in zip(p,z,strict=True))/2
    small -= 200000*sum(abs(v)*m**3 for v,m in zip(p,z,strict=True))/6
    assert small > arb("1.0728e-14")
    assert all(1/(cmin*a) < 17000000 for a,b,w in edges)
    assert all(v > 0 or v < 0 for v in p)
    for left in range(200000,17000000,20000):
        lower = arb(0)
        for v,m in zip(p,z,strict=True):
            at = left+20000 if v > 0 else left
            lower += -v*(-m*at).expm1()/at
        assert lower > arb("2.7594e-9")
    # P(0)=0 且 P 单调，负谱线的最大损失来自其各自的最小可能节点。
    negative = arb(0)
    for a,_,w in edges:
        if w < 0:
            u = 1/(cmin*a)
            value = sum(v*(heat(1,m*u)+arb.const_euler()+(m*u).log())
                        for v,m in zip(p,z,strict=True))
            negative -= w*value
        else:
            assert w >= 0
    common = offset-negative
    assert model["beta"][2] == -7 and model["beta"][1] > -7 and model["beta"][0] > 0
    floor = common-7*L/(2*pi)
    assert common > arb("36.23729") and floor > arb(".22067")
    print("五例外内若匹配牛顿且在指定能标交汇，共同逆耦合与全程逆耦合下界",
          common.str(12),floor.str(12))


if __name__ == "__main__":
    verify()
