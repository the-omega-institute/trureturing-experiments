"""右外壳荷六第 1 和 16 层固定在源低截止时，牛顿与弱差已不能同时匹配。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    n = 2_306_888_100_087
    source = run_path(str(root/'ch39-adjacent-ordered-charge-six-cutoffs.py'))['verify'](root)[2][1]
    mz, mp = arb('91.1876'), arb('1.22089e19')
    assert len(source) == 19 and arb('103.04') < source[0].sqrt()*mp < arb('103.05')
    native = run_path(str(root/'ch39-native-rh-ray-pvm.py'))
    _, _, cross, _ = next(x for x in native['shell_weights'](root,n) if x[0] == 6)
    gravity = run_path(str(root/'ch39-three-index-free-atoms.py'))
    _, _, fields = gravity['broken_fields'](root,(n,))
    thermal = run_path(str(root/'ch39-common-response-required-drawdown.py'))['field_bounds']
    k, _, error, _ = thermal(fields[0])
    ctx.prec = 320
    c = [v.real for v in cross[0]]
    p0,p1 = (mz/mp)**2, (arb('1.176397e14')/mp)**2
    weak = native['weak_series']()
    high = run_path(str(root/'ch39-full-neutral-free-cutoff.py'))['weak_column']
    global_range = arb('-.085').union(arb('-.04'))
    def inside(x, interval):
        return x.lower() > interval.lower() and x.upper() < interval.upper()
    # 对 1 <= x/p0 <= 1024 逐小段使用低核。其上用随 x 增大而减小的统一误差。
    for level in range(10):
        for i in range(16):
            t = arb(2)**level*(1+arb(i)/16).union(1+arb(i+1)/16)
            assert inside(weak(t),global_range)
    assert inside(high(1024*p0,p0,p1),global_range)
    fixed = weak(source[0]/p0)
    cut = arb('1e-14')
    weak18 = high(cut**2,p0,p1)
    free = [c[i] for i in range(19) if i not in (0,15,17)]
    weak_total = (c[0]+c[15])*fixed+c[17]*weak18+sum(v*global_range for v in free)
    assert weak_total > arb('.015'), weak_total
    # 若第十八层也低，K(x) = -k x + eps(x)，|eps| <= error。
    others = sum(max(arb(0),(-v).upper()) for i,v in enumerate(c) if i != 17)
    upper = k*(others+abs(c[17])*cut**2)+sum(abs(v) for v in c)*error
    assert upper < arb('.4'), upper
    print('第十八层振幅至少 1e-14 时弱差下界；否则牛顿项上界',
          weak_total.lower().str(14),upper.upper().str(14),flush=True)


if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)
