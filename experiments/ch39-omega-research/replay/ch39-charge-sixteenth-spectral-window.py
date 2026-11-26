"""荷六零弱差及牛顿目标对第十六层截止和源高层投影概率的必要窗口。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def weak_cells(weak, mz, minimum, steps=128):
    # 从源最小截止开始逐倍频程覆盖。在目标附近细分，直到弱差界可与高端区间接续。
    cells = []
    for power in range(8):
        unit = minimum * 2**power
        for i in range(steps):
            q = (unit*(1+arb(i)/steps)).union(unit*(1+arb(i+1)/steps))
            cells.append((q, weak((q/mz)**2)))
    return cells


def verify(root):
    centered = 2_306_888_100_085
    trained = run_path(str(root / 'ch39-three-shell-two-space-readout.py'))['verify']
    table, _ = trained(root, return_data=True)
    native = run_path(str(root / 'ch39-native-rh-ray-pvm.py'))
    training_values = run_path(str(root / 'ch39-adjacent-pvm-transport.py'))['training_values']
    constants = run_path(str(root / 'ch39-common-response-required-drawdown.py'))['field_bounds']
    fields = run_path(str(root / 'ch39-three-index-free-atoms.py'))['broken_fields']
    high_weak = run_path(str(root / 'ch39-full-neutral-free-cutoff.py'))['weak_column']
    targets = [centered-2,centered+2]
    _, _, data = fields(root, targets)
    ctx.prec = 320
    weak = native['weak_series']()
    mz, mp = arb('91.1876'), arb('1.22089e19')
    p0, p1 = (mz/mp)**2, (arb('1.176397e14')/mp)**2
    for side, n, field in zip((-1,1), targets, data, strict=True):
        _, levels, _, _ = next(z for z in native['shell_weights'](root, centered+side) if z[0] == 6)
        source = training_values(levels, table)
        _, _, cross, _ = next(z for z in native['shell_weights'](root,n) if z[0] == 6)
        ctx.prec = 320
        lo = min(x.lower() for x in source).union(min(x.upper() for x in source))
        hi = max(x.upper() for x in source)
        bright = [x for x in source if x > arb('.001')]
        assert len(bright) == (4 if side == -1 else 2)
        # 其余层在源构造中取同一个低振幅。
        assert all(x > arb('.001') or x.overlaps(lo) for x in source)
        c = [z.real for z in cross[0]]
        assert all(w > 0 or w < 0 for w in c) and c[15] > 0 > c[17]
        k, _, ek, _ = constants(field)
        other_neg = sum(-x for j,x in enumerate(c) if j != 17 and x < 0)
        crit = ((1 - ek*sum(abs(x) for x in c))/k - other_neg*hi**2)/abs(c[17])
        assert crit > 0
        amp18 = crit.sqrt().lower()
        d18 = high_weak(amp18**2,p0,p1)
        minimum = (lo*mp).lower()
        cells = weak_cells(weak, mz, minimum)
        tail = high_weak((minimum*256/mp)**2,p0,p1)
        values = [d for _,d in cells] + [tail]
        entire = min(z.lower() for z in values).union(max(z.upper() for z in values))
        other = c[17]*d18 + sum(x*entire for j,x in enumerate(c) if j not in (15,17))
        assert not (other+c[15]*tail).contains(0)
        kept = [q for q,d in cells if (other+c[15]*d).contains(0)]
        assert kept
        qmin,qmax = min(z.lower() for z in kept),max(z.upper() for z in kept)
        dim_floor = max(z.upper() for z in source if z < arb('.001'))
        probmin = ((qmin/mp-dim_floor)/(hi-dim_floor)).lower()
        probmax = ((qmax/mp-lo)/(min(z.lower() for z in bright)-lo)).upper()
        print(n, '第18层振幅必要下限、第16层截止范围及来自源高空间的概率范围',
              amp18.str(12),qmin.str(12),qmax.str(12),probmin.str(12),probmax.str(12),flush=True)
        if side == -1:
            assert amp18 > arb(".2658958") and arb("92.615") < qmin < qmax < arb("170.983")
            assert arb("1.3786e-19") < probmin < probmax < arb("1.2574e-17")
        else:
            assert amp18 > arb(".2473094") and arb("104.651") < qmin < qmax < arb("156.173")
            assert arb("1.7454e-19") < probmin < probmax < arb("1.0131e-17")


if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)
