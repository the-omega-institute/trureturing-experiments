"""随绕数非增的全局真空与十位参考味码。"""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


DIRECTIONS = ((1, 0), (0, 1), (1, -1))
BOXES = ((".0051", ".0145", "-.1424", "-.1350"), (".1350", ".1423", "-.0147", "-.0053"))
LIMIT = 8000


def geometry(root, exact):
    faces = exact["enumerate_faces"]()
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    lib = run_path(str(root / "ch39-mincut-global-winding.py"))
    tri = lib["terms"](shell)
    data = []
    for k in DIRECTIONS:
        expected = Counter(z for key, z in faces.values() if key == k)
        units = {z: exact["unit"](z) for z in expected}
        actual = Counter()
        for native, spoke in tri[k]:
            label, = [z for z, unit in units.items() if (unit - native).contains(0)]
            assert spoke in (0, 1)
            actual[label, spoke] += 1
        for z, n in expected.items():
            assert n == expected[exact["conjugate"](z)] == 3
            assert n == actual[z, 0] + actual[z, 1]
            data.append((k, z, units[z], actual[z, 1]))
    assert len(data) == 32 and sum(s for _, _, _, s in data) == 6
    return data


def code(c, bits):
    return int((2**bits * c).floor().unique_fmpz())


def physical_cover(root, distance, boxes):
    """荷二具有 π 周期。两种等码角差、两个折叠原像皆随后的能量比较覆盖。"""
    schur = run_path(str(root / "ch39-twentyseven-coefficient-window.py"))["signs"]
    inertia = run_path(str(root / "ch39-neutral-mean-phase-inertia.py"))
    band = run_path(str(root / "ch39-curvature-window-density-cap.py"))["eigen_band"]()
    lower, upper = arb("18.94088212"), arb("18.94414547")
    extra = arb("4.48e-9") - arb(23) / 25 * arb(2)**-20 / 14**2
    assert 0 < extra and lower < band[0] - extra/3 < band[1] + 4*extra < upper
    counts = [0, 0]
    for sign in (-1, 1):
        pending = [(-arb.pi()/2, arb.pi()/2, distance.lower(), distance.upper())]
        while pending:
            left, right, dl, dr = pending.pop()
            t, d = left.union(right), dl.union(dr)
            y = t - d
            if any(a < t < b and c < y < e for a, b, c, e in boxes):
                counts[0] += 1
            elif (schur(inertia, 2, lower, sign*t, sign*y)[0] >= 16
                  or schur(inertia, 2, upper, sign*t, sign*y)[1] < 16):
                counts[1] += 1
            elif right - left > dr - dl:
                mid = (left+right)/2
                pending.extend(((left,mid,dl,dr), (mid,right,dl,dr)))
            else:
                assert dr-dl > arb("1e-7"), (t,d,y)
                mid = (dl+dr)/2
                pending.extend(((left,right,dl,mid), (left,right,mid,dr)))
    print("十位角差闭包、相位框内与惯性排除末格数", distance.str(11), counts, flush=True)


def ten_bits(root, previous, data, flavor):
    ref = flavor["singlet"](11661)[1]
    a, b = [flavor["singlet_at"]({724:z,804:acb(1)})[1] for z in (acb(1), acb(0,1))]
    amplitude, offset = (a*a+b*b).sqrt(), abs((b/a).atan())
    assert a > 0 and code(ref,10) == 1011
    distance = ((arb(1012)/1024/amplitude).acos()-offset).union((arb(1011)/1024/amplitude).acos()+offset)
    assert arb(".145329") < distance < arb(".151932")
    boxes = [tuple(arb(v) for v in row) for row in BOXES]
    physical_cover(root, distance, boxes)
    angles = [(k, z, unit.arg()) for k,z,unit,_ in data]
    locks = previous["LOCKS"]
    u,v = (-next(angle for k,z,angle in angles if k == key and z == locks[k]) for key in DIRECTIONS[:2])
    rival = [(arb(0) if z == locks.get(k) else angle + k[0]*u + k[1]*v, 3) for k,z,angle in angles]
    aligned, rival_error = previous["geometric_tail"](rival, identically_zero=True)
    assert aligned == 6 and rival_error < 6085
    average = previous["kernel_sum"]
    ref_energy = [arb(0)] + [average(rival,n) for n in range(1,LIMIT)]
    target = arb(1)/8 + 12*arb("4.48e-9")
    for box, lift in product(boxes,(0,1)):
        def rows(x,y,d):
            native = {k:t for k,t in zip(DIRECTIONS,(x,y,d),strict=True)}
            result = []
            for k, _, angle in angles:
                shift = arb(0)
                if lift and k != (1,-1):
                    shift = -arb.pi() if angle > 0 else arb.pi()
                result.append((angle+native[k]+shift,3))
            return result
        xl,xr,yl,yr = box
        coarse = rows(xl.union(xr),yl.union(yr),distance)
        near,error = previous["geometric_tail"](coarse)
        assert near <= 3
        # K_N=(sin((N+1/2)α)/sin(α/2)-1)/(2N)；每个非共振项的常数为 -1/(2N)。
        lower = aligned-near-(rival_error+96-aligned+error-(96-near))/(2*LIMIT)
        assert lower > target, (near,error,lower)
        pending=[(xl,xr,yl,yr,range(1,LIMIT))]
        checked=0
        while pending:
            xl,xr,yl,yr,ns = pending.pop()
            x,y = xl.union(xr),yl.union(yr)
            if not (x-y).overlaps(distance):
                continue
            r = rows(x,y,(x-y).intersection(distance))
            nclose,err = previous["geometric_tail"](r)
            unresolved=[]
            for n in ns:
                upper = nclose+(err-(96-nclose))/(2*n)
                if ref_energy[n]-upper > target or ref_energy[n]-average(r,n) > target:
                    checked += 1
                else:
                    unresolved.append(n)
            if unresolved:
                if xr-xl>yr-yl:
                    mid=(xl+xr)/2
                    pending.extend(((xl,mid,yl,yr,unresolved),(mid,xr,yl,yr,unresolved)))
                else:
                    assert yr-yl>arb("1e-8"), (x,y,unresolved[:12])
                    mid=(yl+yr)/2
                    pending.extend(((xl,xr,yl,mid,unresolved),(xl,xr,mid,yr,unresolved)))
        print("相位框首坐标、提升、前缀比较格数与尾部余量",box[0].str(6),lift,checked,lower.lower().str(10),flush=True)


def nine_bits(root, exact, data, flavor):
    """平方帐篷给出无限单调见证。四项低阶校正产生能隙，允许明确截到有限阶。"""
    mul, conj, unit = (exact[k] for k in ("product","conjugate","unit"))
    width = arb(1)/50000
    epsilon = arb(10)**-17
    radius = arb(10)**-12
    last = 10**24
    delta = arb("4.48e-9")
    sets = [{z for key,z,_,_ in data if key == k} for k in DIRECTIONS]
    X,Y,D = sets
    for labels in sets:
        assert all(abs(unit(mul(z,conj(w))).arg()) > 2*width for z,w in combinations(labels,2))
    assert all(abs(unit(mul(mul(x,conj(y)),conj(d))).arg()) > 3*width for x,y,d in product(X,Y,D))
    chosen = ((7514,3),(21097,-103))
    mirror = tuple(conj(z) for z in chosen)
    for a,b in (chosen,mirror):
        assert any(k == (1,0) and z == conj(a) and soft == 0 for k,z,_,soft in data)
        assert any(k == (0,1) and z == conj(b) and soft == 0 for k,z,_,soft in data)
    peaks = {(conj(x),conj(y)) for x,y in product(X,Y)}
    peaks.update((conj(x),mul(conj(x),d)) for x,d in product(X,D))
    peaks.update((mul(conj(y),conj(d)),conj(y)) for y,d in product(Y,D))
    assert len(peaks) == 308
    poly = {1:6,19:-41,21:-69,22:16}
    def score(vertex):
        a,b = vertex
        clocks = dict(zip(DIRECTIONS,(a,b,mul(a,conj(b))),strict=True))
        return sum(3*sum(w*(unit(mul(z,clocks[k]))**q).real for q,w in poly.items())
                   for k,z,_,_ in data)
    target = score(chosen)
    gap = min(target-score(v) for v in peaks if v not in (chosen,mirror))
    assert gap > 70 and (target-score(mirror)).contains(0)
    norm = sum(abs(w) for w in poly.values())
    lip = sum(q*abs(w) for q,w in poly.items())
    # 两张支撑条带相交时，最近顶点的两个钟坐标各偏至多 2L，故各面至多偏 4L。
    other_gap = gap - 12*delta*norm - 4*width*96*lip
    assert other_gap > 49
    # 目标两条条带来自 x,y，各三张单位刚度面。目标框外至少损失 3r/L。
    gain_outside = 2*96*norm
    margin = 40*epsilon
    assert 3-epsilon*gain_outside > margin
    assert 3*radius/width-epsilon*gain_outside > margin
    assert epsilon*other_gap > margin
    # 精确傅里叶系数和低频校正。高于 22 阶的系数严格递减由书面导数证明。
    def coeff(q):
        return 4*(1-(q*width).sinc())/(arb.pi()*width*q*q)+epsilon*poly.get(q,0)
    assert all(coeff(q) > coeff(q+1) > 0 for q in range(1,max(poly)+1))
    tail = 192*8/(arb.pi()*width*last)
    assert tail < margin
    print("九位全环面交点数、单位多项式分差、邻域外能差及截尾差上界",
          len(peaks),gap.str(13),margin.str(8),tail.str(8),flush=True)
    reference = flavor["singlet"](11661)[1]
    physical = run_path(str(root/"ch39-first-three-mode-observables.py"))
    model, initial = physical["model_inputs"]()
    chiral = physical["joint"]["chiral"]
    for a,b in (chosen,mirror):
        x,y = unit(a).arg(),unit(b).arg()
        phases = {724:acb(0,x).exp(),804:acb(0,y).exp()}
        c = flavor["singlet_at"](phases)[1]+arb(0,2*radius)
        assert code(c,9) == code(reference,9) == 505
        assert code(c,10) == 1010 and code(reference,10) == 1011
        rows = chiral["source"]["spectral_supplier"](model)(arb("1e13"),clocks=phases)
        # 最大度至多 18。角漂移改变带荷拉普拉斯每行的绝对值之和至多 18|q|r。
        layers = {}
        for q in (0,*model["charges"]):
            error = 87*(4*delta+18*abs(q)*radius)+12*delta
            layers[q] = chiral["ordered"]([v+arb(0,error) for v,k in rows if k == q])
        assert all(layers[q][2] < 9 and layers[q][3] > 10 for q in model["charges"])
        _,alpha,curve=physical["layer_observables"](model,initial,layers)
        assert arb(".117") < alpha < arb(".119") and arb(".998") < curve < arb("1.002")
        print("九位候选框的余弦坐标、强耦合与曲率比",c.str(14),alpha.str(12),curve.str(12),flush=True)


def verify(root):
    exact = run_path(str(root / "ch39-mincut-real-triangle-tomography.py"))
    previous = run_path(str(root / "ch39-monotone-winding-action.py"))
    flavor = run_path(str(root / "ch39-chiral-neutral-flavor-phase.py"))
    ctx.prec=176
    data=geometry(root,exact)
    ten_bits(root,previous,data,flavor)
    nine_bits(root,exact,data,flavor)


if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
