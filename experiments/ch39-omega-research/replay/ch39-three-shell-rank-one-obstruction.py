"""从第五壳完整见证出发，排除单个秩一酉扰动在三个壳同时消去弱差。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


CANDIDATES = ((0,9,15),(0,9,15),(0,15))
THRESHOLDS = (1294,1041,2241)
HIGH_CORRECTION = (".0000824096",".0001280487",".0000272865")


def adj(a):
    return a.conjugate().transpose()


def minimum_gain(rows, fifth, columns):
    # 预条件器只取所给中点；其酉性由 Hermitian 特征向量的正交性给出。
    y = (rows*(fifth*columns)).mid()
    values,vectors = (adj(y)*y).eig(right=True,algorithm="rump")
    assert all(z.real > 0 and z.imag.contains(0) for z in values)
    assert values[0].real < values[1].real or values[1].real < values[0].real
    right = []
    left = []
    for i in range(2):
        col = acb_mat([[vectors[j,i]] for j in range(2)])
        col = col/(adj(col)*col)[0,0].real.sqrt()
        right.append(col)
        image = y*col
        left.append(adj(image)/(adj(image)*image)[0,0].real.sqrt())
    f = acb_mat([[left[i][0,j] for j in range(3)] for i in range(2)])
    v = acb_mat([[right[j][i,0] for j in range(2)] for i in range(2)])
    # 先线性合并高精度的三个谱行，再传递共同见证的区间误差。
    a = ((f*rows)*fifth)*(columns*v)
    frobenius = sum(abs(x)**2 for x in a.entries()).sqrt()
    return (abs(a.det())/frobenius).lower()


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    ctx.prec = 256
    mp,mz = arb("1.22089e19"),arb("91.1876")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    ref = arb("-.074256577194663365545")
    native = lib("native-rh-ray-pvm")
    finite = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    join = arb(10)**10
    lower = arb("-1.4e-11")
    full = arb(".020483945")

    def cover(start, upper):
        assert start > 1 and lower < large(join*p0,p0,p1)-ref < upper
        stack = [(start.lower(),join)]
        count = 0
        while stack:
            l,r = stack.pop()
            value = finite(l.union(r))-ref if r < 2*l else None
            if value is not None and lower < value < upper:
                count += 1
                continue
            assert r-l > arb("1e-14")
            mid = ((l+r)/2).mid()
            stack.extend(((l,mid),(mid,r)))
        return count

    floor,heavy,top = levels
    assert 0 < floor < heavy < top
    count = cover(floor**2/p0, full)
    shell = lib("three-ray-chiral-boundary")
    bases = []
    for j,(chosen,q,corr) in enumerate(zip(CANDIDATES,THRESHOLDS,HIGH_CORRECTION,strict=True)):
        n = 2_306_888_100_092+j
        nodes,basis = lib("adjacent-pvm-transport")["spectrum"](shell,6,n)
        _,ordered,pivots,_ = next(v for v in native["shell_weights"](root,n) if v[0]==6)
        assert all(a.overlaps(b) for a,b in zip(nodes,ordered,strict=True))
        weights = [x.real for x in pivots[2]]
        assert sum(weights).contains(1) and all(not w.contains(0) for w in weights)
        assert all(weights[k] > 0 for k in chosen)
        local = cover((arb(q)/mz)**2,arb(corr))
        upper = ref+sum(w*(arb(corr) if k in chosen else (full if w>0 else lower))
                        for k,w in enumerate(weights))
        print("壳、能标、闭段数、候选全部达到能标时的弱差上界",n,q,count,local,upper.str(15),flush=True)
        assert upper < 0
        bases.append(adj(basis))

    cap = sum((arb(q)/mp-floor)/(heavy-floor) for q in THRESHOLDS).sqrt().upper()
    gains = {}
    for selection in product(*CANDIDATES):
        rows = acb_mat([[b[i,k] for k in range(19)] for b,i in zip(bases,selection,strict=True)])
        gain = minimum_gain(rows,fifth,columns)
        gains[selection] = gain
        print("三线编号、最小奇异值下界",tuple(k+1 for k in selection),gain.str(14),flush=True)
    minimum = min(gains.values())
    print("十八组三线一致下界、弱差给出的上界",minimum.str(15),cap.str(15),flush=True)
    assert len(gains)==18 and minimum > arb("3.64675622e-8") and cap < arb("2.85236572e-8")
    robust = (minimum-cap)/arb(3).sqrt()
    assert robust > arb("4.58e-9")
    print("每个解的第二奇异值 σ₂(U−I) 严格下界",robust.lower().str(14),flush=True)
    return minimum,cap


def verify(root):
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*source["verify"](root,return_operator="spectral"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
