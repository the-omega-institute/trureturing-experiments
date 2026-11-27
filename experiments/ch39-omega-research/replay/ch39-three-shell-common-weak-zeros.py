"""同一个荷六源谱共轭同时实现第七、八、九壳原生第三枢轴的三个弱端点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


NUMBERS = (2_306_888_100_092, 2_306_888_100_093, 2_306_888_100_094)
LINES = (0, 15, 0)
WIDE_BOX = ((521,524), (894,897), (543,546))
BOX = (("522.398135899","522.398135902"), ("895.563345344","895.563345347"),
       ("544.705683197","544.705683201"))


def adj(matrix):
    return matrix.conjugate().transpose()


def unit(column):
    norm2 = (adj(column)*column)[0,0].real
    assert norm2 > 0
    return column / norm2.sqrt()


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-"+name+".py")))

    shell = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    shells = []
    for n in NUMBERS:
        nodes,basis = spectrum(shell,6,n)
        _,ordered,pivots,_ = next(s for s in native["shell_weights"](root,n) if s[0]==6)
        assert all(a.overlaps(b) for a,b in zip(nodes,ordered,strict=True))
        weights = [w.real for w in pivots[2]]
        assert sum(weights).contains(1)
        shells.append((adj(basis),weights))

    ctx.prec = 256
    # S 是三壳指定的本征行。S* (SS*)^-1 提供其共同右逆。
    rows = acb_mat([[basis[i,j] for j in range(19)] for (basis,_),i in zip(shells,LINES,strict=True)])
    gram = rows*adj(rows)
    print("三个指定行的 Gram 行列式",gram.det().real.str(13),flush=True)
    assert gram.det().real > arb(".005")
    lift = adj(rows)*gram.inv()
    def residual(v):
        return v - lift*(rows*v)
    v = unit(residual(acb_mat([[acb(2*j-13,(j*j+7)%13-6)] for j in range(19)])))
    w = residual(acb_mat([[acb((3*j*j+4)%19-9,j-7)] for j in range(19)]))
    w = unit(w-v*(adj(v)*w))
    source = fifth*columns
    source_gram = adj(source)*source
    assert all(source_gram[i,j].contains(int(i==j)) for i in range(2) for j in range(2))

    mp,mz = arb("1.22089e19"),arb("91.1876")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    floor,heavy,top = levels
    assert 0 < floor < heavy < top
    energy = [arb(a).union(arb(b)) for a,b in WIDE_BOX]
    b = acb_mat([[((q/mp-floor)/(heavy-floor)).sqrt()] for q in energy])
    particular = lift*b
    norm2 = (adj(particular)*particular)[0,0].real
    assert norm2 < arb("1e-12")
    u = particular + v*(1-norm2).sqrt()
    common = acb_mat([[u[i,0],w[i,0]] for i in range(19)])
    target_gram = adj(common)*common
    assert all(target_gram[i,j].contains(int(i==j)) for i in range(2) for j in range(2))
    pinned = rows*common
    assert all(pinned[i,0].overlaps(b[i,0]) and pinned[i,1].contains(0) for i in range(3))

    weak = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    for n,(basis,weights),i,wide,narrow in zip(NUMBERS,shells,LINES,WIDE_BOX,BOX,strict=True):
        assert weights[i] > 0
        projected = basis*common
        floor_q = min((floor+(heavy-floor)*abs(projected[j,0])**2+(top-floor)*abs(projected[j,1])**2
                       ).lower() for j in range(19) if j!=i)
        assert floor_q > arb("1e-10")
        # 18 条余线均在同一尾界上。指定行的截止用 S X=(b,0) 的精确恒等式。
        tail = (large(floor_q**2,p0,p1)-ref).abs_upper()*sum(abs(c).upper() for j,c in enumerate(weights) if j!=i)
        def faces(bounds):
            return [ref+weights[i]*(weak((arb(end)/mz)**2)-ref)+arb(0,tail) for end in bounds]
        outer = faces(wide)
        assert outer[0] > arb(".00024") and outer[1] < arb("-.00024")
        values = faces(narrow)
        print("壳、谱线、正权、其他振幅下界、其弱差总余量、左与右面",
              n,i+1,*[x.str(13) for x in (weights[i],floor_q,tail,*values)],flush=True)
        assert values[0] > 0 > values[1]
    return common


def verify(root):
    joint = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    return experiment(root,*joint["verify"](root,return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
