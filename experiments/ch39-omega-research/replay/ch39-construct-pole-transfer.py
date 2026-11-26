"""由荷六三项公共根显式构造正交谱基及极点搬移概率上界。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


# 每个三元组指定目标行、余空间内的高轴、低轴，索引从零起。
PLANS = {
    -1: ([(0,18,15),(15,6,1),(18,18,13),(17,17,14),(2,6,9),(11,7,5),
          (10,7,12),(16,6,16),(14,6,11),(13,7,4),(6,17,2),(7,18,3),
          (9,7,17),(5,7,0),(1,6,8),(8,7,18),(4,6,10),(3,7,6)], {0,5,12}),
    1: ([(0,17,11),(15,17,4),(18,18,5),(17,17,16),(14,17,9),(3,17,7),
         (10,18,3),(5,18,6),(7,17,14),(4,17,13),(12,17,2),(9,17,15),
         (6,17,0),(2,18,10),(1,18,17),(16,18,1),(8,18,8),(13,18,12)], {0,3,6,9,12,15,17}),
}


def mixing_matrix(values, goals, plan):
    size = len(values)
    eigenvalues = dict(enumerate(values))
    basis = {i: [arb(int(i == j)) for j in range(size)] for i in range(size)}
    rows = [None]*size
    for j, hi, lo in plan:
        a, b, y = eigenvalues[hi], eigenvalues.pop(lo), goals[j]
        assert b < y < a and rows[j] is None
        c, s = ((y-b)/(a-b)).sqrt(), ((a-y)/(a-b)).sqrt()
        v, w = basis[hi], basis.pop(lo)
        rows[j] = [c*x+s*z for x, z in zip(v, w, strict=True)]
        basis[hi] = [s*x-c*z for x, z in zip(v, w, strict=True)]
        eigenvalues[hi] = a-(y-b)
    remaining, = eigenvalues
    j, = [i for i, row in enumerate(rows) if row is None]
    assert eigenvalues[remaining].overlaps(goals[j])
    rows[j] = basis[remaining]
    orthogonal = arb_mat(rows)
    gram = orthogonal * orthogonal.transpose()
    assert all(gram[i, j].contains(int(i == j)) for i in range(size) for j in range(size))
    assert all((sum(a*x*x for a, x in zip(values, row, strict=True))-y).contains(0)
               for row, y in zip(rows, goals, strict=True))
    return orthogonal


def norm_lower(vector):
    return sum(max(arb(0), abs(z).lower())**2 for z in vector).sqrt().lower()


def verify(root):
    center = 2_306_888_100_085
    sh = run_path(str(root/'ch39-schur-horn-gravity-holdout.py'))
    shell = run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    transport = run_path(str(root/'ch39-adjacent-pvm-transport.py'))
    pivots = run_path(str(root/'ch39-native-rh-ray-pvm.py'))
    gravity = run_path(str(root/'ch39-three-index-free-atoms.py'))
    high_weak = run_path(str(root/'ch39-full-neutral-free-cutoff.py'))['weak_column']
    # 两个既有训练构造以同一数值对象填充每个等幅组。
    left = run_path(str(root/'ch39-left-cutoff-two-turns.py'))['verify'](root)
    right = run_path(str(root/'ch39-adjacent-ordered-charge-six-cutoffs.py'))['verify'](root)[2]
    moments, heat, fields = gravity['broken_fields'](root, [center-2, center+2])
    ctx.prec = 320
    weak = pivots['weak_series']()
    for side, (nodes, positions), field in zip((-1, 1), (left, right), fields, strict=True):
        target = center+2*side
        levels, src = transport['spectrum'](shell, 6, center+side)
        _, dest = transport['spectrum'](shell, 6, target)
        assert all(v.overlaps(e) for v, e in zip(nodes, levels, strict=True))
        values = [x.sqrt() for x in positions]
        _, _, cross, _ = next(t for t in pivots['shell_weights'](root, target) if t[0] == 6)
        ctx.prec = 320
        goals = sh['certify'](target, values, [z.real for z in cross[0]], field, moments, heat, weak, high_weak)
        plan, negative = PLANS[side]
        rotation = mixing_matrix(values, goals, plan)
        pole = shell['soft']['pole']
        assert pole == 14
        p = [abs(src[pole, j]) for j in range(19)]
        q = arb_mat([[(-1 if j in negative else 1)*abs(dest[pole, j])] for j in range(19)])
        back = (rotation.transpose()*q).entries()
        medium = {6, 7, 17} if side < 0 else {17}
        groups = [set(range(18))-medium, medium, {18}]
        assert all(len({id(positions[j]) for j in group}) == 1 for group in groups)
        # 组内正交变换保持整个输出对角。把每个源投影转到回拉投影方向即可同相相加。
        keep = sum(norm_lower([p[i] for i in group])*norm_lower([back[i] for i in group]) for group in groups)
        moved = (1-keep**2).upper()
        assert moved < (arb('.12341912') if side < 0 else arb('.04732188'))
        print('源壳、目标壳、极点保留模下界、站外转移概率上界', center+side, target,
              keep.lower().str(13), moved.str(13), flush=True)


if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)
