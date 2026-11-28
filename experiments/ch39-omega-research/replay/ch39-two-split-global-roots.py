"""在冻结的 357 个格点下，认证三个背景在四维参数立方体上的唯一根。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
PARTS = 256
TAIL = 244


def point_matrix(matrix):
    return arb_mat([[matrix[i, j].mid() for j in range(matrix.ncols())]
                    for i in range(matrix.nrows())])


def join_columns(columns):
    return arb_mat([column.entries() for column in columns]).transpose()


def hull(columns):
    values = columns[0].entries()
    for column in columns[1:]:
        values = [v.union(w) for v, w in zip(values, column.entries(), strict=True)]
    return arb_mat([[v] for v in values])


def certify(system, seed, beta, logscale):
    readings, response, derivative, directions, project, target, jacobian = system
    center = [arb(z) for z in seed]
    for _ in range(2):
        correction = jacobian(center).solve(project*readings(center)-target)
        center = [(v-correction[i, 0]).mid() for i, v in enumerate(center)]

    # 先证局部存在，再在初值无关的整个闭立方体上排除第二个解。
    pre = point_matrix(point_matrix(jacobian(center)).inv())
    assert not pre.det().contains(0)
    residual = pre*(project*readings(center)-target)
    radius = arb(2)**-36
    box = [v+arb(0, radius) for v in center]
    assert all(0 < v < 1 for v in box)
    identity = arb_mat([[int(i == j) for j in range(4)] for i in range(4)])
    error = identity-pre*jacobian(box)
    contraction = max(sum(abs(error[i,j]) for j in range(4)).upper() for i in range(4))
    displacement = max(abs(v).upper() for v in residual.entries())
    assert contraction < 1 and displacement + radius*contraction < radius
    root = [v+arb(0, (displacement/(1-contraction)).upper()) for v in center]
    common = readings(root)[4, 0]
    minimum = common + beta[2]*logscale
    assert beta[2] < beta[1] < 0 < beta[0] and logscale > 0 and minimum > 1

    projection = pre*project
    d0, d16 = (projection*v for v in directions)
    unit = arb(0).union(1)
    cells = [(arb(i)/PARTS).union(arb(i+1)/PARTS) for i in range(PARTS)]
    d338 = hull([projection*derivative(338,t) for t in cells])
    # 积分型中值公式。两个份额的导数恒定，最后一列用全区间逐分量包络。
    fixed = residual + d0*(unit-center[0]) + d16*(unit-center[1]) + d338*(unit-center[3])
    reference = response(17,center[2])
    excluded = []
    remainder = []
    for i, t in enumerate(cells):
        slope = projection*derivative(17,t)
        if i >= TAIL:
            remainder.append(slope)
        else:
            mid = t.mid()
            values = fixed + projection*(response(17,mid)-reference) + slope*(t-mid)
            excluded.append(values[1,0] if i < 4 else -values[2,0])
    lower_first = min(v.lower() for v in excluded[:4])
    lower_next = min(v.lower() for v in excluded[4:])
    det = join_columns((d0,d16,hull(remainder),d338)).det()
    assert lower_first > arb("1.092") and lower_next > arb("0.01539")
    assert det > arb("0.93885"), det
    print('排除区间的最小有向残差与余区平均雅可比行列式下界',
          lower_first.str(9), lower_next.str(9), det.lower().str(9))
    print('根坐标、共同逆耦合、全程最低', [v.str(14) for v in root],
          common.str(12), minimum.str(12))
    return root


def verify():
    parent = run_path(str(ROOT/'ch39-two-copy-two-split-lines.py'))
    load = parent['load']
    edges, central, heat, _ = load('two-copy-gram-cutoff-two-exceptions')['model_inputs']()
    _, lib, adjacent = load('two-zero-cutoff-gravity')['gravity_fields'](
        ROOT,(2_345_665,2_345_667),tight_neutral=True)
    model, _, _, initial, _ = lib['chiral']['source']['mean']['inputs']()
    gauge = load('signed-gravity-unification')
    moments = load('signed-cutoff-gravity')['moments']
    fields = gauge['gauge_weights'](model)
    logscale = (arb(10**16)/model['mz']).log()
    beta = [v/(2*arb.pi()) for v in model['beta']]
    labels = ''.join(load('two-copy-four-split-lines')['LABELS'])
    roots = []
    for n, rows in zip((2_345_666,2_345_665,2_345_667),(central,*adjacent),strict=True):
        print('重谱零点编号', n, flush=True)
        system = parent['calibration_system'](edges,rows,heat,moments,gauge['threshold_column'],
                                              fields,labels,arb(7)/1000,initial,beta,logscale)
        roots.append(certify(system,parent['SEED'],beta,logscale))
    for root in roots[1:]:
        drift = [v-w for v,w in zip(root,roots[0],strict=True)]
        assert all(abs(v) > arb('4.3e-7') for v in drift)
        print('邻点四参数减中央', [v.str(12) for v in drift])


if __name__ == '__main__':
    verify()
