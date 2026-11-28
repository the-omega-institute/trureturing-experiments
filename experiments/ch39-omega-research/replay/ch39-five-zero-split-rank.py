"""固定壳谱与 357 个格点时，五枚真实零点要求的四坐标仿射维数。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent


def systems(codes):
    parent = run_path(str(ROOT/'ch39-two-copy-two-split-lines.py'))
    load = parent['load']
    edges,_,heat,_ = load('two-copy-gram-cutoff-two-exceptions')['model_inputs']()
    _,lib,backgrounds = load('two-zero-cutoff-gravity')['gravity_fields'](
        ROOT,codes,tight_neutral=True)
    model,_,_,initial,_ = lib['chiral']['source']['mean']['inputs']()
    beta = [b/(2*arb.pi()) for b in model['beta']]
    logscale = (arb(10**16)/model['mz']).log()
    gauge = load('signed-gravity-unification')
    moments = load('signed-cutoff-gravity')['moments']
    fields = gauge['gauge_weights'](model)
    labels = ''.join(load('two-copy-four-split-lines')['LABELS'])
    for code,bg in zip(codes,backgrounds,strict=True):
        print('重谱零点编号',code,flush=True)
        system = parent['calibration_system'](edges,bg,heat,moments,gauge['threshold_column'],
                                              fields,labels,arb(7)/1000,initial,beta,logscale)
        yield code,system,parent['SEED'],beta,logscale


def verify():
    old = run_path(str(ROOT/'ch39-two-split-global-roots.py'))
    roots = [old['certify'](system,seed,beta,logscale)
             for _,system,seed,beta,logscale in systems(range(2_345_664,2_345_669))]
    differences = arb_mat([[v-c for v,c in zip(root,roots[0],strict=True)] for root in roots[1:]])
    d = differences.det()
    print("五个正交汇根坐标差的行列式",d.str(15),flush=True)
    assert arb("-1.31e-21") < d < arb("-1.29e-21")
    inverse = differences.inv()
    # 五个投影标量减去同一中心，最小的最大绝对偏差是极差的一半。
    # 极差不超过 1 的四个相对坐标是五维单位立方体减去首坐标的线性像。
    norms = []
    for vertex in product((0,1),repeat=4):
        normal = inverse*arb_mat([[v] for v in vertex])
        norms.append(sum(map(abs,normal.entries())))
    largest = max(v.lower() for v in norms).union(max(v.upper() for v in norms))
    tolerance = 1/(2*largest)
    normal = inverse*arb_mat([[0],[0],[0],[1]])
    fifth = 1/sum(map(abs,normal.entries()))
    assert arb("2.56e-9") < tolerance < arb("2.58e-9")
    assert arb("1.73e-8") < fifth < arb("1.75e-8")
    print("以任意三维仿射空间近似五个根的最优最大坐标误差",tolerance.str(9))
    print("第五个根至前四个根仿射包的最大坐标距离",fifth.str(9))


if __name__ == "__main__":
    verify()
