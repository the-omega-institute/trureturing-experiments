"""原十九站壳上，有限个实际零点排除若干按轨道规模等权的正混合。

按六个三站轨道的同时置换分组。每组固定上下子集的三站轨道数、是否含
单站极点及公共三站轨道数。组内对有序子集对等权，组间可取任意概率。
作用是三个字符最低根上的六个参考因子的平方积之和，使用原能量单位。
"""
from collections import defaultdict
from functools import cache
from itertools import product
from math import comb, prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


TARGET = {-4: 2345666, 2: 964890}
MODES = ("D", "I")
MASKS = tuple(mask for mask in range(1, 128) if mask != 64)
MODE_PAIRS = tuple(product(MODES, repeat=2))


def geometry_classes():
    classes = defaultdict(list)
    for up, down in product(MASKS, repeat=2):
        u, d = up & 63, down & 63
        key = (u.bit_count(), bool(up & 64), d.bit_count(), bool(down & 64),
               (u & d).bit_count())
        classes[key].append((up, down))
    assert len(MASKS) == 126 and len(classes) == 284
    for (u, _, d, _, overlap), pairs in classes.items():
        assert len(pairs) == comb(6,u)*comb(u,overlap)*comb(6-u,d-overlap)
    return classes


def verify(replay):
    def load(name):
        return run_path(str(replay / f"ch39-{name}.py"))

    shell = load("three-ray-chiral-boundary")
    pump = load("native-light-projector-pump")
    reference = load("geometric-selector-mixture")
    ctx.prec = 320
    triple = reference["geometry"](shell, pump)
    frames = list(pump["rotation_bases"](shell))
    adjoint = pump["adjoint"]
    roots = {charge: [triple(charge,mask,mode) for mask,mode in product(MASKS,MODES)]
             for charge in TARGET}
    classes = geometry_classes()
    row = {(mask,mode): i for i,(mask,mode) in enumerate(product(MASKS,MODES))}
    assert all(len(items) == 252 for items in roots.values())

    @cache
    def action(charge, index):
        matrix = shell["laplacian"](charge, acb.zeta_zero(index).imag)
        energies = [reference["minimum"](adjoint(frame)*matrix*frame) for frame in frames]
        factors = []
        for source_charge in (-4,2):
            factors.append(arb_mat([[prod((energy-root)**2 for root in ref) for energy in energies]
                                    for ref in roots[source_charge]]))
        # 矩阵的 (i,j) 项是指定上下参考对的三字符作用。
        return factors[0] * factors[1].transpose()

    @cache
    def differences(charge, rival):
        assert rival != TARGET[charge]
        return action(charge,TARGET[charge]) - action(charge,rival)

    def averaged(competitors):
        # 以下权重只组合用于反证的作用差；每项都相对于本荷的中央编号。
        denominator = sum(weight for _,_,weight in competitors)
        assert all(weight > 0 for _,_,weight in competitors)
        difference = sum((weight*differences(charge,index) for charge,index,weight in competitors),
                         arb_mat(252,252)) / denominator
        return {key: {mode: sum(difference[row[u,mode[0]],row[d,mode[1]]]
                                for u,d in pairs)/len(pairs) for mode in MODE_PAIRS}
                for key,pairs in classes.items()}

    def certify(label, values, bound):
        lower = min(value.lower() for value in values)
        assert lower > bound
        print(label, "最小作用差下界", lower.str(14))

    # 每个规模与交叠类均可独立选下列边界概率律，再混合各几何类。
    # 四种模式等权、DD/II 各半、固定 II 三者共享同一组竞争编号。
    balanced = averaged(((2,1,1), (-4,11874,1080), (-4,21276,919)))
    certify("四边界模式各四分之一", (sum(v.values())/4 for v in balanced.values()), 8)
    certify("共享边界且 D/I 各半",
            ((v["D","D"]+v["I","I"])/2 for v in balanced.values()), 8)
    certify("固定 II 与上述规则的共同证书", (v["I","I"] for v in balanced.values()), 8)

    induced = averaged(((-4,11874,1),))
    certify("固定诱导图", (v["I","I"] for v in induced.values()), 900)

    opposite = averaged(((2,1132,1121), (-4,11035,879)))
    certify("DI/ID 各半", ((v["D","I"]+v["I","D"])/2 for v in opposite.values()), 2900)

    # 上下规模与极点标志相同。在每个交叠类内可任意混合 DD、II。
    shared = averaged(((2,1,439), (-4,3200,24561)))
    certify("共用规模且共用边界",
            (v[mode,mode] for (u,pu,d,pd,_),v in shared.items()
             if (u,pu)==(d,pd) for mode in MODES), 3300)

    # 上下子集独立、分别均匀取遍 126 个；边界联合律与掩码独立，允许任意相关。
    uniform = averaged(((-4,24685,1),))
    certify("全部子集分别均匀、任意固定边界联合律",
            (sum(len(classes[key])*v[mode] for key,v in uniform.items())/len(MASKS)**2
             for mode in MODE_PAIRS), 16000000)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
