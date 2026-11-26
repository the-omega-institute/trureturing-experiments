"""两个真实零点对逐边实磁动能系数和逐站势的唯一性与稳定界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


local = run_path(str(Path(__file__).with_name("ch39-local-light-generator.py")))
mix = local["mix"]
shell = mix["shell"]


def leakage_map(number):
    modes = mix["light_modes"](number)
    complement = local["identity"](19) - modes * local["adjoint"](modes)
    lap = shell["laplacian"](1, acb.zeta_zero(number).imag)
    # 先保留所有列，直接核对两个公共零方向，再各删去一个代表。
    site_columns = [[complement[row,i]*modes[i,k] for row in range(19) for k in range(3)]
                    for i in range(19)]
    edge_columns = []
    currents = []
    for i,j,_,_ in shell["geo"]["edges"]:
        weight = abs(lap[i,j])
        phase = -lap[i,j] / weight
        column = [(complement[row,i] - phase.conjugate()*complement[row,j])
                  * (modes[i,k] - phase*modes[j,k]) * weight
                  for row in range(19) for k in range(3)]
        edge_columns.append(column)
        currents.append([acb(0,1)*weight*(phase*complement[row,i]*modes[j,k]
                         - phase.conjugate()*complement[row,j]*modes[i,k])
                         for row in range(19) for k in range(3)])
    for group in (site_columns,edge_columns):
        assert max(abs(sum(entries)) for entries in zip(*group)) < arb("1e-38")
    both = site_columns[1:] + edge_columns[1:]
    rows = list(zip(*(both + currents)))
    return arb_mat([[z.real for z in row] for row in rows] + [[z.imag for z in row] for row in rows])


def sigma_floor(gram):
    inverse = gram.inv()
    return 1/sum(x*x for x in inverse.entries()).sqrt().sqrt()


def verify():
    ctx.prec = 448
    nodes = (28325,48067,58214,22027,11661)
    grams = []
    for number in nodes:
        leakage = leakage_map(number)
        grams.append(leakage.transpose()*leakage)
    two = grams[0] + grams[1]
    real_weights = arb_mat([[two[i,j] for j in range(152)] for i in range(152)])
    floor = sigma_floor(real_weights)
    print("逐边实动能，两枚零点的联合泄漏奇异值下界",floor.str(12))
    assert floor > arb("2.885e-7")
    five = two + grams[2] + grams[3] + grams[4]
    assert five.nrows() == five.ncols() == 287
    floor = sigma_floor(five)
    print("加独立有向边电流，五枚零点的联合泄漏奇异值下界", floor.str(12))
    assert floor > arb("1.879e-11")
    fixed = [(i,j,d) for i,j,d,_ in shell["geo"]["edges"] if i >= 6 and j >= 6]
    assert len(fixed) == 13*12//2
    tuned = set(shell["soft"]["star"]["orbit"])
    pole = shell["soft"]["pole"]
    assert all(d not in (724,804) and not (pole in (i,j) and tuned.intersection((i,j)))
               for i,j,d in fixed)


if __name__ == "__main__":
    verify()
