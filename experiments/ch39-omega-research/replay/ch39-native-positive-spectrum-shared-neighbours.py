"""固定同一味标架和两块基准振幅，判定三枚实际编号的非负谱权可行性。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat, ctx


def positive_fit(columns, target, indices):
    floor = arb("1e-9")
    base = [floor*sum(col[k] for col in columns) for k in range(6)]
    matrix = arb_mat([columns[i-1] for i in indices]).transpose()
    assert not matrix.det().contains(0)
    gain = matrix.solve(arb_mat([[x-y] for x,y in zip(target,base,strict=True)]))
    assert all(x > 0 for x in gain.entries())
    weights = [floor for _ in columns]
    for i,x in zip(indices,gain.entries(),strict=True):
        weights[i-1] += x
    rebuilt = [sum(w*col[k] for w,col in zip(weights,columns,strict=True)) for k in range(6)]
    assert all((x-y).contains(0) for x,y in zip(rebuilt,target,strict=True))
    print("右邻下荷正增量谱序号、数值",indices,[x.str(12) for x in gain.entries()])


def separate(columns, target, offsets, sector, dual):
    # 将各编号的六维块按编号堆叠。共享谱权时，一列必须采用同一个谱序号。
    def evaluate(column):
        return sum(c*x for c,x in zip(dual,column,strict=True))

    col_values = []
    for j in range(19):
        stacked = [x for n in offsets for x in columns[n][sector][j]]
        col_values.append(evaluate(stacked))
    target_value = evaluate(target*len(offsets))
    margin = min(x.lower() for x in col_values)
    assert margin > 0 and target_value < 0
    # 六个实坐标按普通欧氏范数计，同一整数对偶向量给出目标到锥的距离下界。
    distance = -target_value/arb(sum(c*c for c in dual)).sqrt()
    print("偏移、部门、生成列最小值、目标值、距离下界",offsets,sector,
          margin.str(10),target_value.str(12),distance.lower().str(10))


def verify(root):
    mod = run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    ckm = run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    ctx.prec = 320
    center = 2_306_888_100_085
    rows = {k:mod["native"](root,center+k) for k in (-1,0,1)}
    columns = {k:[[mod["coord"](p) for p in sector] for sector in q] for k,q in rows.items()}
    v = acb_mat(ckm["standard_matrix"](*[arb(x).sin() for x in (".22704",".04275",".00377")],arb("1.139")))
    masses = [list(map(arb,row)) for row in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    up,down = [acb_mat([[m[i]/m[2] if i==j else 0 for j in range(3)] for i in range(3)]) for m in masses]
    adj = mod["adj"]
    up = adj(v)*up*v
    basis = mod["target_basis"](up,down)
    targets = []
    for q,h,picks in zip((0,1),(up,down),((1,6,13,16),(1,6,10,16)),strict=True):
        m = adj(basis)*h*basis
        block = mod["decompose"](rows[0][q],m,arb(".02145737"),arb("-.44525"),picks)
        targets.append(mod["coord"](block))

    positive_fit(columns[1][1],targets[1],(1,2,6,9,12,16))
    tests = [
        ((-1,),0,(316940,287236,-696693,575913,1,0)),
        ((-1,),1,(907284,-66930,-409217,70099,1,0)),
        ((1,),0,(532009,153927,-619776,-556026,1,0)),
        ((0,1),1,(119633,-780789,114831,57794,1,0,129597,567353,112498,-90500,1,0)),
    ]
    for offsets,sector,dual in tests:
        separate(columns,targets[sector],offsets,sector,dual)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
