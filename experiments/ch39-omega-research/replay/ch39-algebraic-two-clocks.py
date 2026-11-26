"""四个有理圆点给出的完整壳正权锥与半径九的四个第二陈数。"""
from itertools import permutations, product
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


def algebraic_ends(source):
    far = sorted(set(permutations((1, 5, 29))))
    near = sorted(set(permutations((13, 13, 23))))
    sites = far + near + sorted(
        p for p in product(range(30), repeat=3)
        if sum(v*v for v in p) == 867 and p not in far+near)
    distances = [[sum((a-b)**2 for a,b in zip(p,s)) for p in sites] for s in sites]
    geom = source["reconstruct"](distances, 0)
    assert len(sites) == 19 and len(geom[-1]) == 135
    # 半角 0、1/29、-4/17、1/6，端点符号分别为负、负、正、正。
    phases = ((acb(1), acb(-1), -acb(420,29)/421),
              (acb(1), acb(273,-136)/305, acb(35,12)/37))
    assert all((abs(v)-1).contains(0) for row in phases for v in row)
    def ends(q):
        start, end = [source["matrix"](geom, q, phi, 0) for phi in phases]
        changes = acb_mat(19, 2)
        for k, (i,j) in enumerate(((0,2),(5,14))):
            assert geom[-1][i,j][1] == 0
            changes[i,k] = 1
            changes[j,k] = start[j,i]
        return start, end-start, changes
    return ends


def verify(root):
    source = run_path(str(root/"ch39-unlabeled-native-connection.py"))
    chord = run_path(str(root/"ch39-two-index-chord.py"))
    ctx.prec = 256
    ends = algebraic_ends(source)
    adj, green, norm = (chord[k] for k in ("adj","green","frobenius"))
    start, change, B = ends(-4)
    eye = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    center_t, center_e = arb(".0198278251954"), arb("9.0078347438574")
    for _ in range(3):
        g, J, _ = green(start, change, B, eye, center_t, center_e)
        delta = J.solve(arb_mat([[g[0,1].real],[g[0,1].imag]]))
        center_t = (center_t-delta[0,0]).mid()
        center_e = (center_e-delta[1,0]).mid()
    gc, Jc, _ = green(start, change, B, eye, center_t, center_e)
    pre = Jc.inv()
    residual = max(abs(z).upper() for z in (pre*arb_mat([[gc[0,1].real],[gc[0,1].imag]])).entries())
    # 大盒控制小球唯一性，小盒给能量、边权与旁观谱的精确区间。
    identity = arb_mat([[1,0],[0,1]])
    for r in ("1e-9","1e-25"):
        radius = arb(r)
        t = center_t+arb(0,radius)
        energy = center_e+arb(0,radius)
        g, J, V = green(start, change, B, eye, t, energy)
        defect = identity-pre*J
        rho = max(sum(abs(defect[i,j]) for j in range(2)).upper() for i in range(2))
        assert rho<1 and residual+rho*radius<radius
    increments = [-1/g[i,i].real for i in range(2)]
    assert -1<increments[0]<0<increments[1] and 0<t<1
    below, gap = chord["gap_and_count"](start+center_t*change, center_e)
    assert below==3 and gap>(norm(change)+1)*radius
    v = acb_mat([[V[i,0]] for i in range(19)])
    v = v/(adj(v)*v)[0,0].real.sqrt()
    u = acb_mat([[V[i,1]] for i in range(19)])
    u = u-v*(adj(v)*u)[0,0]
    u = u/(adj(u)*u)[0,0].real.sqrt()
    frame = acb_mat([[v[i,0],u[i,0]] for i in range(19)])
    pauli = []
    tangent = [change]
    for j in range(2):
        col = acb_mat([[B[i,j]] for i in range(19)])
        tangent.append(col*adj(col))
    for M in tangent:
        X = adj(frame)*M*frame
        pauli.append([X[0,1].real,-X[0,1].imag,(X[0,0].real-X[1,1].real)/2])
    speed = arb_mat(pauli).det()
    assert speed>0
    print("根的弦坐标、能量、两条正电导", t.str(18), energy.str(18), *[(1+x).str(18) for x in increments])
    print("小盒压缩率上界、Pauli 速度有向行列式、基准谱隙", rho.str(10), speed.str(15), gap.str(10))
    lows, highs = [], []
    for q in (-4,1,2,-3,6):
        A,D,Q = ends(q)
        assert norm(D)<20
        if q==-4:
            continue
        L = A+t*D+Q*acb_mat([[increments[0],0],[0,increments[1]]])*adj(Q)
        center = acb_mat([[L[i,j].mid() for j in range(19)] for i in range(19)])
        center = (center+adj(center))/2
        error = norm(L-center)
        levels = center.eig(multiple=True)
        assert len(levels)==19 and all(z.imag.contains(0) for z in levels)
        ordered = sorted([z.real for z in levels], key=float)
        assert all(x<y for x,y in zip(ordered,ordered[1:]))
        third, fourth = [z+arb(0,error.upper()) for z in ordered[2:4]]
        assert third<9<fourth
        lows.append(third)
        highs.append(fourth)
        print("旁观荷、第三和第四特征值、半径九处第二陈数",q,third.str(16),fourth.str(16),3)
    low = max(lows, key=float)
    high = min(highs, key=float)
    assert all(low>v for v in lows if v is not low)
    assert all(high<v for v in highs if v is not high)
    assert 8<low<9<high<10
    print("四荷共用的开区间、唯一整数半径", low.str(16), high.str(16), 9)


if __name__=="__main__":
    here = Path(__file__).resolve().parent
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here
    del sys.argv[1:]
    verify(root)
