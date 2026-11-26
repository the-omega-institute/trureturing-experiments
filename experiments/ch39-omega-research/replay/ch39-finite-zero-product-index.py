"""四枚实际零点在五荷的整数面上链，以及两个四面体边界的直积配对。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb_mat, acb, arb, ctx


VERTICES = ((9, ".0005", "27"), (4, ".0005", "27"), (16, ".016", "7"), (10, ".0005", "27"))
FACES = ((1,2,3),(0,3,2),(0,1,3),(0,2,1))
CHARGES = (1,-4,2,-3,6)


def adj(matrix):
    return matrix.conjugate().transpose()


def frame(shell, charge, num, x, y):
    matrix = shell["laplacian"](charge, acb.zeta_zero(num).imag)
    for (i,j), value in zip(((0,2),(5,14)), map(arb, (x,y))):
        matrix[i,i] += value-1
        matrix[j,j] += value-1
        matrix[i,j] *= value
        matrix[j,i] = matrix[i,j].conjugate()
    eigen, vectors = matrix.eig(right=True,algorithm="rump")
    ix = sorted(range(len(eigen)),key=lambda j:float(eigen[j].real))
    assert all(v.imag.contains(0) for v in eigen)
    assert all(eigen[i].real < eigen[j].real for i,j in zip(ix[:3],ix[1:4]))
    columns = []
    for k in ix[:3]:
        col = [vectors[i,k] for i in range(19)]
        length = sum(abs(z)**2 for z in col).sqrt()
        columns.append([z/length for z in col])
    f = acb_mat(columns).transpose()
    minor = acb_mat([[f[i,j] for j in range(3)] for i in (0,1,2)]).det()
    assert abs(minor) > 0
    # 使前三行的行列式成为正实数，固定第三外幂的顶点规范。
    phase = minor.conjugate()/abs(minor)
    for i in range(19):
        f[i,0] *= phase
    return f


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 224
    chern = []
    data = []
    for q in CHARGES:
        frames = [frame(shell,q,*v) for v in VERTICES]
        angles, links, mods = {}, {}, []
        for i in range(4):
            for j in range(i+1,4):
                z = (adj(frames[j])*frames[i]).det()
                mods.append(abs(z).lower())
                angle = z.arg()
                assert abs(z)>0 and abs(angle)<arb.pi()
                angles[i,j], angles[j,i] = angle, -angle
                links[i,j], links[j,i] = z, z.conjugate()
        integers, phases = [], []
        for i,j,k in FACES:
            keys=((i,j),(j,k),(k,i))
            f=(links[keys[0]]*links[keys[1]]*links[keys[2]]).arg()
            assert abs(f)<arb.pi()
            n=((f-sum(angles[u] for u in keys))/(2*arb.pi())).unique_fmpz()
            assert n is not None
            integers.append(int(n))
            phases.append(f)
        d=(sum(phases)/(2*arb.pi())).unique_fmpz()
        assert d is not None and int(d)==sum(integers)
        chern.append(int(d))
        data.append(integers)
        print("荷、六边最小行列式模、四面相角、整数上链与总次数",
              q,min(mods).str(10),[p.str(10) for p in phases],integers,int(d))
    assert chern == [0,0,0,1,0]
    assert data == [[1,0,0,-1],[0,0,0,0],[-1,0,0,1],[0,0,0,1],[0,0,0,0]]
    products = [[sum(i*j for i in a for j in b) for b in data] for a in data]
    assert products == [[a*b for b in chern] for a in chern]
    print("两因子各自取五荷时的四维 Dirac 指标矩阵",products)
    multiplet = (6,3,3,2,1)
    for name,counts in (("五类型",(1,)*5),("完整一代",multiplet)):
        diagonal=sum(n*d*d for n,d in zip(counts,chern))
        twist=sum(n*d*q for n,d,q in zip(counts,chern,CHARGES))
        cx=sum(n*d for n,d in zip(counts,chern))
        cy=3*sum(n*q for n,q in zip(counts,CHARGES))
        print(name,"同荷外积指标、超荷线幂指标及后二者中线幂协议的第一陈类系数",
              diagonal,twist,(cx,cy))
    for a,b in (((1,-1,0),(1,-2,1)),((1,1,-1,-1,0),(1,1,-1,0,-1))):
        assert sum(a)==sum(b)==0 and sum(i*j for i,j in zip(a,b))==3
        print("行列式恒一的外积线和、指数",a,b,3)


if __name__=="__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
