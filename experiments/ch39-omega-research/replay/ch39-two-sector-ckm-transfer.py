"""两个一跳轻空间等距共同承担参考 CKM 的轮换偏差。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    light = load("local-yukawa-plucker")
    cyclic = load("joint-zero-cyclic-light-transfer")
    prev = load("joint-zero-diagonal-light-filters")
    ckm = load("shell-clock-ckm-blocks")
    ctx.prec = 320
    adj = cyclic["adjoint"]
    def identity(n):
        return acb_mat([[int(i==j) for j in range(n)] for i in range(n)])
    def frobenius(a):
        return sum(abs(z)**2 for z in a.entries()).sqrt()
    # V 的主平方根用单根处的谱投影多项式求出，既容许任意列相位，
    # 此例固定原参考的标准参数化。下界覆盖任意外部行列重定相。
    v = acb_mat([list(row) for row in ckm["standard_matrix"](*ckm["centers"])])
    eig = v.eig()
    assert len(eig)==3 and all(z.real>0 for z in eig)
    half = acb_mat(3,3)
    for i,z in enumerate(eig):
        projector = identity(3)
        for j,other in enumerate(eig):
            if i!=j:
                assert not (z-other).contains(0)
                projector = (projector*(v-other*identity(3)))/(z-other)
        half += z.sqrt()*projector
    assert all(z.contains(0) for z in (half*half-v).entries())
    shell = light["shell"]
    height = acb.zeta_zero(2_306_888_100_085).imag
    delta = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    center = 15-delta/3
    rot,edges = cyclic["generator"](shell)
    blocked = [(i,j) for i in range(19) for j in range(19)
               if i!=j and tuple(sorted((i,j))) not in edges]
    u = light["light_basis"](1,height,center)
    chars = cyclic["character_order"](u,rot)
    assert chars == [0,2,1]
    omega = acb(-arb(1)/2,arb(3).sqrt()/2)
    rotation = acb_mat([[([1,omega,omega.conjugate()][chars[i]] if i==j else 0)
                        for j in range(3)] for i in range(3)])
    # 允许三条质量线重新标号时，旧六置换证书已证参考中心最佳匹配为对角。
    load("ckm-cyclic-covariance-cost")["verify"]()
    defect = (3*sum(abs(v[i,j])**2 for i in range(3) for j in range(3) if i!=j)).sqrt()
    lower = (defect**2/4+defect**4/192).sqrt()
    upper = frobenius(rotation*half-half*rotation)
    print("参考中心最佳双扇区最大偏差的下界与主平方根实现",lower.str(18),upper.str(18))
    assert arb(".28117247")<lower<upper<arb(".28208788")
    comp_source = identity(19)-u*adj(u)
    for q,m in ((-4,adj(half)),(2,half)):
        target = light["light_basis"](q,height,center)
        assert cyclic["character_order"](target,rot)==chars
        comp_target = identity(19)-target*adj(target)
        h = acb_mat([[comp_target[i,k]*comp_source[l,j] for k,l in blocked] for i,j in blocked])
        prev["certify_positive"](h)
        base = target*m*adj(u)
        rhs = acb_mat([[base[i,j]] for i,j in blocked])
        norm2 = (adj(rhs)*h.solve(rhs))[0,0].real
        assert norm2>0 and norm2<1
        print("目标荷、一跳补块的 Frobenius 范数",q,norm2.sqrt().str(16))


if __name__=="__main__":
    verify()
