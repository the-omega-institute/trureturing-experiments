"""用高斯有理矩阵与符号二阶块核验共同三代压缩的 CP 交换子分解。"""
from sympy import I, Matrix, diag, eye, expand, symbols, zeros


def commutator(a, b):
    return a*b-b*a


def verify():
    # 二阶代数恒等式对任意复矩阵成立。此处自由符号互不相关。
    m = Matrix(2, 2, symbols("m:4"))
    b = Matrix(2, 2, symbols("b:4"))
    c = commutator(m, b)
    assert (m*c*m-m.det()*c).applyfunc(expand) == zeros(2)
    a, b1, b2, r1, r2, l1, l2 = symbols("a b1 b2 r1 r2 l1 l2")
    c2 = Matrix([[a,b1],[b2,-a]])
    right = Matrix([r1,r2])
    left = Matrix([[l1,l2]])
    c3 = c2.row_join(right).col_join((-left).row_join(zeros(1)))
    assert expand(c3.det()+(left*c2*right)[0]) == 0

    # 以非单位法向量形成精确有理正交投影。三阶无迹交换子的行列式等于三次迹除三。
    tu = Matrix([[1, 2+I], [2-I, 4]])
    td = Matrix([[3, 1-2*I], [1+2*I, -2]])
    real_cases = ((2,7,-1,5), (2,7,5,5), (7,7,-1,5), (7,7,5,5))
    normals = (Matrix([1,I,2-I,3]), Matrix([I,3-I,0,2]), Matrix([0,0,2,3]))
    for ua,ub,da,db in real_cases:
        u,d = diag(tu,ua,ub), diag(td,da,db)
        delta = (da-db)*(tu-ua*eye(2))+(ub-ua)*(td-da*eye(2))
        for normal in normals:
            norm = (normal.H*normal)[0].expand()
            p = (eye(4)-normal*normal.H/norm).applyfunc(expand)
            h,k = [(p*a*p).applyfunc(expand) for a in (u,d)]
            cubic = commutator(h,k).applyfunc(expand)
            det = ((cubic**3).trace()/3).expand()
            v = normal[:2, :]
            q = normal[2]*normal[2].conjugate()*normal[3]*normal[3].conjugate()
            predicted = -q*delta.det()*(v.H*commutator(tu,td)*v)[0]/norm**3
            assert expand(det-predicted) == 0
    print("两个符号块恒等式和十二个有理压缩均精确成立。")


if __name__ == "__main__":
    verify()
