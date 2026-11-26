"""三代折叠随 Higgs 背景的空间、有限有理检验及 Cayley 质量斜率。"""
import sympy as sp


def rational_cayley():
    t = sp.symbols("t", real=True)
    # 一个非对角的矩形物理块，同时检验左右补空间的次序。
    Y = sp.Matrix([[1,sp.Rational(1,3),0],[0,sp.Rational(1,2),sp.Rational(1,5)]])
    il, ir = sp.eye(Y.cols), sp.eye(Y.rows)
    a, b = t**2*Y.T*Y, t**2*Y*Y.T
    fl, fr = (il+a).inv(), (ir+b).inv()
    physical = 2*t*Y*fl
    U = physical.row_join((ir-b)*fr).col_join(((il-a)*fl).row_join(-physical.T))
    assert (U.T*U-sp.eye(5)).applyfunc(sp.cancel) == sp.zeros(5)
    assert U.subs(t,0)[:Y.rows,:Y.cols] == sp.zeros(2,3)
    assert (physical.diff(t).subs(t,0)-2*Y) == sp.zeros(2,3)
    y = 2*t/(1+t**2)
    sensitivity = sp.cancel(t*y.diff(t)/y)
    assert sp.cancel(sensitivity**2+y**2-1)==0 and sensitivity.subs(t,1)==0


def degree_two_witness():
    # 四个不同有理背景上的严格折叠见证。Gram 的正定性由有理主子式给出。
    x=sp.symbols("x", real=True)
    a,b=sp.Rational(1,96),sp.Rational(1,48)
    t=x+sp.Rational(1,2)
    product4=(x*x-a*a)*(x*x-b*b)
    polynomial = 1-t**2+256*product4
    c=sp.Poly(polynomial,x)
    gram=sp.Matrix([[c.nth(0),c.nth(1)/2,(c.nth(2)-1)/2],
                    [c.nth(1)/2,1,0],[(c.nth(2)-1)/2,0,256]])
    powers=sp.Matrix([1,x,x*x])
    assert sp.expand((powers.T*gram*powers)[0]-polynomial)==0
    assert all(gram[:k,:k].det()>0 for k in (1,2,3))
    for sample in (-a,-b,b,a):
        assert (t*t+polynomial).subs(x,sample)==1
    assert (t*t+polynomial).subs(x,0)!=1


def sm_sizes():
    # 场的维数含颜色重数。H=0 时最少辅助表示是对侧物理表示的精确镜像。
    left=3*(3*2+2)
    right=3*(3+3+1)
    assert (left,right,left+right)==(24,21,45)
    assert left+right-3==42
    assert left+right+3==48


if __name__=="__main__":
    sm_sizes()
    rational_cayley()
    degree_two_witness()
    print("H=0 每侧最少维数 45；加入三个右手中微子时为 48")
    print("已核验矩形有理 Cayley 恒等式、质量斜率及四个有理背景上的二次多项式")
