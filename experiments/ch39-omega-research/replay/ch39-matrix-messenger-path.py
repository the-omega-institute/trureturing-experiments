"""一般复矩阵中介的两条路径边可按奇异值精确消去。"""
from sympy import Matrix, Rational, cancel, expand, symbols


def verify():
    ell, mu, t, s, v = symbols("l mu t s v", real=True)
    x, y = symbols("x y", real=True)
    edge = (x*x-1)**2 + (y*y-1)**2
    assert expand(edge - ((v-1)**2 + 1 - 2*t**2).subs({v:x*x+y*y,t:x*y})) == 0
    upper = 2*ell*(t-1)**2 + mu*(t-s)**2
    lower = ell*(1-2*t**2) + mu*(t-s)**2
    half = Rational(1,2)
    assert expand(upper.subs(t,half)-lower.subs(t,half)) == 0
    assert expand(lower.diff(t) - (2*(mu-2*ell)*t-2*mu*s)) == 0
    star = (2*ell+mu*s)/(2*ell+mu)
    d = 2*ell*mu/(2*ell+mu)
    assert cancel(star-half - (2*ell-mu+2*mu*s)/(2*(2*ell+mu))) == 0
    assert cancel(upper-(2*ell+mu)*(t-star)**2-d*(s-1)**2) == 0
    # 对任意二维奇异值，精确分离均方根模长与偏离缩放酉方向的代价。
    loss = ((x-1)**2+(y-1)**2)/2
    polar_dist = ((x-s)**2+(y-s)**2)/2
    assert cancel((s*(loss-(s-1)**2)-polar_dist).expand()
                  -(s-1)*(x*x+y*y-2*s*s)/2) == 0
    # 两种原边使用同一中间点时，恒等式并不要求它们的 Gram 矩阵可交换。
    aa, bb, cc, dd, ee, ff = symbols("a b c d e f", real=True)
    m = Matrix([[aa,bb],[bb,cc]])
    n = Matrix([[dd,ee],[ee,ff]])
    identity = Matrix.eye(2)
    assert expand(((m-identity)**2+(n-identity)**2-(m+n-identity)**2
                   -identity+2*m*n).trace()) == 0
    print("Exact reduced penalty d times the sum of squared singular-value shifts for mu<=2l.")


if __name__ == "__main__":
    verify()
