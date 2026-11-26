"""带荷矢量对的有限谱积分给十九维中性自能的欧氏相对算子界。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx

census = run_path(str(Path(__file__).with_name("ch39-shell-charged-channel-census.py")))


def subtracted_integral(a,b,steps):
    """在 Q²=1 对两次减法的核作上下 Darboux 包围，提出规范和站点因子。"""
    first,second = (a+b)**2,(a-b)*(a-b)
    if first > 1:
        return None
    assert first < 1
    x,z = a**2,b**2
    def emission(y):
        beta2 = (1-first/y)*(1-second/y)
        return beta2.sqrt()**3*(y**2+10*y*(x+z)+x**2+10*x*z+z**2)
    step = (1-first)/steps
    left,previous = first,arb(0)
    lower,upper = arb(0),arb(0)
    for j in range(1,steps+1):
        right = first+j*step
        current = emission(right)
        lower += previous/(right*(1+right))
        upper += current/(left*(1+left))
        left,previous = right,current
    factor = step/(192*arb.pi()**2*x*z)
    return (lower*factor).lower(),(upper*factor).upper()


def positive_definite(matrix):
    # Hermiticity follows from the construction. Sylvester's criterion is
    # invariant under the nonorthogonal coordinates used for the heavy block.
    for k in range(1,matrix.nrows()+1):
        minor = acb_mat([[matrix[i,j] for j in range(k)] for i in range(k)])
        determinant = minor.det().real
        if determinant < 0:
            return False
        assert determinant > 0
    return True


def verify():
    ctx.prec = 256
    c,n,r = census["mass_matrices"]()
    modes = census["modes"](c,r)
    upper,lower = acb_mat(19,19),acb_mat(19,19)
    for i,left in enumerate(modes):
        for j,right in enumerate(modes[:i+1]):
            interval = subtracted_integral(left.mass,right.mass,512)
            if interval is None:
                continue
            h = acb_mat([[a.conjugate()*b] for a,b in zip(left.profile,right.profile)])
            gram = 19*h*h.conjugate().transpose()
            if i != j:
                gram = gram+gram.conjugate()
            lower += gram*interval[0]
            upper += gram*interval[1]
    metric = n/(r*r)
    for i in range(19):
        metric[i,i] += 1
    heavy = acb_mat([[int(i==j)-int(i==18) for j in range(18)] for i in range(19)])
    for name,basis,lo,hi in (
        ("all",None,arb(".358"),arb(".367")),
        ("heavy",heavy,arb(".185"),arb(".190"))):
        def restrict(m):
            return basis.transpose()*m*basis if basis is not None else m
        assert not positive_definite(restrict(lo*metric-lower))
        assert positive_definite(restrict(hi*metric-upper))
        print("Relative transverse insertion / (g_D^2 T), Q=mu",name,lo.str(6),hi.str(6))


if __name__ == "__main__":
    verify()
