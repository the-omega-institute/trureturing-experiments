"""两枚相邻零点共用的普朗克有界截止混合所需负量的双侧区间。"""
from pathlib import Path
from runpy import run_path
from math import factorial
import sys

from flint import arb, arb_mat


def small_derivative(x, rows):
    if x.contains(0) or x.upper()>2*x.lower():
        result=arb(0)
        for z,k,w in rows:
            if z.is_zero():
                exponential=arb(1)
            elif x.contains(0):
                exponential=arb(0).union((-z/x.upper()).exp())
            else:
                exponential=(-z/x).exp()
            result+=exponential*(k+w*x)
        return result

    # 在 y=1/x 的中点合并各场的四阶多项式，然后才对区间取范围。
    y=1/x
    mid=y.mid()
    radius=abs(y-mid).upper()
    coeff=[arb(0) for _ in range(5)]
    remainder=arb(0)
    for z,k,w in rows:
        exp0=(-z*mid).exp()
        intercept=exp0*(k*mid+w)
        slope=k*exp0
        coeff[0]+=intercept
        power=arb(1)
        for n in range(1,5):
            coeff[n]+=slope*power
            power*=-z/n
            coeff[n]+=intercept*power
        u=z*radius
        remainder+=(abs(k*y+w).upper()*(-z*y.lower()).exp()*u**5
                    +5*abs(k)*radius*exp0*u**4)/120
    v=coeff[-1]
    offset=arb(0,radius)
    for c in coeff[-2::-1]:
        v=c+v*offset
    return x*(v+arb(0,remainder.upper()))


def long_series(rows, degree):
    coefficients=[sum(w for _,_,w in rows)]
    for n in range(1,degree+1):
        coefficients.append(sum((-z)**(n-1)*(k-w*z/n) for z,k,w in rows)/factorial(n-1))
    remainder=sum(abs(k)*z**degree/factorial(degree)+abs(w)*z**(degree+1)/factorial(degree+1)
                  for z,k,w in rows)
    def evaluate(x):
        v=coefficients[-1]
        for coefficient in coefficients[-2::-1]:
            v=coefficient+v/x
        return x*v+arb(0,remainder/x.lower()**degree)
    return evaluate


def monotone(lo, hi, sign, derivative):
    pending=[(lo,hi)]
    leaves=0
    while pending:
        l,r=pending.pop()
        if sign*derivative(l.union(r))>0:
            leaves+=1
        else:
            assert r-l>arb("1e-22"),(l,r,sign)
            mid=(l+r)/2
            pending.extend(((l,mid),(mid,r)))
    return leaves


def verify(root):
    source = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))
    single, lib, first, second = source["adjacent_fields"](root)
    epsilon=arb("1e-122")
    target=arb_mat([[1],[1],[epsilon],[0],[0]])
    def column(x):
        k,w=single["moments"](first,lib["expint"],x)
        other_k,other_w=single["moments"](second,lib["expint"],x)
        return [arb(1),k,w,other_k-k,other_w-w]

    zero,a,b,c,one=map(arb,("0","3.52557830920e-6",".0274297501576",".541174933638","1"))
    columns=arb_mat([column(x) for x in (zero,a,b,c,one)])
    assert not columns.det().contains(0)
    dual=columns.solve(arb_mat([[-1],[0],[-1],[0],[-1]])).entries()
    _,alpha,beta,gamma,eta=dual
    rows=[(z,field_k*ck/(2*arb.pi()),-field_g*cw/(4*arb.pi()))
          for fields,ck,cw in ((first,alpha-gamma,beta-eta),(second,gamma,eta))
          for field_g,field_k,z in fields]
    def small(x):
        return small_derivative(x,rows)
    large=long_series(rows,19)
    cut=arb("1e-5")
    da,db,dc=map(arb,("2e-14","2e-8","8e-8"))
    pieces=[(zero,a-da,1,small),(a+da,cut,-1,small),(cut,b-db,-1,large),
            (b+db,c-dc,1,large),(c+dc,one,-1,large)]
    leaves=[monotone(*piece) for piece in pieces]
    errors=[width*abs(derivative((x-width).union(x+width))).upper()
            for x,width,derivative in ((a,da,small),(b,db,large),(c,dc,large))]
    e=max(v.upper() for v in errors)
    lower=(dual[0]+alpha+epsilon*beta-e)/(1+2*e)
    assert e<arb("8.209e-14") and lower>arb("337760.30905784")

    # 原子均严格为正，首点的零仅用于上述对偶插值和覆盖闭端。
    sites=(arb("1e-21"),a,b,c,one)
    basis=arb_mat([column(x) for x in sites]).transpose()
    assert not basis.det().contains(0)
    weights=basis.solve(target).entries()
    assert all(sign*p>0 for sign,p in zip((-1,1,-1,1,-1),weights,strict=True))
    upper=-weights[0]-weights[2]-weights[4]
    assert upper<arb("337760.30905793")
    print("单调性认证闭段数、各转折区间最大余量",leaves,*[v.str(12) for v in errors])
    print("两高度共用截止的负系数总量下、上界",lower.str(20),upper.str(20))
    print("五个严格正截止的有符号权",*[v.str(14) for v in weights])


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
