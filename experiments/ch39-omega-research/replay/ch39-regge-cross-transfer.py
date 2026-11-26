"""交叉 Regge 偏振的两条时间传递根与八片有理反射检验。"""
from pathlib import Path
from runpy import run_path

from flint import fmpq, fmpq_mat
import sympy as sp


def verify():
    root = Path(__file__).parent
    regge = run_path(str(root/"ch39-regge-transverse-polarization.py"))
    T, V = sp.symbols("T V")
    squared = sp.sqrt(V**2+96)
    a, b = (-V+squared)/4, (-V-squared)/4
    kinetic = -regge["cross_symbol"]()
    assert sp.simplify(kinetic-(a-T)*(T-b)/10) == 0
    pieces = 10/(a-b)*(1/(a-T)+1/(T-b))
    assert sp.simplify(kinetic*pieces-1) == 0
    assert a.subs(V,2)==2 and b.subs(V,2)==-3
    assert a.subs(V,-2)==3 and b.subs(V,-2)==-2
    # 对任意两个实传递根，半片反射的二阶 Hankel 子式可直接因式分解。
    x,y,c,d = sp.symbols("x y c d")
    moment = lambda n: c*x**n+d*y**n
    assert sp.factor(moment(1)*moment(3)-moment(2)**2) == c*d*x*y*(x-y)**2
    n = 8
    shift = fmpq_mat([[int((i+1)%n==j) for j in range(n)] for i in range(n)])
    identity = shift**0
    time = shift+shift.transpose()
    matrix = (12*identity-2*time*time)/20
    covariance = matrix.inv()
    assert [covariance[0,j] for j in range(5)] == [
        fmpq(35,12),0,fmpq(5,6),0,fmpq(5,12)]
    half = fmpq_mat([[covariance[0,i+j+1] for j in range(2)] for i in range(2)])
    onsite = fmpq_mat([[covariance[0,i+j] for j in range(2)] for i in range(2)])
    assert half.det()==-fmpq(25,36)
    assert (half[0,0]+half[1,1]-2*half[0,1]) == -fmpq(5,3)
    assert onsite[0,0]>0 and onsite.det()>0
    print("八片周期格的前三个协方差、半片反射子式及两片之差的反射二次型",
          covariance[0,0],covariance[0,1],covariance[0,2],half.det(),-fmpq(5,3))


if __name__ == "__main__":
    verify()
