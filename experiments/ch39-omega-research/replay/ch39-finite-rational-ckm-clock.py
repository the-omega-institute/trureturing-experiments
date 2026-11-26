"""只用有限整数与高斯有理矩阵重放最短三段时钟。在 --source 下核对零点区间来源。"""
import sys
from math import factorial
from pathlib import Path
from runpy import run_path

from sympy import I, Matrix, Rational, eye, expand, zeros

# 按列向量的实际演化先后排列。前两项是 X,Z 的 80 位小数，后四项是
# E*H E 上两个对角及上三角元的实、虚部的 18 位小数，H 为逆凯莱坐标。
records = (
    (11410,10985,(3916189176622806672740,3916189179731174079171),(181839,193406277,8660,2643)),
    (11431,676974,(0,0),(-101968456,-12530752429,48213349,41102986)),
    (11432,327502,(-131355955457923331743,-131355955454426259888),(6910008,8689838409,1776439,1434998)),
)
basis = Matrix([[5,29],[-1,145],[0,-26]])
radial = Matrix([1,5,29])


def verify_rational():
    dual_basis = basis*(basis.T*basis).inv()
    wave = eye(3)
    for _,_,_,(a,b,x,y) in records:
        block = Matrix([[a,x+I*y],[x-I*y,b]])/2**18
        generator = dual_basis*block*dual_basis.T
        gate = ((eye(3)+I*generator)*(eye(3)-I*generator).inv()).applyfunc(expand)
        assert (gate.H*gate-eye(3)).applyfunc(expand)==zeros(3)
        assert (gate*radial-radial).applyfunc(expand)==zeros(3,1)
        wave = (gate*wave).applyfunc(expand)
    assert (wave.H*wave-eye(3)).applyfunc(expand)==zeros(3)
    a,b,u = [expand(wave[i,j]*wave[i,j].conjugate()) for i,j in ((0,1),(1,2),(0,2))]
    for v,lo,hi in ((a/(1-u),"0.224955","0.225065"),(b/(1-u),"0.041775","0.041885"),
                    (u,"0.003677","0.003787")):
        assert Rational(lo)**2<v<Rational(hi)**2
    phase = expand(wave[0,1]*wave[1,2]*wave[0,2].conjugate()*wave[1,1].conjugate()+a*b*u/(1-u))
    real,imag = phase.as_real_imag()
    assert real>0 and imag>0
    # sin、cos 的交错项相继减小。16、17 项分别给下界和上界。
    def tangent_bounds(x):
        sine = sum((-1)**k*x**(2*k+1)/factorial(2*k+1) for k in range(16))
        cosine = sum((-1)**k*x**(2*k)/factorial(2*k) for k in range(16))
        assert cosine>0
        return sine/(cosine+x**32/factorial(32)), (sine+x**33/factorial(33))/cosine
    assert tangent_bounds(Rational("1.141"))[1] < imag/real < tangent_bounds(Rational("1.153"))[0]
    total = sum(row[1] for row in records)
    assert total==1015461 and 3*records[1][1]==2*total
    for k in (0,1):
        mean = sum(Rational(row[1]*row[2][k],total*2**80) for row in records)
        assert abs(mean)+Rational(1,2**81)<Rational("2.18768e-14")
    print("Three Gaussian-rational unitaries, radial axis, four CKM windows and two budgets verified.")


def verify_source():
    from flint import acb, acb_mat, arb
    ordered = run_path(str(Path(__file__).with_name("ch39-ordered-flavor-unitary.py")))
    cost,_,scale,*_ = ordered["previous"]["setup"]()
    flavor,study = ordered["flavor"],ordered["study"]
    coeff = ordered["previous"]["box"]["centers"]
    b = flavor["basis"]
    e = acb_mat([[int(v) for v in row] for row in basis.tolist()])
    ident = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    for n,ticks,readings,integers in records:
        height = acb.zeta_zero(n).imag
        response = study["native"]["plane_response"](height)
        x = study["center"]/(height**2+study["offset"])-1 if n!=11431 else arb(0)
        z = x+scale*cost(x)[0] if n!=11431 else arb(0)
        wave=flavor["radial"]+b*ordered["evolution"](response,arb(ticks)/1015461,coeff)*b.transpose()
        generator=acb(0,-1)*(wave-ident)*(wave+ident).inv()
        block=e.transpose()*generator*e
        entries=(block[0,0].real,block[1,1].real,block[0,1].real,block[0,1].imag)
        for values,stored,digits in (((x,z),readings,80),(entries,integers,18)):
            for value,integer in zip(values,stored):
                assert abs(value*2**digits-integer)<arb(1)/2
    print("All dyadic inputs match the isolated zero data at the central flavor coefficients.")


if __name__=="__main__":
    verify_rational()
    if "--source" in sys.argv:
        verify_source()
