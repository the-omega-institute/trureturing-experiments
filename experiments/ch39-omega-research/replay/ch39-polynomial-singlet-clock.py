"""以一个局部多项式接触项和中性链制备非阿贝尔单态。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
from sympy import cancel, symbols

root = Path(__file__).parent
q = arb


def verify():
    s,t,c = symbols("s t c", real=True)
    contact = s*t-c + ((s-1)**2+(t-1)**2)/2
    difference = s+t-1-c
    assert cancel(contact-difference-(s+t-2)**2/2)==0
    matrix = run_path(str(root/"ch39-quantum-matrix-weak-shell.py"))
    count_power,label_power,e,reverse = matrix["verify"](singlet_clocks=True)
    ctx.prec = 128
    length = 23226
    h_min = q("5e-18")
    heat_low,heat = q("1.25e-18"),q("5e-18")
    box = q(".05")
    # 每个矩阵相对酉矩阵中心的 L²、L⁴ 范数，单位 sqrt(H)。
    sigma = q(".25").sqrt()
    d2 = 2*box+2*sigma
    d4 = 2*box+q(20).root(4)*sigma
    path2 = 2*d2+q(2).sqrt()*d2**2*heat.sqrt()
    path4 = 2*d4+q(2).sqrt()*d4**2*heat.sqrt()
    gram = 2*d2+q(2).sqrt()*d4**2*heat.sqrt()
    # 控制两个路径的外积、内积以及单条路径模长平方的偏差。
    scalar = 2*path2+path2**2*heat.sqrt()
    modulus = 2*path2+path4**2*heat.sqrt()
    radial_cost = 15*6*gram**2
    closure_cost = 90*path2**2
    contact_cost = 3*(scalar**2+modulus**2)
    chain_cost = 3*((1+q(1)/length)*scalar**2+2*(e/(length*heat_low)).sqrt()*scalar)
    classical_gap = 3*e**2/(8*heat_low)
    lower_cost = 960+radial_cost+closure_cost+contact_cost+chain_cost+classical_gap
    span = 15*q(2).log()-60*(box**2*heat_low).log()+lower_cost
    energy = 3*e+2*heat*(span+60*q(2).log()+6*(length-1)*(2*q(2).cosh()).log())
    charge = (4*energy*(q(45)/2+energy/8)/h_min**2).sqrt()
    radius = 1/(28*charge)
    assert 0 < radius < q.pi()/2
    loss = q(2).log()-14*(8*radius**3/(3*q.pi()**3)).log()
    slope = 3*(reverse-e)/(2*(1+e/8)*(1+reverse/8))
    print("多项式单态链的错计次、错编号指数下界",(count_power-span-loss).str(12),
          (label_power-span-loss).str(12),"目标取向均值缺额上界",(2*heat*(span+loss)/slope).str(12))
    assert count_power-span-loss>13500 and label_power-span-loss>790000
    assert 2*heat*(span+loss)/slope < q("3.5e-10")
    assert (1+e/8).sqrt()-1<q("2e-15")


if __name__ == "__main__":
    verify()
