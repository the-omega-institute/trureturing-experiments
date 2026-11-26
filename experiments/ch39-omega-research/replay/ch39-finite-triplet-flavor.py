"""Finite local multiplication vertices on the three chiral zero modes of the permutation shell."""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from sympy import Matrix, I, Rational as Q, diag, eye, simplify, sqrt, zeros

shell = run_path(str(Path(__file__).with_name("ch39-finite-shell-first-triplet.py")))


def cs(t):
    return (1-t*t)/(1+t*t),2*t/(1+t*t)


def main():
    alpha = (sqrt(5)-1)/2
    seed = (1,alpha,-alpha,-1)
    pts = list(permutations(range(4)))
    V = Matrix([[seed[k] for k in p[:3]] for p in pts])
    L = shell["matrix"]({1:1,2:1,4:1})
    assert (L*V-(3-sqrt(5))*V).applyfunc(simplify)==zeros(24,3)
    selected = (0,1,2,3,4,6)
    table = Matrix([[V[k,i]*V[k,j] for k in selected] for i in range(3) for j in range(i,3)])
    assert simplify(table.det(method="domain-ge")+80-36*sqrt(5))==0
    # Real invertible normalization to an orthonormal frame is a congruence
    # on the six-dimensional space of complex symmetric bilinear vertices.
    c12,s12 = cs(Q(57,500))
    c23,s23 = cs(Q(107,5000))
    c13,s13 = cs(Q(377,200000))
    real,imag = cs(Q(16,25))
    phase = real+I*imag
    mix = Matrix([[1,0,0],[0,c23,s23],[0,-s23,c23]])*Matrix(
        [[c13,0,s13*phase.conjugate()],[0,1,0],[-s13*phase,0,c13]])*Matrix(
        [[c12,s12,0],[-s12,c12,0],[0,0,1]])
    assert (mix*mix.H-eye(3)).applyfunc(simplify)==zeros(3)
    J = c12*c23*c13**2*s12*s23*s13*imag
    assert Q(3208,10**8)<J<Q(3210,10**8)
    assert abs(mix[0,2])<Q(38,10000) and abs(mix[2,0])>Q(87,10000)
    up = diag(*(Q(x) for x in (".00000615",".00311",".8616")))
    down = diag(*(Q(x) for x in (".0000135",".000268",".01401")))
    Yu,Yd = up,mix*down*mix.T
    assert Yu==Yu.T and (Yd-Yd.T).applyfunc(simplify)==zeros(3)
    assert (Yd*Yd.H-mix*down**2*mix.H).applyfunc(simplify)==zeros(3)
    print("Six selected site vertices span Sym(3), determinant -80+36 sqrt(5).",
          "Rational-unitary CKM witness J",float(J))


if __name__=="__main__":
    main()
