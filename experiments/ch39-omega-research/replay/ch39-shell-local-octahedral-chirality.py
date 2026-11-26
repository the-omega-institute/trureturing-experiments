"""A charge-one, nearest-neighbor Ginsparg-Wilson operator on six octahedral vertices."""
import sympy as sp


def make_operator():
    vectors = [(0,0,1),(0,0,-1),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0)]
    q = sp.sqrt(2)
    spinors = [(sp.S.One,sp.S.Zero),(sp.S.Zero,sp.S.One),*((1/q,w/q) for w in (1,-1,sp.I,-sp.I))]
    A = sp.Matrix(vectors)/q
    B = sp.Matrix([[(a*a-b*b)/2,sp.I*(a*a+b*b)/2,-a*b] for a,b in spinors])
    assert A.H*A==B.H*B==sp.eye(3)
    assert sp.ones(1,6)*A==sp.zeros(1,3)
    P0 = sp.ones(6)/6
    unit = sp.eye(6)
    mix = -sp.Rational(1,3)
    cross = sp.sqrt(1-mix**2)*A*B.H
    plus = 2*P0-unit+(1+mix)*A*A.H
    minus = unit-(1+mix)*B*B.H
    H = plus.row_join(cross).col_join(cross.H.row_join(minus))
    gamma = sp.diag(*([1]*6+[-1]*6))
    H = H.applyfunc(sp.expand)
    return vectors,H,gamma


if __name__ == "__main__":
    vectors,H,gamma = make_operator()
    eye = sp.eye(12)
    assert H.H==H and (H*H-eye).applyfunc(sp.expand)==sp.zeros(12)
    D = eye-gamma*H
    assert (D*gamma+gamma*D-D*gamma*D).applyfunc(sp.expand)==sp.zeros(12)
    assert (gamma+H).trace()/2==1
    assert D*sp.Matrix([1]*6+[0]*6)==sp.zeros(12,1)
    for i in range(6):
        j = i^1
        assert all(D[i+6*a,j+6*b]==0 for a in (0,1) for b in (0,1))
    parameter=sp.Symbol("z")
    squared=(D.H*D).applyfunc(sp.expand)
    assert sp.factor(squared.charpoly(parameter).as_expr())==sp.factor(
        parameter*(parameter-sp.Rational(8,3))**6*(parameter-4)**5)
    print("one positive zero, squared singular spectrum 0:1, 8/3:6, 4:5")
    print("all six directed antipodal 2-by-2 blocks vanish exactly")
