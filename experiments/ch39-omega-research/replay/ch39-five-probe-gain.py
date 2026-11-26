"""Rational primal and dual bounds for the best exact, phase-uniform five-probe effect."""
from itertools import product
from sympy import Matrix, Poly, Rational as Q, diag, eye, symbols


G = diag(1,3,3,1)
COEFF = Matrix([18,-19,-3,4])


def positive(A):
    return A == A.T and all(A[:k,:k].det() > 0 for k in range(1,5))


def fourier(A, k):
    return sum(2**j*A[j,j+k] for j in range(4-k))


def factors():
    z = symbols("z")
    rank_one = []
    interior = G*0
    for r, flip in product(range(3), range(2)):
        first = (z-1)**(2-r)*(z-2)**r
        last = 9*z+8 if flip else 4*z+9
        co = Matrix(Poly(first*last, z).all_coeffs()[::-1])
        factor = Q(2)**(1-r-flip) * co*co.T
        assert all(fourier(factor,k) == fourier(COEFF*COEFF.T,k) for k in range(4))
        rank_one.append((G.inv()*factor).trace())
        interior += factor/6
    assert min(rank_one) == Q(584,3)
    assert positive(interior)
    return interior


def verify():
    interior = factors()
    moments = [fourier(COEFF*COEFF.T, k) for k in range(4)]
    # Six free entries are rounded. Reconstruct the first row from the four exact moments.
    X = Matrix([[0,0,0,0], [0,106374339072,32082853105,-62671525319],
                [0,32082853105,99671836673,-68930674898],
                [0,-62671525319,-68930674898,64734660057]])/10**9
    for k,b in enumerate(moments):
        X[0,k] = X[k,0] = b-fourier(X,k)
    X = X + (interior-X)/10**8
    upper = Q("189.60582")
    assert all(fourier(X,k) == b for k,b in enumerate(moments))
    assert positive(X) and positive(upper*G-X)

    y = -Matrix([353570004591,996097826459,1506967976825,-445376643199])/10**12
    Y = Matrix([[423903524977,-120116872788,-136972356379,380316786232],
                [-120116872788,34036195218,38812347948,-107766178718],
                [-136972356379,38812347948,44258717629,-122888542580],
                [380316786232,-107766178718,-122888542580,341211736572]])/10**12
    Y += eye(4)/10**9
    adjoint = Matrix(4,4,lambda i,j: 2**min(i,j)*y[abs(i-j)] / (1 if i==j else 2))
    assert positive(Y) and positive(Y-adjoint)
    lower = y.dot(Matrix(moments))/(G*Y).trace()
    assert Q("189.60581") < lower < upper
    lo_gain, hi_gain = 970/upper, 970/lower
    assert Q("5.115876") < lo_gain < hi_gain < Q("5.115878")
    print("Rank-one optimal gain", 970/Q(584,3))
    print("All effects, strict bounds",float(lo_gain),float(hi_gain))


if __name__ == "__main__":
    verify()
