# SPDX-License-Identifier: Apache-2.0
"""Exact checks for the spacing mechanism; no floating point and no limit certification."""
from fractions import Fraction as F
from math import factorial
import sympy as S


def matmul(a, b):
    return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in zip(*b)] for row in a]


def quadratic(a, c):
    return sum((c[i]*a[i][j]*c[j] for i in range(len(c)) for j in range(len(c))),F(0))


def main():
    # Source, paragraph following Lemma 3.1: uniform-order-statistic inverse covariance.
    for n in range(1, 25):
        covariance=[[F(min(i,j)*(n+1-max(i,j)),(n+1)**2*(n+2)) for j in range(1,n+1)] for i in range(1,n+1)]
        inv=[[(n+1)*(n+2)*(2 if i==j else -1 if abs(i-j)==1 else 0) for j in range(n)] for i in range(n)]
        assert matmul(covariance,inv)==[[int(i==j) for j in range(n)] for i in range(n)]
    print('FIDELITY: uniform covariance inverse, dimensions 1 through 24: exact')

    # Integrating a monomial on the n-dimensional simplex under Dirichlet(1,...,1).
    for n in range(1, 33):
        den=(n+1)*(n+2)
        assert F(factorial(n)*factorial(2),factorial(n+2)) == F(2,den)
        assert F(factorial(n),factorial(n+2)) == F(1,den)
    print('SPACINGS: diagonal 2/((n+1)(n+2)), off-diagonal 1/((n+1)(n+2)), n=1..32: exact')

    packets=0
    for n in range(3, 25):
        second=[[F(min(i,j)*(max(i,j)+1),(n+1)*(n+2)) for j in range(1,n+1)] for i in range(1,n+1)]
        for length in range(2,n,2):
            b=[(-1)**j for j in range(length)]
            for start in range(n-length):
                c=[0]*n
                for j in range(length):
                    c[start+j]-=b[j]
                    c[start+j+1]+=b[j]
                assert sum(b)==sum(c)==0
                assert sum(x*x for x in c)==4*length-2
                assert quadratic(second,c)==F(length,(n+1)*(n+2))
                assert n*n*quadratic(second,c)/sum(x*x for x in c)==F(n*n*length,(n+1)*(n+2)*(4*length-2))
                packets+=1
    print('PACKETS:',packets,'all admissible positions and even lengths for n=3..24: exact')

    n,L=S.symbols('n L',positive=True,integer=True)
    g=S.symbols('g',positive=True)
    rayleigh=n**2*L*g**2/((n+1)*(n+2)*(4*L-2))
    assert S.limit(S.limit(rayleigh,n,S.oo),L,S.oo)==g**2/4
    assert S.simplify((S.sqrt(2*S.pi))**2/4-S.pi/2)==0
    assert S.simplify((S.sqrt(S.pi/2))**2/4-S.pi/8)==0
    for c,target in [(1/S.sqrt(2*S.pi),S.pi/2),(S.sqrt(2/S.pi),S.pi/8)]:
        assert S.simplify(S.limit(n**2/(4*c**2*(n+1)*(n+2)),n,S.oo)-target)==0
        assert S.simplify(S.sqrt(1/target)-(S.sqrt(2/S.pi) if target==S.pi/2 else S.sqrt(8/S.pi)))==0
    print('CONSTANTS: iterated packet bound g^2/4; normal pi/2; half-normal pi/8; distortion factors: exact')
    print('SCOPE: these algebra checks support the outline; Gaussian local L2 convergence and the spectral limit are not certified here')
    print('ALL_CHECKS_OK')


if __name__=='__main__':
    main()
