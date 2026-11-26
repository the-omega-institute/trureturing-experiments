"""Unique least-degree tetrahedral fixed-ray update and its finite section count."""
from sympy import Matrix, Poly, Rational as Q, cancel, diff, eye, factor, gcd, symbols

z,w,t,a,omega = symbols("z w t a omega")


def main():
    # Orthogonal projections of the four coordinate vectors to the sum-zero hyperplane.
    projection = eye(4)-Matrix.ones(4,4)/4
    assert projection**2==projection
    assert all(projection[i,j]==(Q(3,4) if i==j else -Q(1,4))
               for i in range(4) for j in range(4))
    metric = Matrix.diag(1,2)
    half_turn = Matrix([[1,2],[1,-1]])
    assert half_turn.T*metric*half_turn==3*metric
    quartic = w*(z**3-w**3)
    P,S = -diff(quartic,w),diff(quartic,z)
    assert gcd(P,S)==1
    assert (z*S-w*P-4*quartic).expand()==0
    R = (4-t**3)/(3*t**2)
    M = (t+2)/(t-1)
    assert cancel(R.subs(t,M)-M.subs(t,R))==0
    num = cancel(R.subs(t,omega*t)-omega*R).as_numer_denom()[0]
    assert Poly(num,omega).rem(Poly(omega**2+omega+1,omega)).is_zero
    # For fixed pole [1:0] and degree 3, rotation covariance restricts the
    # numerator to t^3 and 1, and the denominator to a multiple of t^2.
    family = a*t+(1-a)/t**2
    assert cancel(family.subs(t,1))==1
    pole_chart = cancel(1/family.subs(t,1/z))
    assert cancel(diff(pole_chart,z).subs(z,0)-1/a)==0
    assert cancel(diff(family,t).subs(t,1)-(3*a-2))==0
    assert factor(a*(3*a-2)-1)==(a-1)*(3*a+1)
    assert cancel(family.subs(a,-Q(1,3))-R)==0
    assert cancel(diff(R,t).subs(t,1))==-3
    assert (P-S+(z-w)*(z+2*w)**2).expand()==0
    assert S==3*z**2*w
    # For target [0:1] all three points are simple and have nonzero denominator.
    assert gcd(4-t**3,3*t**2)==1
    assert gcd(4-t**3,diff(4-t**3,t))==1
    # Degree-two homogeneous sections of O(2); optional finite sampling frame.
    eval3 = Matrix([[1,k,k*k] for k in (0,1,-1)])
    assert abs(eval3.det())==2
    print("Four tetrahedral fixed rays, least nonidentity degree 3; unique update",
          "[4w^3-z^3:3z^2w]; twisted holomorphic spin sections have basis w², zw, z².")


if __name__=="__main__":
    main()
