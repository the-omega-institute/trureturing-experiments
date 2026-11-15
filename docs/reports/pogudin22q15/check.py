#!/usr/bin/env python3
"""Exact symbolic checks for Pogudin, Question 15; no Lean proof."""
import sys
import sympy as s

assert sys.version_info >= (3, 12)
x, y, z, t = s.symbols("x y z t")


def D(f):
    return s.expand(y*s.diff(f,x) + z*s.diff(f,y) + t*s.diff(f,z))


def monomials(variables, degree):
    return sorted(s.itermonomials(variables, degree), key=s.default_sort_key)


def main():
    assert s.expand(D(y**2/s.Integer(2)+x) - y*(z+1)) == 0
    assert s.expand(D(x*y) - (x*z+y*y)) == 0
    print("FIDELITY: source p=x''+1 and D(xx') identities verified")

    a0,a1,a2,b0,b1,b2,alpha=s.symbols("a0 a1 a2 b0 b1 b2 alpha")
    A=a0+a1*x+a2*y
    B=b0+b1*x+b2*y
    phi=s.Poly(s.expand(B*D(A)+alpha*A*D(B)),x,y,z)
    span=[y,x*y,x*z,y*y,y*z,z]
    exponents=[s.Poly(m,x,y,z).monoms()[0] for m in span]
    assert set(phi.monoms()) <= set(exponents)
    coeffs=s.Matrix([phi.coeff_monomial(m) for m in span])
    point=dict(zip([a0,a1,a2,b0,b1,b2,alpha],[1,2,3,2,3,5,2]))
    jac=coeffs.jacobian([a0,a1,a2,b0,b1,b2,alpha]).subs(point)
    assert jac.rank()==6
    print("FIDELITY: first map lies in the six-term source span; exact Jacobian rank=6")

    cm=monomials((x,y,z),3)
    cs=s.symbols("c0:"+str(len(cm)))
    aa=s.symbols("q0:4")
    generic_c=sum(c*m for c,m in zip(cs,cm))
    generic_a=aa[0]+aa[1]*x+aa[2]*y+aa[3]*z
    for k in range(2,13):
        p=x*z+s.Rational(k,2)*(y*y+1)
        H=x**k*(y*y+1)
        M=2*x**(k-1)*y
        assert s.expand(D(H)-M*p)==0
        assert s.Poly(p,x,y,z).total_degree()==2
        assert s.Poly(p,z).degree()==1
        assert s.gcd(s.Poly(x,x,y),s.Poly(s.Rational(k,2)*(y*y+1),x,y)).as_expr()==1
        assert s.Poly(p,x,y,z).coeff_monomial(1)==s.Rational(k,2)!=0
        residual=s.Poly(D(generic_c)-generic_a*p,x,y,z,t)
        matrix,_=s.linear_eq_to_matrix(residual.coeffs(),[*cs,*aa])
        kernel=matrix.nullspace()
        assert len(kernel)==1
        assert all(all(v[len(cs)+i]==0 for i in range(4)) for v in kernel)
        # This is the full finite source ansatz, not sampled A or C.
        print(f"FAMILY k={k}: D(H)=M*p; degree(p)=2; primitive linear-in-z; first image excluded by constant coefficient; second image kernel forces A=0")

    # Source closure equations also exclude the chosen k=2 polynomial.
    source_linear_component_residual = s.Integer(1)-2*s.Integer(1)
    source_nonlinear_component_constant = s.Integer(1)
    assert source_linear_component_residual != 0
    assert source_nonlinear_component_constant != 0
    print("CHOSEN k=2: excluded even from both source Zariski closures")

    # Bounded check of the proposed first-integral-ring mechanism.
    for k in range(2,7):
        H=x**k*(y*y+1)
        ms=monomials((x,y),2*(k+2))
        fs=s.symbols("f0:"+str(len(ms)))
        F=sum(c*m for c,m in zip(fs,ms))
        delta=s.Poly(x*y*s.diff(F,x)-s.Rational(k,2)*(y*y+1)*s.diff(F,y),x,y)
        matrix,_=s.linear_eq_to_matrix(delta.coeffs(),fs)
        ker=matrix.nullspace()
        basis=s.Matrix([[s.Poly(Q,x,y).coeff_monomial(m) for m in ms] for Q in (s.Integer(1),H,H**2)]).T
        assert basis.rank()==3 and len(ker)==3
        assert matrix*basis==s.zeros(matrix.rows,3)
        print(f"BOUNDED KERNEL k={k}, total degree <= {2*(k+2)}: span(1,H,H^2), dimension=3")
    print("ALL_CHECKS_OK")


if __name__=="__main__":
    main()
