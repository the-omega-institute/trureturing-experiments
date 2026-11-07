#!/usr/bin/env python3
"""Exact checks for the singleton-loop subproblem of arXiv:2604.21552v1.

All matrix entries, polynomial coefficients, ranks and determinants are exact.
Finite checks support, but do not prove, the quantified classification.
"""
import sympy as S

x = S.Symbol('x')

def path(n, k, upper=None, lower=None):
    """One-based loop position; exact support, arbitrary directed edge weights."""
    A = S.zeros(n)
    A[k - 1, k - 1] = 1
    for i in range(n - 1):
        A[i, i + 1] = 1 if upper is None else upper[i]
        A[i + 1, i] = 1 if lower is None else lower[i]
    return A

def zero(M):
    return all(S.simplify(e) == 0 for e in M)

def supported(A, k):
    return all((A[i,j] != 0) == (abs(i-j) == 1 or i == j == k-1)
               for i in range(A.rows) for j in range(A.cols))

def forbidden(A, X):
    assert X != S.zeros(A.rows)
    assert zero(A.multiply_elementwise(X))
    assert zero(A * X.T - X.T * A)
    assert not any(e.has(S.Float) for e in list(A) + list(X))

def verification_rank(A):
    """Full arbitrary-X commutator map, after eliminating supported X entries."""
    n = A.rows
    columns = []
    for i in range(n):
        for j in range(n):
            if A[i,j] == 0:
                E = S.zeros(n); E[i,j] = 1
                columns.append(S.Matrix(list(A * E.T - E.T * A)))
    return (S.Matrix.hstack(*columns).rank() if columns else 0, len(columns))

def continuant(q, n):
    p0, p1 = S.Poly(1,x), S.Poly(x,x)
    if n == 0:
        return p0
    for i in range(n - 1):
        p0, p1 = p1, S.Poly(x,x) * p1 - q * p0
    return p1

def coefficient_jacobian(q, n, k):
    """Differentiate the full tridiagonal recurrence in n-1 products and loop a."""
    Z, one = S.Poly(0,x), S.Poly(1,x)
    p0, p1 = one, S.Poly(x - (1 if k == 1 else 0),x)
    d0 = [Z for _ in range(n)]
    d1 = [-one if j == n-1 and k == 1 else Z for j in range(n)]
    for vertex in range(1,n):
        factor = S.Poly(x - (1 if vertex == k-1 else 0),x)
        p2 = factor*p1-q[vertex-1]*p0
        d2 = [factor*d1[j]-q[vertex-1]*d0[j]
              - (p0 if j == vertex-1 else Z)
              - (p1 if j == n-1 and vertex == k-1 else Z)
              for j in range(n)]
        p0,p1,d0,d1 = p1,p2,d1,d2
    return S.Matrix(n,n,lambda i,j:d1[j].nth(i)), p1

def fidelity():
    A6 = S.Matrix([[0,1,0,0,0,0],[1,0,-1,0,0,0],[0,1,1,1,0,0],
                   [0,0,1,0,1,0],[0,0,0,-1,0,1],[0,0,0,0,2,0]])
    X6 = S.Matrix([[-1,0,1,-1,2,-2],[0,-2,0,3,-1,4],[-1,0,0,0,-1,0],
                   [1,-3,0,2,0,2],[2,-1,1,0,-1,0],[-1,2,0,-1,0,-2]])
    assert supported(A6,3)
    forbidden(A6,X6)
    rank6,cols6=verification_rank(path(6,3)); assert rank6 == cols6 == 25
    print('FIDELITY Example 5.9: displayed A,X satisfy both equations; adjacency rank=25/25')
    A5 = S.Matrix([[0,1,0,0,0],[1,0,-1,0,0],[0,1,1,1,0],
                   [0,0,1,0,-1],[0,0,0,1,0]])
    assert supported(A5,3)
    rank5,cols5=verification_rank(A5); assert rank5 == cols5 == 16
    print('FIDELITY Example 5.10: displayed A verification rank=16/16')

def allowed_cases():
    count=0
    for n in range(1,13):
        for k in range(1,n+1):
            r,s=k-1,n-k
            if r%2 == s%2 == 1:
                continue
            for t in range(1,1001):
                P,Q=continuant(1,r),continuant(t,s)
                if S.gcd(P,Q).degree() == 0:
                    break
            else:
                raise AssertionError(('coprime search exhausted',n,k))
            # Edges indexed zero-based: left internals precede r-1, right internals follow r.
            q=[S.Integer(t if i > r else 1) for i in range(n-1)]
            J,F=coefficient_jacobian(q,n,k)
            det=J.det(); assert det != 0
            # Direct characteristic and arbitrary-X checks independently of the recurrence/Jacobian.
            A=path(n,k,lower=q); assert supported(A,k)
            assert S.Poly(A.charpoly(x).as_expr(),x) == F
            if n <= 7:
                rank,cols=verification_rank(A); assert rank == cols
            print(f'ALLOW n={n} k={k} right_internal_product={t} jacobian_det={det}')
            count+=1
    assert count == 63
    print('ALLOW_SUMMARY: 63/63 exact nonzero determinants, n=1..12; direct arbitrary-X ranks n<=7')

def obstruction_cases():
    count=0
    for n in range(3,12,2):
        for k in range(2,n+1,2):
            A=path(n,k,upper=[(-1)**i*(i+1) for i in range(n-1)],
                       lower=[i+2 for i in range(n-1)])
            v,w=S.zeros(n,1),S.zeros(n,1);v[0]=w[0]=1
            for i in range(1,n-1,2):
                v[i+1]=-A[i,i-1]*v[i-1]/A[i,i+1]
                w[i+1]=-A[i-1,i]*w[i-1]/A[i+1,i]
            assert zero(A*v) and zero(A.T*w)
            forbidden(A,(v*w.T).T); assert supported(A,k)
            print(f'OBSTRUCTION n={n} k={k} X11=1')
            count+=1
    assert count == 15
    print('OBSTRUCTION_SUMMARY: 15/15 exact nonsymmetric signed-weight witnesses, odd n<=11, even k')

def unit_eigen_branch(m):
    """Positive symmetric branch with eigenvalue 1, using only exact radicals."""
    b=[S.Rational(1,2)]*(m-2)+[S.sqrt(S.Rational(m,2*(m-1)))]
    B=S.zeros(m)
    for i,t in enumerate(b): B[i,i+1]=B[i+1,i]=t
    v=S.zeros(m,1);v[0]=1;v[1]=1/b[0]
    for i in range(1,m-1):v[i+1]=S.simplify((v[i]-b[i-1]*v[i-1])/b[i])
    assert zero(B*v-v) and all(t.is_positive for t in v)
    return B,v

def nonrequire_cases():
    branches={m:unit_eigen_branch(m) for m in range(2,7)}
    count=0
    for r in range(2,7):
        for s in range(2,7):
            B,v=branches[r];D,w=branches[s]
            n,k=r+s+1,r+1
            A=S.diag(B,S.ones(1),D)
            A[r-1,r]=A[r,r-1]=A[r,r+1]=A[r+1,r]=1
            z=S.Matrix(list(w[0]*v)+[0]+list(-v[r-1]*w))
            alt=S.diag(*[(-1)**i for i in range(n)]); u=alt*z
            assert zero(A*z-z) and zero(A*u+u)
            C=z*z.T+u*u.T
            assert supported(A,k);forbidden(A,C)
            print(f'NONREQUIRE r={r} s={s} exact_radical_commutant=OK')
            count+=1
    assert count == 25
    print('NONREQUIRE_SUMMARY: 25/25 connected exact-support symmetric witnesses, r,s=2..6')

if __name__ == '__main__':
    print('Exact arithmetic: SymPy '+S.__version__+'; no floating point evidence')
    fidelity()
    allowed_cases()
    obstruction_cases()
    nonrequire_cases()
    print('ALL_CHECKS_OK')
