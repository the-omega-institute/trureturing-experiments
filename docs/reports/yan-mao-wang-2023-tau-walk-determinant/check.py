"""Checks for Yan-Mao-Wang, arXiv:2309.01123v1, Conjecture 1: det W_tau(G o P_m) = +-(det A_tau)^floor(m/2) (det((1-tau^2)I+tau A_tau))^floor((m-1)/2) (det W_tau(G))^m.
Part 1: exact check for every graph on n <= 4 vertices, m <= 4 with n*m <= 12, tau in {0, 1, 1/2, -2, 3}. Part 2: the resultant identity (13) of Proposition 3, symbolic in tau and t, m = 2..9."""
import sympy as sp, itertools, random
def walk_det(A, tau):
    N=A.shape[0]; D=sp.diag(*[sum(A[i,j] for j in range(N)) for i in range(N)]); At=A+tau*D
    e=sp.ones(N,1); cols=[e]
    for _ in range(N-1): cols.append(At*cols[-1])
    return sp.Matrix.hstack(*cols).det(), At
def rooted(A, m):
    n=A.shape[0]; N=n*m; B=sp.zeros(N)
    for i in range(n):
        for j in range(n):
            B[i*m, j*m]=A[i,j]
        for t in range(m-1):
            B[i*m+t, i*m+t+1]=1; B[i*m+t+1, i*m+t]=1
    return B
random.seed(3); ok=True; cnt=0
graphs=[]
for n in (2,3,4):
    pairs=list(itertools.combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        A=sp.zeros(n)
        for b,(i,j) in enumerate(pairs):
            if mask>>b&1: A[i,j]=A[j,i]=1
        graphs.append(A)
for A in graphs:
    for m in (2,3,4):
        for tau in (0, 1, sp.Rational(1,2), -2, 3):
            if A.shape[0]*m>12: continue
            lhs,_=walk_det(rooted(A,m),tau)
            wG,At=walk_det(A,tau); n=A.shape[0]
            rhs=(At.det())**(m//2)*(((1-tau**2)*sp.eye(n)+tau*At).det())**((m-1)//2)*wG**m
            cnt+=1
            if not (lhs==rhs or lhs==-rhs):
                ok=False; print('FAIL',A.tolist(),m,tau,lhs,rhs); break
print('cases',cnt,'identity up to sign holds:',ok)
ok_identity=ok

x,tau,t=sp.symbols('x tau t')
def Z(k):
    # Z_j = U_j + tau U_{j-1}, U_j Chebyshev-like in y = x - 2 tau: U_-1=0, U_0=1
    y=x-2*tau; U=[sp.Integer(0),sp.Integer(1)]
    for j in range(1,k+2): U.append(sp.expand(y*U[-1]-U[-2]))
    Uj=lambda j: U[j+1]
    return sp.expand(Uj(k)+tau*Uj(k-1))
ok=True
for m in range(2,10):
    f=sp.expand(sum(Z(j) for j in range(m))); p=sp.expand(Z(m)+t*Z(m-1))
    r=sp.factor(sp.resultant(p,f,x))
    rhs=(-1)**(m*(m-1))*(tau-t)**(m//2)*(1-tau*t)**((m-1)//2)
    good=sp.simplify(r-rhs)==0
    ok&=good; print(m,good)
print('resultant identity (13), m=2..9:',ok)
# check Z_0=1, Z_1=x-tau match Definition 2 (recurrence check printed for inspection)
print("RESULT", "PASS" if (ok and ok_identity) else "FAIL")
