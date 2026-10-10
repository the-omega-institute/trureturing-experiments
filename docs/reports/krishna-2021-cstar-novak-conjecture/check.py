# Exact check of a counterexample to Krishna's C*-algebraic Novak conjecture (arXiv:2108.06662, Conjecture 4.3; JKMS 59 (2022) Conjecture 4.2)
# in A = M_2(C), n = 3, any d >= 2 (the same coordinate in every position l), with cos x := (exp(ix) + exp(-ix))/2.
import sympy as sp
pi=sp.pi; r15=sp.sqrt(15); I2=sp.eye(2)
A=[3*pi*I2, pi*sp.Matrix([[5,0],[0,1]]), (pi/4)*sp.Matrix([[19,r15],[r15,5]])]
for X in A: assert X.is_hermitian
def cosm(X):
    E=lambda Y: sp.simplify((Y).exp())
    return sp.simplify((E(sp.I*X)+E(-sp.I*X))/2)
C=[[cosm(A[j]-A[k]) for k in range(3)] for j in range(3)]
for j in range(3):
    for k in range(3): print(j,k,'(A_j-A_k)^2 =',sp.simplify((A[j]-A[k])**2).tolist(),' cos =',C[j][k].tolist())
n=3
for d in (2,3,5):
    N=[[sp.simplify(((I2+C[j][k])/2)**d - I2/n) for k in range(n)] for j in range(n)]
    v=[-2,1,1]
    # <N x, x> with x_k = v_k * 1_A and the standard A-valued inner product sum_j y_j x_j^*
    form=sp.simplify(sum((N[j][k]*v[k]*v[j] for j in range(n) for k in range(n)),sp.zeros(2)))
    blocks=sp.Matrix([[N[j][k][0,0] for k in range(n)] for j in range(n)])
    print('d =',d,' scalar block of N =',blocks.tolist(),' <Nx,x> =',form.tolist())
    assert form==-2*I2
print('OK')
