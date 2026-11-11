"""Exact finite linear certificate for the full-kernel Sym_4(R) PB obstruction.

Python >= 3.12; build needs numpy/python-flint; verify needs sympy/python-flint.
No floating arithmetic. All test frames are products of rational Householders.
The certificate proves injectivity of a finite subsystem of necessary linear
conditions, rather than estimating a PB law from samples.
"""
import argparse
import hashlib
import json
import random
import sys
from itertools import combinations
from pathlib import Path

assert sys.version_info >= (3, 12)
ROOT = Path(__file__).resolve().parent
COORDS = [(i, i) for i in range(4)] + list(combinations(range(4), 2))
PAIRS = list(combinations(range(10), 2))
PRIME = 1000003


def digest(rows):
    return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()


def scalar_determinant(rows, p):
    """Separate plain-Python field elimination, without FLINT or NumPy."""
    a=[r[:] for r in rows]; determinant=1
    for k in range(len(a)):
        pivot=next(i for i in range(k,len(a)) if a[i][k] % p)
        if pivot!=k:
            a[k],a[pivot]=a[pivot],a[k]
            determinant=-determinant
        pk=a[k][k]; determinant=determinant*pk%p
        inv=pow(pk,-1,p); tail=a[k][k+1:]
        for i in range(k+1,len(a)):
            factor=a[i][k]*inv%p
            if factor:
                a[i][k+1:]=[(x-factor*y)%p for x,y in zip(a[i][k+1:],tail)]
                a[i][k]=0
    return determinant%p


def build():
    import numpy as np
    from flint import nmod_mat
    p = PRIME
    basis = []
    for i, j in COORDS:
        b = np.zeros((4, 4), dtype=np.int64)
        b[i, j] = b[j, i] = 1
        basis.append(b)

    def hh(a):
        a = np.array(a, dtype=np.int64)
        return (np.eye(4, dtype=np.int64) - 2 * pow(int(a @ a), -1, p) % p * np.outer(a, a)) % p

    def get_rows(test):
        U = hh(test['a']) @ hh(test['b']) % p
        u = np.array(test['u'], dtype=np.int64)
        lam = np.diag([1, 2, 3, 0])
        C = np.array([[2 * pow(i + j + 2, -1, p) % p for j in range(3)] for i in range(3)], dtype=np.int64)
        X = np.zeros((4, 4), dtype=np.int64)
        X[:3, :3] = (np.outer(u, u) * C) % p
        c = U @ lam % p @ U.T % p
        x = U @ X % p @ U.T % p
        w = np.array([int(c[i,j])*int(x[k,l])-int(c[k,l])*int(x[i,j]) for (i,j),(k,l) in [(COORDS[a],COORDS[b]) for a,b in PAIRS]], dtype=np.int64) % p
        w1 = np.array([int(u[1]), -int(u[0]), 0, 0], dtype=np.int64)
        w2 = np.array([int(u[2]), 0, -int(u[0]), 0], dtype=np.int64)
        n = np.array([0,0,0,1], dtype=np.int64)
        vs = [w1,w2,w1+w2,w1+n,w2+n]
        rows = []
        for v0 in vs:
            v = U @ v0 % p
            xv = x @ v % p
            d = np.array([int(xv @ (b @ v % p)) % p for b in basis], dtype=np.int64)
            rows.append((np.outer(w,d) % p).reshape(-1).tolist())
        return rows

    rng = random.Random(731924)
    tests, rows, ranks = [], [], []
    for block in range(12):
        for _ in range(40):
            def normal():
                while True:
                    a = [rng.randint(-3,3) for _ in range(4)]
                    if any(a):
                        return a
            t = {'a':normal(), 'b':normal(), 'u':[rng.choice([-3,-2,-1,1,2,3]) for _ in range(3)]}
            tests.append(t)
            rows.extend(get_rows(t))
        mat = nmod_mat(rows,p)
        rank = mat.rank()
        ranks.append({'tests':len(tests),'rows':len(rows),'rank':rank})
        print(json.dumps(ranks[-1]),flush=True)
        if rank == 450:
            rr, trank = mat.transpose().rref()
            selected = []
            for i in range(trank):
                selected.append(next(j for j in range(rr.ncols()) if int(rr[i,j]) != 0))
            square = [rows[i] for i in selected]
            det = int(nmod_mat(square,p).det())
            assert det != 0
            selected_tests = sorted(set(i//5 for i in selected))
            remap = {old:new for new,old in enumerate(selected_tests)}
            cert = {
                'statement':'The full-kernel compression condition for every legal real c,x forces every coefficient of D to vanish.',
                'prime':p,'dimension':450,'coordinate_basis':COORDS,
                'column_order':'input pairs lexicographic in combinations(range(10),2), then output coordinate 0..9',
                'test_definition':'U=H(a)H(b), H(a)=I-2aa^T/(a^Ta); c=U diag(1,2,3,0) U^T; X_ij=2u_i u_j/(i+j+2) for i,j=0..2; x=UXU^T',
                'kernel_vector_definition':'v=U v0; selectors 0..4: (u2,-u1,0,0),(u3,0,-u1,0),their sum,first+e4,second+e4',
                'row_definition':'(c_a x_b-c_b x_a) * (xv)^T B_gamma v, for a<b and symmetric coordinate basis B_gamma',
                'tests':[tests[i] for i in selected_tests],
                'selected_rows':[[remap[i//5],i%5] for i in selected],
                'determinant_mod_prime':det,'matrix_sha256':digest(square),
                'rank_data':ranks
            }
            (ROOT/'kernel_certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
            return
    raise RuntimeError('Full rank not established; no no-go conclusion')


def verify():
    """Reconstruct selected rows over QQ, then independently reduce modulo p."""
    import sympy as s
    from flint import nmod_mat
    cert = json.loads((ROOT/'kernel_certificate.json').read_text())
    p = cert['prime']
    assert s.isprime(p)
    assert cert['coordinate_basis'] == [list(x) for x in COORDS]
    C = s.Matrix(3,3,lambda i,j:s.Rational(2,i+j+2))
    assert C.inv()*s.ones(3,1) == s.Matrix([6,-30,30])
    assert [C[:i,:i].det() for i in range(1,4)] == [1,s.Rational(1,18),s.Rational(1,5400)]
    cache = {}
    denominator_nonzero = True
    for ti,t in enumerate(cert['tests']):
        def hh(a):
            a = s.Matrix(a)
            return s.eye(4)-2*a*a.T/a.dot(a)
        U = hh(t['a'])*hh(t['b'])
        assert U.T*U == s.eye(4)
        u = s.Matrix(t['u']); du = s.diag(*u)
        X = s.zeros(4); X[:3,:3] = du*C*du
        c = U*s.diag(1,2,3,0)*U.T; x = U*X*U.T
        P = U*s.diag(1,1,1,0)*U.T
        assert P*x*P == x
        r = (c*x+x*c)/2
        u4 = U*s.Matrix([*u,0])
        assert r == u4*u4.T
        w1 = s.Matrix([u[1],-u[0],0,0]); w2 = s.Matrix([u[2],0,-u[0],0]); n = s.eye(4)[:,3]
        vs = [U*v for v in [w1,w2,w1+w2,w1+n,w2+n]]
        assert all(r*v == s.zeros(4,1) for v in vs)
        # Coordinates of c and x are the actual symmetric entries, not traces.
        def mod(q):
            nonlocal denominator_nonzero
            q=s.Rational(q); num,den=map(int,q.as_numer_denom())
            assert den % p != 0
            denominator_nonzero &= den % p != 0
            return (num % p)*pow(den,-1,p)%p
        cc = [mod(c[i,j]) for i,j in COORDS]
        xx = [mod(x[i,j]) for i,j in COORDS]
        ww = [(cc[a]*xx[b]-cc[b]*xx[a])%p for a,b in PAIRS]
        ds = []
        for v in vs:
            z=x*v
            # Direct QQ output functional: z_i v_i, or z_i v_j+z_j v_i.
            ds.append([mod(z[i]*v[j]+(z[j]*v[i] if i!=j else 0)) for i,j in COORDS])
        cache[ti] = (ww,ds)
    square = []
    for ti,selector in cert['selected_rows']:
        ww,ds=cache[ti]; d=ds[selector]
        square.append([a*b%p for a in ww for b in d])
    assert len(square)==450 and all(len(r)==450 for r in square)
    assert digest(square)==cert['matrix_sha256']
    mat=nmod_mat(square,p); rank=mat.rank(); det=int(mat.det())
    scalar_det=scalar_determinant(square,p)
    assert rank==450 and det==cert['determinant_mod_prime'] and det==scalar_det and det!=0
    return {'rational_legality_checks':True,'all_denominators_invertible_mod_prime':bool(denominator_nonzero),'selected_constraints':450,'distinct_rational_frames_and_actors':len(cache),'prime':p,'rank':rank,'determinant_mod_prime':det,'plain_python_determinant_mod_prime':scalar_det,'matrix_sha256':digest(square)}


def fixed_corner_checks():
    import sympy as s
    C=s.Matrix(3,3,lambda i,j:s.Rational(2,i+j+2)); Lam=s.diag(1,2,3)
    K=s.Matrix([[0,1,0],[-1,0,0],[0,0,0]])
    w=s.Matrix([1,-1,0]); z=s.Matrix([0,1,-1])
    def q(y): return (y.T*C*y)[0]
    def f(y): return ((C*y).T*(Lam*K*C-C*K*Lam)*y)[0]
    assert C.inv()*s.ones(3,1)==s.Matrix([6,-30,30])
    assert [(C*y).dot(C*y)/q(y) for y in [w,z]]==[s.Rational(67,75),s.Rational(19,15)]
    assert [f(w),f(z)]==[s.Rational(1151,900),-s.Rational(11,20)]
    # Cyclic jk=(23,31,12) orientation, including its minus sign for (3,1).
    permutation_vectors=[]
    for lam in [(1,2,3),(2,1,3)]:
        cl=s.Matrix(3,3,lambda i,j:s.Rational(2,lam[i]+lam[j])); d=cl.inv()*s.ones(3,1)
        vec=[s.cancel(d[i]*s.Rational(lam[j]+lam[k],2)/(lam[j]-lam[k])) for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]]
        permutation_vectors.append(vec)
    assert permutation_vectors==[[-15,-30,-45],[30,15,45]]
    assert s.Matrix.hstack(*map(s.Matrix,permutation_vectors)).rank()==2
    # Whole 4x4 products for D=0, using one unnormalised positive functional.
    c=s.diag(1,2,3,0); x=s.diag(100,0,0,0); v=s.Matrix([s.Rational(1,10),0,0,1])
    t=s.zeros(4); t[0,3]=t[3,0]=1; t[1,3]=t[3,1]=s.Rational(1,10)
    def jordan(a,b): return (a*b+b*a)/2
    r=jordan(c,x); k=jordan(c,t); L=2*jordan(x,k)-2*jordan(t,r)
    phi=lambda a:(v.T*a*v)[0]
    assert phi(L)==-10
    budgets=[]; readouts=[]; wedges=[]
    for delta in [1,-1]:
        R=x-delta*k; S=t+delta*r
        A=R*v; B=S*v
        budget=phi(R*R)*phi(S*S)-phi(jordan(R,S))**2
        exterior=[A[i]*B[j]-A[j]*B[i] for i,j in combinations(range(4),2)]
        expected=[s.Rational(2)+s.Rational(delta,20),0,s.Rational(3,2),0,-s.Rational(delta,200),0]
        assert exterior==expected
        assert budget==sum(a*a for a in exterior)
        assert budget==(2+s.Rational(delta,20))**2+s.Rational(9,4)+s.Rational(1,40000)
        budgets.append(budget); wedges.append(list(map(str,exterior)))
        readouts.append({'delta':delta,'Rv':list(map(str,A)),'Sv':list(map(str,B)),'phi_R2':str(phi(R*R)),'phi_S2':str(phi(S*S)),'phi_RcircS':str(phi(jordan(R,S)))})
    mean=sum(budgets)/2
    assert mean==s.Rational(250101,40000)
    assert 4*mean<100
    return {'C_inverse_ones':[6,-30,30],'permutation_vectors':[list(map(int,a)) for a in permutation_vectors],'ratio_QC2_over_Q':[str(s.Rational(67,75)),str(s.Rational(19,15))],'skew_quadratic_values':['1151/900','-11/20'],'phi_L':-10,'functional_mass':str(v.dot(v)),'full_product_readouts':readouts,'wedge_order':['12','13','14','23','24','34'],'wedges':wedges,'budgets':list(map(str,budgets)),'budget_mean':str(mean),'root_sum_upper_bound':'sqrt(250101)/100 < 10'}


def universal_identity_checks():
    """Free associative algebra checks, with no selected matrix realization."""
    import sympy as s
    c,x,t=s.symbols('c x t',commutative=False)
    delta=s.symbols('delta',real=True)
    jordan=lambda a,b:(a*b+b*a)/2
    checks=[]
    for sign in [1,-1]:
        D=lambda a,b:sign*s.I*(a*b-b*a)/2
        h=D(c,x); g=D(c,t); r=jordan(c,x); k=jordan(c,t)
        R=x+delta*(h-k); S=t+delta*(r+g); Z=R+s.I*S
        expected=(1+delta*s.I*c)*(x+s.I*t) if sign==1 else (x+s.I*t)*(1+delta*s.I*c)
        assert s.expand(Z-expected)==0
        # Formal adjoint reverses order and conjugates i, with c,x,t fixed.
        Zstar=R-s.I*S
        product=Zstar*Z if sign==1 else Z*Zstar
        expected_product=(x-s.I*t)*(1+c*c)*(x+s.I*t) if sign==1 else (x+s.I*t)*(1+c*c)*(x-s.I*t)
        for d in [1,-1]:
            assert s.expand(product.subs(delta,d)-expected_product)==0
        L=-2*jordan(x,h-k)-2*jordan(t,r+g)
        comm=lambda a,b:a*b-b*a
        commdiff=comm(R.subs(delta,1),S.subs(delta,1))-comm(R.subs(delta,-1),S.subs(delta,-1))
        assert s.expand(L-sign*s.I*commdiff/2)==0
        # The unstarred display as a candidate is false as a free identity.
        assert s.expand(Z.subs(delta,1)**2-Z.subs(delta,-1)**2)!=0
        checks.append({'D_sign':sign,'factorization':True,'correct_adjoint_product_independent_of_delta':True,'L_commutator_identity':True,'unstarred_product_independence':False})
    # The ell(X)c term vanishes on 1-perp, identically in both coordinates.
    C=s.Matrix(3,3,lambda i,j:s.Rational(2,i+j+2)); Lam=s.diag(1,2,3)
    a,b=s.symbols('a b'); y=s.Matrix([a,b,-a-b])
    assert s.expand((y.T*C*Lam*y)[0])==0
    # Independent sign averaging on all three skew coordinates.
    from itertools import product
    m12,m13,m23=s.symbols('m12 m13 m23')
    M=s.Matrix([[0,m12,m13],[-m12,0,m23],[-m13,-m23,0]])
    avg=s.zeros(3); weighted=s.zeros(3)
    for signs in product([-1,1],repeat=3):
        T=s.diag(*signs); N=T*M*T
        avg+=N/8; weighted+=signs[0]*signs[1]*N/8
    assert avg==s.zeros(3)
    assert weighted==m12*s.Matrix([[0,1,0],[-1,0,0],[0,0,0]])
    # Spectral permutations span diagonal coordinates, without covariance.
    from itertools import permutations
    assert s.Matrix(list(permutations([1,2,3]))).rank()==3
    return {'complex_optional_free_algebra':checks,'ell_X_term_zero_on_kernel':True,'unweighted_skew_sign_average_zero':True,'weighted_average_keeps_only_M12':True,'spectral_permutations_span_all_diagonal_inputs':True}


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('mode',choices=['build','verify']); args=ap.parse_args()
    if args.mode=='build': build()
    else:
        result={'linear_certificate':verify(),'fixed_corner_exact_checks':fixed_corner_checks(),'universal_identity_checks':universal_identity_checks(),'status':'PASS','arithmetic':'rational, prime-field and free-algebra exact; no floating point','not_lean_kernel_certified':True}
        (ROOT/'exact_results.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))
