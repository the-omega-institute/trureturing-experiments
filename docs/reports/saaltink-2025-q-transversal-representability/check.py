"""Check of a representation of q-transversal matroids (arXiv:2503.12201v1, Conjecture 3) over q = 2.

For a family X_1, ..., X_k of subspaces of F_2^n, A_i = ann(X_i) with basis l_i1, ..., l_ia_i, B_i = prod_{h<i}(a_h+1),
D = prod_i (a_i+1) and alpha a root of an irreducible polynomial of degree D over F_2, the matrix G with rows
g_i = sum_j alpha^(j B_i) l_ij over F_(2^D) is compared with the q-transversal matroid: for every subspace U of F_2^n,
rank(G S_U) (S_U a basis of U as columns) must equal the largest dimension of a partial q-transversal inside U.
Independence is computed both from Definition 14 (every basis admits an avoiding injection) and from Theorem 18.
The families are random subspaces (not coordinate subspaces) for (n, k) = (3, 3) and (4, 2).
"""
import itertools, random, sys
# GF(2)^n vectors as ints (bitmasks). Subspaces as frozenset of vectors.
def span(vs):
    S={0}
    for v in vs:
        S|={s^v for s in S}
    return frozenset(S)
def rank2(vs):
    basis=[]
    for v in vs:
        for b in basis:
            v=min(v,v^b)
        if v: basis.append(v)
    return len(basis)
def dim(S): return len(S).bit_length()-1
def all_subspaces(n):
    seen=set(); out=[]
    vecs=list(range(1,2**n))
    for t in range(n+1):
        for combo in itertools.combinations(vecs,t):
            if rank2(list(combo))==t:
                S=span(combo)
                if S not in seen: seen.add(S); out.append(S)
    return out
def basis_of(S):
    b=[]
    for v in sorted(S):
        if v and rank2(b+[v])>len(b): b.append(v)
    return b
def all_bases(S):
    t=dim(S); vs=[v for v in S if v]
    return [c for c in itertools.combinations(vs,t) if rank2(list(c))==t]
def has_injection(B,X):
    # bipartite matching: b -> i with b not in X[i]
    k=len(X); match={}
    def try_(b,seen):
        for i in range(k):
            if b not in X[i] and i not in seen:
                seen.add(i)
                if i not in match or try_(match[i],seen): match[i]=b; return True
        return False
    return all(try_(b,set()) for b in B)
def indep_def(U,X): return all(has_injection(B,X) for B in all_bases(U))
def indep_thm18(U,X):
    k=len(X)
    for r in range(k+1):
        for J in itertools.combinations(range(k),r):
            XJ=frozenset(range(0,2**N)) if not J else frozenset.intersection(*[X[i] for i in J])
            if dim(U & XJ)+len(J)>k: return False
    return True
# GF(2^D)
def pmod(a,m):
    dm=m.bit_length()-1
    while a and a.bit_length()-1>=dm: a^=m<<(a.bit_length()-1-dm)
    return a
def pmul(a,b,m):
    r=0
    while b:
        if b&1: r^=a
        b>>=1; a<<=1
        if a>>(m.bit_length()-1)&1: a^=m
    return r
def ppow(a,e,m):
    r=1
    while e:
        if e&1: r=pmul(r,a,m)
        a=pmul(a,a,m); e>>=1
    return r
def pgcd(a,b):
    while b: a,b=b,pmod(a,b)
    return a
def irreducible(D):
    # Rabin: x^(2^D) = x mod f and gcd(x^(2^(D/p)) - x, f)=1 for primes p|D
    ps=[p for p in range(2,D+1) if D%p==0 and all(p%r for r in range(2,p))]
    while True:
        f=(1<<D)|random.getrandbits(D)|1
        x=2; y=x
        for _ in range(D): y=pmul(y,y,f)
        if y!=pmod(x,f): continue
        ok=True
        for p in ps:
            z=x
            for _ in range(D//p): z=pmul(z,z,f)
            if pgcd(z^x,f)!=1: ok=False; break
        if ok: return f
def pinv(a,m): return ppow(a,(1<<(m.bit_length()-1))-2,m)
def rank_gf(M,m):
    M=[row[:] for row in M]; r=0; cols=len(M[0]) if M else 0
    for c in range(cols):
        piv=next((i for i in range(r,len(M)) if M[i][c]),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; iv=pinv(M[r][c],m)
        M[r]=[pmul(v,iv,m) for v in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[M[i][j]^pmul(f,M[r][j],m) for j in range(cols)]
        r+=1
    return r
def annihilator_basis(X,n):
    # functionals l in F2^n (as bitmask) with l.x = 0 for all x in X
    ann=[l for l in range(2**n) if all(bin(l&x).count('1')%2==0 for x in X)]
    return basis_of(frozenset(ann))

def run(seed, n, k, trials):
    global N
    random.seed(seed); N = n
    subs = all_subspaces(n)
    allok = True
    for tr in range(trials):
        X = [random.choice(subs) for _ in range(k)]
        A = [annihilator_basis(Xi, n) for Xi in X]; a = [len(Ai) for Ai in A]
        B = [1] * k
        for i in range(1, k):
            B[i] = B[i - 1] * (a[i - 1] + 1)
        D = B[-1] * (a[-1] + 1)
        m = irreducible(D); alpha = 2
        G = []
        for i in range(k):
            row = [0] * n
            for j, l in enumerate(A[i], start=1):
                c = ppow(alpha, j * B[i], m)
                for t in range(n):
                    if l >> t & 1:
                        row[t] ^= c
            G.append(row)
        ok = True
        for U in subs:
            if indep_def(U, X) != indep_thm18(U, X):
                ok = False
            qr = max(dim(W) for W in subs if W <= U and indep_thm18(W, X))
            bU = basis_of(U)
            gr = 0
            if bU:
                GS = [[0] * len(bU) for _ in range(k)]
                for i in range(k):
                    for c, u in enumerate(bU):
                        s = 0
                        for t in range(n):
                            if u >> t & 1:
                                s ^= G[i][t]
                        GS[i][c] = s
                gr = rank_gf(GS, m)
            if gr != qr:
                ok = False
        print(f"seed {seed} n={n} k={k} trial {tr}: dims X={[dim(x) for x in X]} ann dims={a} D={D} "
              f"subspaces={len(subs)} Definition 14 = Theorem 18 and rank(G S_U) = q-rank for all U: {ok}")
        allok &= ok
    return allok


ok = run(1, 3, 3, 6) & run(2, 4, 2, 3)
print("ALL_CHECKS", ok)
raise SystemExit(0 if ok else 1)
