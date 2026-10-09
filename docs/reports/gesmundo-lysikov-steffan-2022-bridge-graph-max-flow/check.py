import json, numpy as np, sympy as sp
P = 1_000_003

def rank_mod(a):
    a = np.array(a, dtype=np.int64, copy=True) % P
    nr, nc = a.shape
    r = 0
    for c in range(nc):
        piv = np.flatnonzero(a[r:, c])
        if piv.size == 0:
            continue
        s = r + int(piv[0])
        if s != r:
            a[[r, s]] = a[[s, r]]
        inv = pow(int(a[r, c]), P - 2, P)
        a[r, c:] = (a[r, c:] * inv) % P
        if r + 1 < nr:
            f = a[r + 1:, c].copy()
            nz = f != 0
            if np.any(nz):
                a[r + 1:, c:] = (a[r + 1:, c:] - f[:, None] * a[r, c:]) % P
        r += 1
        if r == nr:
            break
    return r

def rectangular(rows, cols):
    z = np.zeros((rows, cols), dtype=np.int64)
    for i in range(min(rows, cols)):
        z[i, i] = 1
    return z

def shift_block(x):
    s1 = np.zeros((x, x + 1), dtype=np.int64)
    s2 = np.zeros_like(s1)
    for i in range(x):
        s1[i, i] = 1
        s2[i, i + 1] = 1
    return s1, s2

def transposed_shift_block(x):
    a, b = shift_block(x)
    return a.T.copy(), b.T.copy()

def cyclic_data(A, B, G, D):
    """Return X,Y,B,D over F_P for the normalized short-short orientation."""
    U = np.zeros((A, A), dtype=np.int64)
    V = np.zeros((G, G), dtype=np.int64)
    for i in range(A):
        U[(i - 1) % A, i] = 1
    for j in range(G - 1):
        V[j + 1, j] = 1
    V[0, G - 1] = (P + 1) // 2
    Bm = np.zeros((B, A), dtype=np.int64)
    Dm = np.zeros((G, D), dtype=np.int64)
    for k in range(B):
        Bm[k, (k * A) // B] = 1
    for h in range(D):
        Dm[(h * G) // D, h] = 1
    return (-U) % P, V, Bm, Dm

def base_slices(p, alpha, beta, gamma, delta):
    a = p * alpha + (p + 1) * beta
    b = (p + 1) * alpha + (p + 2) * beta
    c = p * gamma + (p + 1) * delta
    d = (p + 1) * gamma + (p + 2) * delta
    M1 = np.zeros((a, b), dtype=np.int64)
    M2 = np.zeros((a, b), dtype=np.int64)
    N1 = np.zeros((d, c), dtype=np.int64)
    N2 = np.zeros((d, c), dtype=np.int64)
    row_m = [i * p for i in range(alpha)] + [alpha * p + i * (p + 1) for i in range(beta)]
    col_m = [i * (p + 1) for i in range(alpha)] + [alpha * (p + 1) + i * (p + 2) for i in range(beta)]
    row_n = [i * (p + 1) for i in range(gamma)] + [gamma * (p + 1) + i * (p + 2) for i in range(delta)]
    col_n = [i * p for i in range(gamma)] + [gamma * p + i * (p + 1) for i in range(delta)]
    for i in range(alpha):
        r1, r2 = shift_block(p)
        M1[np.ix_(range(row_m[i], row_m[i] + p), range(col_m[i], col_m[i] + p + 1))] = r1
        M2[np.ix_(range(row_m[i], row_m[i] + p), range(col_m[i], col_m[i] + p + 1))] = r2
    for i in range(beta):
        r1, r2 = shift_block(p + 1)
        rr = range(row_m[alpha + i], row_m[alpha + i] + p + 1)
        cc = range(col_m[alpha + i], col_m[alpha + i] + p + 2)
        M1[np.ix_(rr, cc)] = r1
        M2[np.ix_(rr, cc)] = r2
    for j in range(gamma):
        r1, r2 = transposed_shift_block(p)
        rr = range(row_n[j], row_n[j] + p + 1)
        cc = range(col_n[j], col_n[j] + p)
        N1[np.ix_(rr, cc)] = r1
        N2[np.ix_(rr, cc)] = r2
    for j in range(delta):
        r1, r2 = transposed_shift_block(p + 1)
        rr = range(row_n[gamma + j], row_n[gamma + j] + p + 2)
        cc = range(col_n[gamma + j], col_n[gamma + j] + p + 1)
        N1[np.ix_(rr, cc)] = r1
        N2[np.ix_(rr, cc)] = r2
    M3 = np.zeros((a, b), dtype=np.int64)
    N3 = np.zeros((d, c), dtype=np.int64)
    def put_m(U, r, s):
        for i in range(U.shape[0]):
            for j in range(U.shape[1]):
                M3[row_m[alpha + i] + r, col_m[j] + s] = U[i, j]
    def put_n(V, r, s):
        for i in range(V.shape[0]):
            for j in range(V.shape[1]):
                N3[row_n[i] + r, col_n[gamma + j] + s] = V[i, j]
    if alpha >= beta and gamma >= delta:
        X, Y, Bm, Dm = cyclic_data(alpha, beta, gamma, delta)
        for i in range(alpha):
            for j in range(alpha):
                M3[row_m[i], col_m[j]] = X[i, j]
        for i in range(gamma):
            for j in range(gamma):
                N3[row_n[i], col_n[j]] = Y[i, j]
        for r in range(p + 1):
            put_m(Bm, r, p - r)
            put_n(Dm, p - r, r)
        reason = 'short-short cyclic'
    elif alpha <= beta and gamma <= delta:
        X, Y, Bm, Dm = cyclic_data(beta, alpha, delta, gamma)
        X, Y, Bm, Dm = X.T.copy(), Y.T.copy(), Bm.T.copy(), Dm.T.copy()
        for i in range(beta):
            for j in range(beta):
                M3[row_m[alpha + i], col_m[alpha + j]] = X[i, j]
        for i in range(delta):
            for j in range(delta):
                N3[row_n[gamma + i], col_n[gamma + j]] = Y[i, j]
        for r in range(p + 1):
            put_m(Bm, r, p - r)
            put_n(Dm, p - r, r)
        reason = 'long-long cyclic'
    elif (alpha <= beta and delta <= gamma) or (alpha >= beta and gamma <= delta):
        put_m(rectangular(beta, alpha), p, 0)
        put_n(rectangular(gamma, delta), 0, p)
        reason = 'single'
    else:
        raise AssertionError((alpha, beta, gamma, delta))
    return (M1, M2, N1, N2, M3, N3, (a, b, c, d), reason)

def flow_rank(slices, dims):
    M1, M2, N1, N2, M3, N3 = slices
    a, b, c, d = dims
    F = (np.kron(M1, N1) + np.kron(M2, N2) + np.kron(M3, N3)) % P
    return rank_mod(F)


A,B,G,D = 3,2,4,3
X,Y,Bm,Dm = cyclic_data(A,B,G,D)
def rational(m):
    return sp.Matrix([[sp.Rational(1,2) if int(v)==(P+1)//2 else -1 if int(v)==P-1 else int(v) for v in row] for row in m])
X,Y,Bm,Dm=map(rational,[X,Y,Bm,Dm])
core=sp.kronecker_product(Bm,sp.eye(G))*(sp.eye(A*G)+sp.kronecker_product(X,Y)).inv()*sp.kronecker_product(sp.eye(A),Dm)
assert core.rank()==7 and min(B*G,A*D)==8
pairs=[(a,b) for a in range(1,15) for b in range(a+1,3*a) if a*a+b*b<=3*a*b]
count=failures=0
for a,b in pairs:
    m=b-a;p=(a-1)//m;beta=a-p*m;alpha=m-beta
    if p!=0:continue
    for c,d in pairs:
        n=d-c;q=(c-1)//n;delta=c-q*n;gamma=n-delta
        if q!=0 or min(alpha*delta,beta*gamma)==0:continue
        slices=base_slices(0,alpha,beta,gamma,delta)
        rank=flow_rank(slices[:6],(a,b,c,d))
        count+=1;failures+=rank!=min(a*d,b*c)
assert (failures,count)==(261,3364),(failures,count)
print(json.dumps({"kappa0":{"multiplicities":[A,B,G,D],"field":"Q","rank":core.rank(),"target":8},"p0":{"max_a_c":14,"scope":"positive ordered cone pairs with both depths zero and positive width-two defect","field_prime":P,"pairs":count,"deficient":failures}}))

exact_slices=base_slices(0,3,5,4,7)
assert exact_slices[6]==(5,13,7,18)
qM1,qM2,qN1,qN2,qM3,qN3=map(rational,exact_slices[:6])
exact_p0=sp.kronecker_product(qM1,qN1)+sp.kronecker_product(qM2,qN2)+sp.kronecker_product(qM3,qN3)
assert exact_p0.rank()==88
print(json.dumps({"p0_exact":{"dimensions":[5,13,7,18],"field":"Q","rank":88,"target":90}}))
