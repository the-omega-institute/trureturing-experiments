from math import comb


def add(p, q):
    n = max(len(p), len(q))
    out = [0] * n
    for i, a in enumerate(p): out[i] += a
    for i, a in enumerate(q): out[i] += a
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out


def neg(p): return [-a for a in p]


def sub(p, q): return add(p, neg(q))


def mul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): out[i+j] += a * b
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out


def det3(M):
    a,b,c = M[0]
    d,e,f = M[1]
    g,h,i = M[2]
    return add(sub(mul(a, sub(mul(e, i), mul(f, h))),
                  mul(b, sub(mul(d, i), mul(f, g)))),
               mul(c, sub(mul(d, h), mul(e, g))))


def subset_table(r, nmax):
    S = [[0] * (nmax + 1) for _ in range(nmax + 1)]
    S[0][0] = 1
    for n in range(1, nmax + 1):
        for k in range(1, n + 1):
            S[n][k] = comb(n + (r-1)*k - 1, r-1) * S[n-1][k-1] + k * S[n-1][k]
    return S


def row_poly(S, n):
    return S[n][:n+1]


def minor(r, m):
    S = subset_table(r, m + 4)
    M = [[row_poly(S, m+i+j) for j in range(3)] for i in range(3)]
    return det3(M)


def second_order_data(r, nmax):
    S = subset_table(r, nmax)
    a = [S[n][2] if n >= 2 else 0 for n in range(nmax+1)]
    b = [a[n+2] - 2*a[n+1] + a[n] for n in range(nmax-2)]
    d = [comb(n + 2*r - 2, r-3) for n in range(nmax-1)]
    return a,b,d


def main():
    # Fidelity with Appendix tables: S^(4) row n=2 is (0,1,35).
    S4 = subset_table(4, 2)
    assert S4[2][:3] == [0, 1, 35]
    # The source's Appendix C witness is the m=1, r=3 minor.
    witness = minor(3, 1)
    assert witness[:5] == [0, 0, 0, 0, 0]
    assert witness[5] < 0
    print('fidelity r=4,n=2 coefficients:', S4[2][:3])
    print('source witness r=3,m=1 determinant coefficients:', witness)

    # Verify the claimed all-r, all-m obstruction and the finite-difference identities.
    for r in range(3, 16):
        for m in range(1, 16):
            D = minor(r, m)
            assert all(c == 0 for c in D[:5]), (r, m, D[:5])
            assert D[5] < 0, (r, m, D[5])
            a,b,d = second_order_data(r, m + 5)
            bm = b[m]
            dm = d[m]
            assert D[5] == bm * (d[m+1] - 2*dm) - dm*dm
            for n in range(1, m + 2):
                assert b[n] == a[n] + comb(n + 2*r - 1, r-1)
                assert b[n+1] == 2*b[n] + d[n]
            if r == 3:
                assert d[m] == 1
            else:
                assert d[m+1] - 2*dm < 0
    print('checked exact determinant coefficients for r=3..15, m=1..15')
    print('checked [x^5] identity and strict negativity')
    print('ALL_CHECKS_OK')

if __name__ == '__main__':
    main()
