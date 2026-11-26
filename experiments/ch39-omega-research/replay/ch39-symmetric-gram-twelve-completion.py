"""An explicit coordinate-permutation invariant rank-12 completion at cutoff 56.

Trivial block of rank 3, sign block of rank 1 and standard block of rank 4 give
rank W = 3 + 1 + 2*4 = 12.  All parameters are rational except q1 = sqrt(q1^2) and
one root t of a quadratic with rational discriminant, so the completion is a real
algebraic matrix.  Its positivity and rank follow from exact inertia counts; a
192-bit ball enclosure of W bounds the row-sum norm for the contraction criterion.
With 19 + 12 auxiliary states the three clocks are realised: n_aux <= 31.
"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

low = run_path(str(Path(__file__).with_name("ch39-symmetric-gram-lower-bound.py")))
K1, K2 = low["K1"], low["K2"]
ALPHA, S1 = fmpq(3, 50000000), fmpq(1, 4000000)


def largest(values):
    values = list(values)
    best = values[0]
    for value in values[1:]:
        best = best.max(value) if isinstance(best, arb) else max(best, value)
    return best


def row_norm(matrix):
    return largest(sum(abs(a) for a in row) for row in matrix.tolist())


def zcol(z):
    return fmpq_mat([[c] for c in z])


def ball(q):
    return arb(q.p)/arb(q.q)


def verify():
    ctx.prec = 192
    data = low["setup"]()
    dist = data["dist"]
    far = [(i, j) for i in range(19) for j in range(i+1, 19) if dist[i][j] > 56]
    assert len(far) == 153
    fs = data["fs"]
    inverse_norm = row_norm(data["g0"].inv())
    delta, sym, keys = data["delta"], data["sym"], data["keys"]
    gm = low["gamma"](data, delta)
    t0, s0, _ = low["blocks"](data, delta)
    r = [t0[0, j] for j in (1, 2, 3)]
    e, h = t0[0, 4], [t0[j, 4] for j in (1, 2, 3)]
    sigma = s0[0, 1]
    g = [h[j]-e/ALPHA*r[j] for j in range(3)]
    s = [S1, sigma/S1]
    u = [ALPHA]+r
    def gpp(rho):
        m = fmpq_mat(gm)
        gh = [0]+g
        for i in range(4):
            for j in range(4):
                m[i, j] += K1*(u[i]*u[j]/ALPHA+rho*gh[i]*gh[j])
        for i in range(2):
            for j in range(2):
                m[4+i, 4+j] += K2*s[i]*s[j]
        return m
    d0, d1, d2 = gpp(0).det(), gpp(1).det(), gpp(2).det()
    assert d2-d1 == d1-d0 != 0                          # det is affine in rho
    rho = -d0/(d1-d0)
    assert 2968*10**4 < rho < 2969*10**4
    gp = gpp(rho)
    assert gp.det() == 0 and low["inertia"](gp) == (4, 1, 1)
    # Tie: the trivial and sign diagonals at orbit B share one correction x2.
    kb = keys.index((4, 4))
    tri_b, sgn_b = [low["blocks"](data, sym[kb])[t][1, 1] for t in (0, 1)]
    for k, m in enumerate(sym):
        if k != kb:
            tt, ss, _ = low["blocks"](data, m)
            assert tt[1, 1] == 0 and ss[1, 1] == 0
    x2 = (s0[1, 1]-s[1]**2)/sgn_b
    tie = t0[1, 1]-tri_b*x2
    q1sq = tie-r[0]**2/ALPHA-rho*g[0]**2
    assert fmpq(2434, 10**9) < q1sq < fmpq(2435, 10**9)
    num, _ = gp.numer_denom()
    ker, dim = num.nullspace()
    assert dim == 1
    z = [fmpq(ker[i, 0]) for i in range(6)]
    assert z[2] != 0
    # A rational basis of z-perp: inertia (4, 1, 0) of Gamma'' restricted there.
    piv = max(range(6), key=lambda i: z[i] != 0)
    perp = [[(fmpq(1) if i == j else 0) - (z[j]/z[piv] if i == piv else 0)
             for j in range(6) if j != piv] for i in range(6)]
    nb = fmpq_mat(perp)
    assert nb.rank() == 5 and (zcol(z).transpose()*nb) == fmpq_mat(1, 5)
    assert low["inertia"](nb.transpose()*gp*nb) == (4, 1, 0)
    zz = fmpq_mat([[z[i]*z[j] for j in range(6)] for i in range(6)])
    x = (gp+zz).inv()
    ph = fmpq_mat([[0], [1], [-z[1]/z[2]], [0], [0], [0]])
    dh = fmpq_mat([[0], [0], [z[3]], [-z[2]], [0], [0]])
    zc = fmpq_mat([[c] for c in z])
    assert (zc.transpose()*ph)[0, 0] == 0 and (zc.transpose()*dh)[0, 0] == 0
    form = lambda a, b: (a.transpose()*x*b)[0, 0]
    # q = q1*ph + t*dh;  1 + k1 q^T (Gamma'')^+ q = 0 is a2 t^2 + q1 b t + c0 = 0.
    a2, b, c0 = K1*form(dh, dh), 2*K1*form(ph, dh), 1+K1*q1sq*form(ph, ph)
    disc = q1sq*b**2-4*a2*c0
    assert a2 > 0 and disc > 0                         # root=+1 is the larger root
    def complete(root):
        q1 = ball(q1sq).sqrt()
        t = (-q1*ball(b)+root*ball(disc).sqrt())/(2*ball(a2))
        q = [q1*ball(ph[i, 0])+t*ball(dh[i, 0]) for i in (1, 2, 3)]
        # Trivial block: Q = [[alpha, e], [e, alpha_E]], middle = K^T Q^{-1} K + q q^T.
        alpha_e = e**2/ALPHA+1/rho
        qm = fmpq_mat([[ALPHA, e], [e, alpha_e]])
        km = fmpq_mat([r, h])
        mid = km.transpose()*qm.inv()*km
        assert mid == fmpq_mat([[r[i]*r[j]/ALPHA+rho*g[i]*g[j] for j in range(3)] for i in range(3)])
        tm = arb_mat(5, 5)
        order = [0, 4]
        for i in range(2):
            for j in range(2):
                tm[order[i], order[j]] = ball(qm[i, j])
        for i in range(3):
            tm[0, 1+i] = tm[1+i, 0] = ball(r[i])
            tm[4, 1+i] = tm[1+i, 4] = ball(h[i])
            for j in range(3):
                tm[1+i, 1+j] = ball(mid[i, j])+q[i]*q[j]
        assert alpha_e*ALPHA-e**2 > 0                        # Q is positive definite
        # Solve for the nine invariant corrections from nine entries of the two blocks.
        entries = [(0, 0, 0), (1, 0, 0), (0, 1, 1), (0, 1, 2), (0, 1, 3), (0, 2, 2), (0, 2, 3),
                   (0, 3, 3), (0, 4, 4)]
        shifts = [low["blocks"](data, m) for m in sym]
        coef = fmpq_mat([[shifts[k][b_][i, j] for k in range(9)] for b_, i, j in entries])
        assert coef.det() != 0
        sm = [[ball(s[0]**2), ball(sigma)], [ball(sigma), ball(s[1]**2)]]
        target = [tm[i, j] if b_ == 0 else sm[i][j] for b_, i, j in entries]
        base = [t0[i, j] if b_ == 0 else s0[i, j] for b_, i, j in entries]
        inv = coef.inv()
        xs = [sum((ball(inv[k, l])*(ball(base[l])-target[l]) for l in range(9)), arb(0))
              for k in range(9)]
        w = arb_mat([[ball(delta[i, j]) for j in range(19)] for i in range(19)])
        for k in range(9):
            for i in range(19):
                for j in range(19):
                    if sym[k][i, j] != 0:
                        w[i, j] -= xs[k]
        # Consistency: the enclosure reproduces the sign B diagonal (the tie) and the
        # trivial and sign blocks.
        fg, ft = data["fg"], data["ft"]
        fga = arb_mat([[ball(fg[i, j]) for j in range(2)] for i in range(19)])
        fta = arb_mat([[ball(ft[i, j]) for j in range(5)] for i in range(19)])
        sgn_w, tri_w = fga.transpose()*w*fga, fta.transpose()*w*fta
        assert all((sgn_w[i, j]-sm[i][j]).contains(0) for i in range(2) for j in range(2))
        assert all((tri_w[i, j]-tm[i, j]).contains(0) for i in range(5) for j in range(5))
        norm = row_norm(w)
        # Consistency of the standard block: the enclosure of Lambda = FS^T W FS contains
        # Gamma'' + (2/9) qbar qbar^T and annihilates z up to the enclosure radius.  Its
        # inertia (4, 0, 2) follows from the written Haynsworth argument; the balls only
        # check that the matrix built from W is the one that argument is about.
        fsa = arb_mat([[ball(fs[i, j]) for j in range(6)] for i in range(19)])
        lam = fsa.transpose()*w*fsa
        qbar = [arb(0)]+q+[arb(0), arb(0)]
        assert all((lam[i, j]-ball(gp[i, j])-ball(K1)*qbar[i]*qbar[j]).contains(0)
                   for i in range(6) for j in range(6))
        lz = lam*arb_mat([[ball(c)] for c in z])
        assert all(lz[i, 0].contains(0) for i in range(6))
        return norm
    bound = (1-arb(fmpq(98, 1000)))**2
    assert inverse_norm < 315000 and bound > arb("0.8136")
    norm = complete(1)            # q1 > 0 and the larger root of the quadratic
    other = complete(-1)          # q1 > 0 and the smaller root
    assert arb("8.88e-7") < norm < arb("8.90e-7")
    assert arb("0.2797") < 315000*norm < arb("0.2799")
    assert arb("8.64e-7") < other < arb("8.65e-7")
    assert arb("0.2721") < 315000*other < arb("0.2725")
    print("153 far entries; rho", float(rho), "q1^2", float(q1sq), "discriminant", ball(disc).str(6))
    print("block ranks (3, 1, 4): Gram rank 12; row-sum norms", norm.str(8), other.str(8),
          "contraction products", (315000*norm).str(6), (315000*other).str(6))
    print("Cutoff 56: 19 + 12 =", 19+12, "auxiliary states suffice; interval 29 <= n_aux <= 31")


if __name__ == "__main__":
    verify()
