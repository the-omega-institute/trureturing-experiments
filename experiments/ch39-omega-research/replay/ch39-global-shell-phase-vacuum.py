"""全壳物质数真空的全局能量比较：余弦二次下界与整格饱和。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, arb_poly, fmpz_mat

phase = run_path(str(Path(__file__).with_name("ch39-matter-number-three-spokes.py")))
soft = phase["soft"]
pole = soft["pole"]
primary, mirror = soft["star"]["charged"]["orbits"][:2]


def harmonics(matrix, weight):
    omega = acb(-arb(1)/2, arb(3).sqrt()/2)
    roots = (acb(1),omega,omega.conjugate())
    constant = arb_poly()
    cosine, sine = ([arb_poly() for _ in range(3)] for _ in range(2))
    for a,b in product(roots, repeat=2):
        m = acb_mat(matrix)
        for group, unit, stiffness in ((primary,a,weight),(mirror,b,1)):
            for j in group:
                m[pole,pole] += stiffness-1
                m[j,j] += stiffness-1
                m[pole,j] = -stiffness*unit
                m[j,pole] = m[pole,j].conjugate()
        coeff = (-m).charpoly().coeffs()
        assert all(z.imag.contains(0) for z in coeff)
        p = arb_poly([z.real/9 for z in coeff])
        constant += p
        for i,z in enumerate((a*b.conjugate(), b, a.conjugate())):
            cosine[i] += p*(2*z.real)
            sine[i] += p*(2*z.imag)
    for group in (cosine,sine):
        for i,p in enumerate(group):
            assert all(p[k].contains(0) for k in (17,18,19))
            group[i] = arb_poly(p.coeffs()[:17])
    return constant, cosine, sine


def energy_bound(matrices, weight, c):
    x = (c,arb(1),c)
    Q = arb_mat(3,3)
    gy_minus_gz = arb(0)
    integrate = phase["integrate_force"]
    for q,n in ((0,33),(5,18),(10,3)):
        a,b,d = harmonics(matrices[q],weight)
        p = a+sum((b[i]*x[i] for i in range(3)),arb_poly())
        def integral(top, bottom):
            return n*integrate(top,bottom,"1e-7")/(8*arb.pi()**2)
        if q == 0:
            assert all(v < 0 for poly in b for v in poly.coeffs())
            maximum = a-sum(b,arb_poly())
            gy_minus_gz -= integral(b[1]-b[2],p)
            for i in range(3):
                for j in range(i,3):
                    Q[i,j] += integral(b[i]*b[j]*arb(".5"),maximum**2)
            continue
        # sin(theta_i) sin(theta_j) = cos(theta_i) cos(theta_j) - cos(theta_k).
        mat = [[b[i]*b[j]+(1 if i==j else -1)*d[i]*d[j] for j in range(3)] for i in range(3)]
        at_root = p**2-(1-c*c)*(d[0]-d[2])**2
        linear = [2*(a*b[i]+d[(i+1)%3]*d[(i+2)%3]) for i in range(3)]
        grad = [linear[i]+2*sum((m*t for m,t in zip(row,x)),arb_poly()) for i,row in enumerate(mat)]
        gy_minus_gz -= integral(grad[1]-grad[2],at_root)
        for i in range(3):
            for j in range(i,3):
                Q[i,j] -= integral(mat[i][j],at_root)
    for i in range(3):
        for j in range(i):
            Q[i,j] = Q[j,i]
    transform = arb_mat([[1,1,0],[0,0,1],[1,0,0]])
    K = transform.transpose()*Q*transform
    assert K[0,0] > arb(".006")
    reduced = [[K[i,j]-K[0,i]*K[0,j]/K[0,0] for j in (1,2)] for i in (1,2)]
    # Both error coordinates have square at most 2(1-cos beta).
    loss = 2*sum(abs(v) for row in reduced for v in row)
    assert -gy_minus_gz-loss > arb(".147")
    print("Global coefficient of 1-cos(beta), retaining a positive square",
          (-gy_minus_gz-loss).lower().str(10),K[0,0].lower().str(10))


def flat_torus():
    edges = [(i,j) for i,j,_,_ in soft["geo"]["edges"] if pole not in (i,j)]
    index = {pair:j for j,pair in enumerate(edges)}
    # A spanning tree fixes the gauge among the 18 outer sites.
    reached = {edges[0][0]}
    tree = set()
    while len(reached)<18:
        pair = next(pair for pair in edges if len(reached.intersection(pair))==1)
        reached.update(pair)
        tree.add(pair)
    extra = [pair for pair in edges if pair not in tree]
    triangles = []
    for i,j,k in combinations(sorted(reached),3):
        if all(pair in index for pair in ((i,j),(j,k),(i,k))):
            cycle = {(i,j):1,(j,k):1,(i,k):-1}
            triangles.append([cycle.get(pair,0) for pair in extra])
    normal = fmpz_mat(triangles).snf()
    assert len(extra)==100 and all(normal[i,i]==1 for i in range(100))
    print("One hundred integral Smith factors equal to one")


def verify():
    matrices,weight,c,_,_ = phase["verify"]()
    energy_bound(matrices,weight,c)
    flat_torus()


if __name__ == "__main__":
    verify()
