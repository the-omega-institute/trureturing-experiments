"""Top-only rank-one vacuum using just binary relative signs on shortest scalar edges."""
from itertools import combinations, permutations
from sympy import Matrix, diag, simplify, sqrt, zeros, eye, Rational as Q


def shell_data():
    pts = list(permutations((0,1,3,4)))
    a,b = pts[:2]
    alpha = (sqrt(5)-1)/2
    labels = dict(zip(a,(1,alpha,-alpha,-1)))
    verts,chi,signs = [],[],[]
    for p in pts:
        z = [labels[i] for i in p]
        u = Matrix([z[0]+z[1],z[0]+z[2],z[0]+z[3]])
        verts.append((u*u.T/(8*(2+2*alpha**2))).applyfunc(simplify))
        chi.append(1 if sum(x*(i-j) for x,i,j in zip(p,a,b))>0 else -1)
        signs.append(-1 if p[0]+p[1]==4 else 1)
    return pts,verts,chi,signs


def main():
    pts,V,chi,s = shell_data()
    total = sum(V,zeros(3)).applyfunc(simplify)
    assert total==eye(3)
    upper = sum((z*h for z,h in zip(s,V)),zeros(3)).applyfunc(simplify)
    lower = sum((z*b*h for z,b,h in zip(s,chi,V)),zeros(3)).applyfunc(simplify)
    assert upper==diag(1,0,0) and lower==zeros(3)
    counts = {d:0 for d in (1,2,4)}
    for i,j in combinations(range(24),2):
        d = 26-sum(x*y for x,y in zip(pts[i],pts[j]))
        if d in counts and s[i]!=s[j]:
            counts[d] += 1
    assert counts == {1:16,2:0,4:0}
    print("Unit scalar-site signs produce one up mass and zero down masses; only 16 shortest edges change.")


if __name__=="__main__":
    main()
