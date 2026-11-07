#!/usr/bin/env python3
"""Exact finite-grid checks for arXiv:2610.03438v1, Definition 3.

No floating point. Sector sums come from the literal ordered distinct-site
operator, not from the proposed off-diagonal formulas. Active-face enumeration
is performed in the original Bell coefficient space modulo its computed kernel.
"""
from itertools import combinations
from math import comb, factorial, gcd, lcm
from functools import reduce
from random import Random
import sympy as S

X, Z, u, v = S.symbols('X Z u v')
INDEX = [(j, s) for j in range(1, 4) for s in range(j + 1)]
R = S.Rational


def unit_circle(t):
    return 2*t/(1+t*t), (1-t*t)/(1+t*t)  # sin, cos


def primitive(xs):
    xs = [R(x) for x in xs]
    den = lcm(*(x.q for x in xs))
    ints = [int(x*den) for x in xs]
    g = reduce(gcd, (abs(x) for x in ints), 0)
    return tuple(x//g for x in ints) if g else tuple(ints)


def block(n, settings, j, s):
    """<sum wt k | S_(0^(j-s)1^s) | sum wt l>, unnormalized.

    A selected site has polynomial M00+M01*v+M10*u+M11*u*v.
    Unselected sites contribute 1+u*v. Permutation symmetry makes every
    ordered distinct-site placement contribute the same sector sum.
    """
    local = [c + a*(u+v) - c*u*v for a, c in settings]
    poly = S.Poly(S.expand(local[0]**(j-s)*local[1]**s*(1+u*v)**(n-j)), u, v)
    factor = factorial(n)//factorial(n-j)
    return S.Matrix(n+1, n+1, lambda k, l: factor*poly.coeff_monomial(u**k*v**l))


def bell_blocks(n, settings, order=3):
    return [block(n, settings, j, s) for j in range(1, order+1) for s in range(j+1)]


def offdiag_matrix(blocks):
    n = blocks[0].rows - 1
    return S.Matrix([[b[k, l] for b in blocks] for k in range(n+1) for l in range(k+1, n+1)])


def change(settings):
    cols = []
    for j, s in INDEX:
        pp = S.Poly(S.expand((settings[0][0]*X+settings[0][1]*Z)**(j-s)
                            *(settings[1][0]*X+settings[1][1]*Z)**s), X, Z)
        cols.append(S.Matrix([pp.coeff_monomial(X**t*Z**(j-t)) if h == j else 0
                              for h, t in INDEX]))
    return S.Matrix.hstack(*cols)


def pauli(z1=0, a=0, z2=0, b=0, e=0, z3=0, c=0, f=0, d=0):
    return S.Matrix([z1, a, z2, b, e, z3, c, f, d])


def generators(n):
    m = n-1
    return ([pauli(a=4*i*(m-i-1)-m*(m-1), b=2*(m-2*i-1), c=-1) for i in range(m)]
            + [pauli(a=-m*(m-1), c=1), pauli(e=-(m-1), f=1),
               pauli(e=-(m-1), f=-1), pauli(a=3*m*(m-1), c=-3, d=-4)])


def enumerate_rays(F):
    """Enumerate 1-dimensional faces after removing lineality by pivot columns.

    Every quotient extreme ray is an intersection of dim-1 independent
    active rows. All such subsets are examined; feasibility is exact.
    Return positive projective classes of their full off-diagonal images.
    """
    pivots = F.rref()[1]
    B = F[:, list(pivots)]
    rows = sorted({primitive(row) for row in B.tolist() if any(row)})
    rays = set()
    for chosen in combinations(rows, len(pivots)-1):
        null = S.Matrix(chosen).nullspace()
        if len(null) != 1:
            continue
        r = null[0]
        ev = [sum(a*b for a,b in zip(row,r)) for row in rows]
        if all(x <= 0 for x in ev):
            pass
        elif all(x >= 0 for x in ev):
            r = -r
        else:
            continue
        rays.add(primitive(B*r))
    return rays, len(pivots), F.cols-len(pivots), B


def fidelity(n, settings):
    blocks = bell_blocks(n, settings, 2)
    s0,c0 = settings[0]; s1,c1 = settings[1]; m=n-1
    for k in range(n):
        expected = [s0,s1,2*s0*c0*(m-2*k),(s0*c1+c0*s1)*(m-2*k),2*s1*c1*(m-2*k)]
        den = comb(n,k)*(n-k)
        assert [bb[k,k+1]/den for bb in blocks] == expected
    for k in range(n-1):
        den = comb(n,k)*(n-k)*(n-k-1)
        assert [bb[k,k+2]/den for bb in blocks] == [0,0,s0*s0,s0*s1,s1*s1]
    F=offdiag_matrix(blocks)
    rays,rank,nullity,_ = enumerate_rays(F)
    assert (len(rays),rank,nullity)==(3,3,2)


def generic(n, settings):
    blocks=bell_blocks(n,settings)
    F=offdiag_matrix(blocks)
    P=change(settings)
    assert P.det()!=0
    inv=P.inv()
    # Verify the inverse substitution in the original nine coefficient coordinates.
    s0,c0=settings[0];s1,c1=settings[1];delta=s1*c0-c1*s0
    # input Z -> (sin(theta)*U-sin(phi)*V)/delta,
    # input X -> (-cos(theta)*U+cos(phi)*V)/delta.
    U,V=S.symbols('U V')
    cols=[]
    for j,t in INDEX:
        pp=S.Poly(S.expand(((s1*U-s0*V)/delta)**(j-t)
                           *((-c1*U+c0*V)/delta)**t),U,V)
        cols.append(S.Matrix([pp.coeff_monomial(U**(j-s)*V**s) if h==j else 0 for h,s in INDEX]))
    assert S.Matrix.hstack(*cols)==inv
    # Independently computed physical blocks must have the proposed h_l.
    pauli_blocks=[sum((inv[q,p]*blocks[q] for q in range(9)),S.zeros(n+1)) for p in range(9)]
    m=n-1
    for k in range(n):
        den=comb(n,k)*(n-k)
        assert [bb[k,k+1]/den for bb in pauli_blocks] == [0,1,0,m-2*k,0,0,(m-2*k)**2-m,0,3*k*(m-k)]
    for k in range(n-1):
        den=comb(n,k)*(n-k)*(n-k-1)
        assert [bb[k,k+2]/den for bb in pauli_blocks] == [0,0,0,0,1,0,0,m-1-2*k,0]
    for k in range(n-2):
        den=comb(n,k)*(n-k)*(n-k-1)*(n-k-2)
        assert [bb[k,k+3]/den for bb in pauli_blocks] == [0,0,0,0,0,0,0,0,1]
    for k in range(n+1):
        for l in range(k+4,n+1):assert all(bb[k,l]==0 for bb in blocks)
    expected={primitive(F*inv*g) for g in generators(n)}
    rays,rank,nullity,_=enumerate_rays(F)
    assert rays==expected, (n,len(rays),len(expected))
    assert (rank,nullity,len(rays))==(6,3,n+3)
    D=S.Matrix.hstack(*[inv*S.eye(9)[:,q] for q in [0,2,5]])
    assert D.rank()==3 and F*D==S.zeros(F.rows,3)
    print(f'GENERIC n={n} angles_sin_cos={settings} rank={rank} lineality={nullity} rays={len(rays)} expected={n+3} exact_ray_set_match=True',flush=True)


def collinear(n, sc, eps):
    s,c=sc
    settings=[sc,(eps*s,eps*c)]
    F=offdiag_matrix(bell_blocks(n,settings))
    if s==0:
        assert F==S.zeros(F.rows,9)
        print(f'DIAGONAL n={n} epsilon={eps} cone=R^9 lineality=9',flush=True)
        return
    m=n-1; b=abs(c/s)
    xyz=[(1,0,0),(2*b*m,1,0),(3*b*b*m*(m-1),3*b*(m-1),1)]
    gs=[pauli(z1=-x/s,z2=-y/s**2,z3=-z/s**3) for x,y,z in xyz]
    rays,rank,nullity,_=enumerate_rays(F)
    assert rays=={primitive(F*g) for g in gs}
    assert (rank,nullity,len(rays))==(3,6,3)
    # Check the stated beta-kernel equals the actual off-diagonal kernel.
    beta=S.Matrix([[eps**q if j==h else 0 for j,q in INDEX] for h in range(1,4)])
    assert beta.rank()==3
    assert all(F*w==S.zeros(F.rows,1) for w in beta.nullspace())
    print(f'COLLINEAR n={n} sin_cos={sc} epsilon={eps} rank=3 lineality=6 rays=3 exact_ray_set_match=True',flush=True)


def symbolic():
    k,m,A2,A1,A0,w,u0,v0=S.symbols('k m A2 A1 A0 w u v')
    d=-w;c=-(A2+3*w)/4;b=(A1+m*A2)/2
    a=-A0-m*A1/2-m*(m+1)*A2/4+3*m*(m-1)*w/4
    h=a+b*(m-2*k)+c*((m-2*k)**2-m)+3*d*k*(m-k)
    assert S.expand(h+A2*k*k+A1*k+A0)==0
    e=-(u0+v0)/2;f=(v0-u0)/(2*(m-1))
    assert S.cancel(e+f*(m-1)+u0)==0
    assert S.cancel(e-f*(m-1)+v0)==0
    x,y,z,bb=S.symbols('x y z b')
    p=x+2*bb*y*(m-2*k)+3*z*(bb**2*((m-2*k)**2-m)+k*(m-k))
    diff=S.expand(p-p.subs(k,m))
    assert S.expand(diff-(m-k)*(4*bb*y+3*z*k*(1-4*bb**2)))==0
    assert S.expand(diff.subs(y,3*bb*(m-1)*z)-3*z*(m-k)*(4*bb**2*(m-1-k)+k))==0
    for mm in range(2,13):
        Q=-S.Matrix([[1,j,j*j] for j in range(mm+1)])
        rays,rank,nullity,_=enumerate_rays(Q)
        gs=[S.Matrix([i*(i+1),-2*i-1,1]) for i in range(mm)]+[S.Matrix([0,mm,-1])]
        assert rays=={primitive(Q*g) for g in gs}
        assert len(rays)==mm+1
    print('SYMBOLIC inverse_coordinates=True collinear_difference=True quadratic_cones_m=2..12 complete_ray_sets=True',flush=True)


def main():
    print('Python exact arithmetic; SymPy',S.__version__,flush=True)
    rng=Random(261003438)
    pairs=[]
    while len(pairs)<2:
        t0=R(rng.randint(1,13),rng.randint(2,17));t1=R(rng.randint(1,13),rng.randint(2,17))
        if t0!=t1:
            pairs.append((unit_circle(t0),unit_circle(t1)))
            print(f'RANDOM_ANGLE_PAIR half_angle_tangents=({t0},{t1}) sin_cos={pairs[-1]}',flush=True)
    symbolic()
    for n in range(3,7):
        for settings in pairs:
            fidelity(n,settings)
            generic(n,settings)
    print('FIDELITY source_Section_VI_two_body_hyperplanes=True n=3..6 both_angle_pairs rays=3 lineality=2',flush=True)
    for n in range(3,7):
        for sc in [unit_circle(R(1,3)),unit_circle(R(2)),unit_circle(R(-4,7)),unit_circle(R(1))]:
            for eps in [1,-1]:collinear(n,sc,eps)
        for eps in [1,-1]:collinear(n,(R(0),R(1)),eps)
    print('ALL_CHECKS_OK',flush=True)

if __name__=='__main__':main()
