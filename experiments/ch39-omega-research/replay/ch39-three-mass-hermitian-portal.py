"""Exact spectral and bubble checks for one Hermitian three-layer scalar sector."""

import sympy as sp


def check_invariants():
    a, b, c, x, s = sp.symbols("a b c x s", real=True)
    d = sp.symbols("d1:4", real=True)
    p, q, h, k = sp.symbols("p q h k", real=True)
    r = h+sp.I*k
    portal = sp.Matrix([[d[0], p, sp.conjugate(r)], [p, d[1], q], [r, q, d[2]]])
    masses = (a, b, c)
    linear = [x-m-s*f for m, f in zip(masses, d)]
    expected = (sp.prod(linear)-s**2*(p**2*linear[2]+(h*h+k*k)*linear[1]+q**2*linear[0])
                -2*s**3*p*q*h)
    assert sp.expand((x*sp.eye(3)-sp.diag(*masses)-s*portal).det()-expected) == 0
    normal = sp.Matrix([b-a, a-c, c-b])
    assert sp.simplify(sp.Matrix([[1, 1, 1], [c, b, a]])*normal) == sp.zeros(2, 1)
    t = sp.Symbol("t", positive=True)
    cleared = normal.dot(sp.Matrix([(t+c)**2, (t+b)**2, (t+a)**2]))
    assert sp.expand(cleared-(b-a)*(c-a)*(c-b)) == 0
    z, m = sp.symbols("z m", positive=True)
    assert sp.integrate(z**2, (z, 0, 1))/m == 1/(3*m)
    # In the radial integral, substitute t=m z/(1-z).
    tx = m*z/(1-z)
    assert sp.cancel(tx**2*sp.diff(tx, z)/(tx+m)**4-z**2/m) == 0


def check_examples():
    x, s = sp.symbols("x s")
    masses = sp.diag(1, 2, 3)
    y = x-2-10*s
    shared = y*(y**2-1-99*s**2)
    determinant = (2+10*s)*(3+40*s+s**2)
    strengths = []
    for q, p in ((3, 9), (5, 7), (7, 1)):
        f = sp.Matrix([[10, q, -sp.I*p], [q, 10, q], [sp.I*p, q, 10]])
        assert f == f.H and f.det() == 10
        assert sp.expand((x*sp.eye(3)-masses-s*f).det()-shared) == 0
        assert sp.expand((masses+s*f).det()-determinant) == 0
        strengths.append(sp.Matrix([q*q, p*p, q*q]))
    assert sp.discriminant(determinant, s) > 0
    assert len(sp.polys.polytools.intervals(determinant, eps=sp.Rational(1, 100))) == 3
    for v in strengths[1:]:
        assert v-strengths[0] == (v[0]-strengths[0][0])*sp.Matrix([1, -2, 1])
    print("Equal characteristic curves for (q,p) = (3,9), (5,7), (7,1).")
    print("All three zeros of det(D+sF) are simple, and the bubble separates all three examples.")


def check_resonance():
    # Compare F=±diag(3,4,-5); near s=0 the positive logarithms split.
    a, b, c = map(sp.Integer, (3, 4, 5))
    s, x = sp.symbols("s x")
    assert a*a+b*b == c*c
    q, r = sp.symbols("q r", positive=True)
    g = lambda m: m*m*(sp.log(m)-sp.Rational(3, 2))
    first = g(a*q)+g(b*q)+g(c*r)
    second = g(a*r)+g(b*r)+g(c*q)
    contrast = a*a*sp.log(a/c)+b*b*sp.log(b/c)
    assert sp.expand(sp.expand_log(first-second-contrast*(q*q-r*r), force=True)) == 0
    assert sp.expand((q*q-r*r).subs({q: 1+s, r: 1-s})-4*s) == 0
    f = sp.diag(a, b, -c)
    d = sp.diag(a, b, c)
    bubbles = [sum(portal[j, j]**2/(6*d[j, j]) for j in range(3)) for portal in (f, -f)]
    assert bubbles == [2, 2]
    assert sp.expand((x*sp.eye(3)-d-s*f).det()-(x*sp.eye(3)-d+s*f).det()) != 0
    print("Resonant masses (3,4,5): equal B=2; the two potentials differ by 4 s K, K<0.")


if __name__ == "__main__":
    check_invariants()
    check_examples()
    check_resonance()
