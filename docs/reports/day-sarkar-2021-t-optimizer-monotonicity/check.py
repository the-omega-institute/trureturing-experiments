"""Check of the proposed refutation of Day-Sarkar, 'On a Conjecture of Nagy on Extremal Densities'
(arXiv:1910.13465, SIAM J. Discrete Math. 35 (2021) 294-306), Section 5, Conjecture 2:
'For each graph G_l, we have that q(G_l, beta) is an increasing function of beta.'

T(q) has classes Y (clique, size sqrt(beta) q), R (clique, size 1 - sqrt(1 - beta(1-q^2))) and B
(independent, the rest); R is joined to Y and B; there are no Y-B edges.  By the source's Lemma 2.1,
t(G_l, T(q)) = sum over {0,1/2,1}-fractional independence weightings phi of y^{y_phi} r^{r_phi} b^{b_phi}.
H_k: triangle a,b,c; vertex d joined to a; k leaves joined to d.

(1) Closed form D_k = r z^3 + 3 b r^3 + y z^(k+3) + 2 b y r^2 z^k + b r^(k+1) beta (z = r + y) equals the
    weighting sum (enumerated over all 3^(k+4) colourings), k = 2..5, at sample (beta, q).
(2) Exact certificate for H_3 in the coordinate x = z/sqrt(beta): derivative factorisation
    F'(x) = (s x^2 + s - 2x) Q_s(x) / (16 (s x - 1)^4), and sign changes of Q_s at the isolating rationals
    for s = 1/1000 (beta = 10^-6) and s = 1/100 (beta = 10^-4); the induced q-intervals are disjoint and
    reversed (q(beta = 10^-6) > q(beta = 10^-4)).
(3) Numerical global maximiser over a 4000-point grid refined by Newton, k = 2..5: for k >= 3 the
    maximiser matches sqrt(3/5) - (7/25) sqrt(beta) to O(beta) and decreases on (0, 10^-2]; for k = 2 (the
    source's G_6) it increases there.
"""
import itertools
import mpmath as mp
import sympy as sp
mp.mp.dps = 40
bad = 0

def parts(beta, q):
    s = mp.sqrt(beta); y = s*q; r = 1 - mp.sqrt(1 - beta*(1 - q*q)); return y, r, 1 - r - y

def D(k, beta, q):
    y, r, b = parts(beta, q); z = r + y
    return r*z**3 + 3*b*r**3 + y*z**(k+3) + 2*b*y*r**2*z**k + b*r**(k+1)*beta

def Dsum(k, beta, q):
    y, r, b = parts(beta, q)
    E = [(0, 1), (1, 2), (2, 0), (0, 3)] + [(3, 4+i) for i in range(k)]
    w = {0: r, 1: y, 2: b}; tot = mp.mpf(0)
    for ph in itertools.product((0, 1, 2), repeat=4+k):
        if all(ph[u] + ph[v] <= 2 for u, v in E):
            t = mp.mpf(1)
            for c in ph: t *= w[c]
            tot += t
    return tot

for k in (2, 3, 4, 5):
    for beta, q in (("0.01", "0.6"), ("0.0001", "0.77"), ("0.3", "0.2")):
        d = abs(D(k, mp.mpf(beta), mp.mpf(q)) - Dsum(k, mp.mpf(beta), mp.mpf(q)))
        if d > mp.mpf(10)**-30: bad += 1
print("(1) closed form checked, bad so far", bad)

s, x = sp.symbols('s x', positive=True)
r = s**2*(1 - x**2)/(2*(1 - s*x)); z = s*x; y = z - r; b = 1 - z
F = sp.together((r*z**3 + 3*b*r**3 + y*z**6 + 2*b*y*r**2*z**3 + b*r**4*s**2)/s**5)
Q = (28*s**5*x**8 - 104*s**4*x**7 + (27*s**5 + 132*s**3)*x**6 - (74*s**4 + 48*s**2)*x**5
     + (-3*s**5 + 44*s**3 - 28*s)*x**4 + (28*s**4 + 12*s**2 + 20)*x**3 - (7*s**5 + 16*s**3)*x**2
     - (2*s**4 + 12*s**2 + 12)*x + 3*s**5 + 12*s)
if sp.simplify(sp.diff(F, x) - (s*x**2 + s - 2*x)*Q/(16*(s*x - 1)**4)) != 0: bad += 1
qs = x - s*(1 - x**2)/(2*(1 - s*x))
iv = {}
for sv, a, c in ((sp.Rational(1, 1000), sp.Rational(77451, 100000), sp.Rational(77452, 100000)),
                 (sp.Rational(1, 100), sp.Rational(77384, 100000), sp.Rational(77385, 100000))):
    if not (Q.subs({s: sv, x: a}) < 0 < Q.subs({s: sv, x: c})): bad += 1
    iv[sv] = (qs.subs({s: sv, x: a}), qs.subs({s: sv, x: c}))
lo1, hi1 = iv[sp.Rational(1, 1000)]; lo2, hi2 = iv[sp.Rational(1, 100)]
if not lo1 > hi2: bad += 1
print("(2) q(1e-6) in", (float(lo1), float(hi1)), " q(1e-4) in", (float(lo2), float(hi2)), " bad so far", bad)

def argmax(k, beta):
    f = lambda q: D(k, beta, q); N = 4000
    i = max(range(N+1), key=lambda j: f(mp.mpf(j)/N))
    if i in (0, N): return mp.mpf(i)/N
    q = mp.findroot(lambda q: mp.diff(f, q), mp.mpf(i)/N)
    return q if 0 <= q <= 1 and f(q) >= f(mp.mpf(i)/N) else mp.mpf(i)/N
betas = ["1e-8", "1e-6", "1e-4", "1e-3", "1e-2"]
for k in (2, 3, 4, 5):
    qv = [argmax(k, mp.mpf(bt)) for bt in betas]
    dec = all(qv[i] > qv[i+1] for i in range(len(qv)-1)); inc = all(qv[i] < qv[i+1] for i in range(len(qv)-1))
    pred = mp.sqrt(mp.mpf(3)/5) - mp.mpf(7)/25*mp.mpf("1e-3")
    print("(3) k", k, [mp.nstr(v, 10) for v in qv], "decreasing" if dec else ("increasing" if inc else "mixed"))
    if k >= 3 and (not dec or abs(qv[1] - pred) > mp.mpf("1e-5")): bad += 1
    if k == 2 and not inc: bad += 1
print("bad=", bad)
