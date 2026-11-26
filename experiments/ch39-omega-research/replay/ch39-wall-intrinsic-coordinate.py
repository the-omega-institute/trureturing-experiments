"""Schwarz-Christoffel coordinates of the seven-vertex flat cone sphere."""
from flint import arb, ctx

ctx.prec = 192
pi = arb.pi()
alpha = (arb(18)/68).sqrt().acos()/pi
gamma = (arb(2)/3).sqrt().acos()/pi + arb(1)/6
delta = arb(1)/3
beta = 2 - alpha - gamma - delta
assert all(0 < a < 1 for a in (alpha, beta, gamma, delta))


def euler_beta(a, b):
    return a.gamma()*b.gamma()/(a+b).gamma()


def side_pole(t):
    return (t**(gamma-1) * euler_beta(alpha, beta)
            * (1/t).hypgeom_2f1(1-gamma, alpha, alpha+beta))


def side_wall(t):
    return ((t-1)**(beta+gamma-1) * t**(alpha-1) * euler_beta(gamma, beta)
            * (1-1/t).hypgeom_2f1(1-alpha, gamma, beta+gamma))


def isolate(func, a, b, steps):
    a, b = arb(a), arb(b)
    assert func(a) < 0 < func(b)
    for _ in range(steps):
        m = (a+b)/2
        value = func(m)
        if value > 0:
            b = m
        elif value < 0:
            a = m
        else:
            raise ArithmeticError("insufficient precision to isolate the coordinate")
    return a.union(b)


target_ratio = (arb(12)/68).sqrt()
t = isolate(lambda t: side_wall(t)/side_pole(t)-target_ratio, "1.0001", 5, 90)
scale = arb(68).sqrt()/side_pole(t)

# On the other meridian replace v by s=v/(1+v). The two binomial series
# in s have positive coefficients, so their product has coefficient <= n+1.
size = 150
def binomial_series(a, c):
    terms = [arb(1)]
    for j in range(1, size):
        terms.append(terms[-1]*c*(a+j-1)/j)
    return terms


first = binomial_series(1-delta, 1)
second = binomial_series(1-gamma, 1-1/t)
integrated = [sum(first[j]*second[n-j] for j in range(n+1))/(n+alpha)
              for n in range(size)]


def opposite_length(x):
    y = x/(1+x)
    value = arb(0)
    for a in reversed(integrated):
        value = value*y+a
    tail = 2*y**size/(1-y)
    return scale*t**(gamma-1)*y**alpha*(value+arb(0).union(tail))


mid_distance = arb(18).sqrt()
centroid_distance = mid_distance + arb(19)/18
x = isolate(lambda x: opposite_length(x)-centroid_distance, ".001", 3, 48)
radius = x**(arb(1)/3)
print("outer vertex w, trapezoid centroid |z|", t.str(18), radius.str(15))
print("cone angle fractions at pole, inner vertex, outer vertex",
      (3*alpha).str(12), beta.str(12), gamma.str(12))
assert 1 < t and 0 < radius < 1
