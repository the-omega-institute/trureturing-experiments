from fractions import Fraction as F
from itertools import permutations
import sympy as sp


class Interval:
    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = self.lo if hi is None else F(hi)
        assert self.lo <= self.hi

    def __add__(self, other):
        other = as_interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)
    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-as_interval(other))

    def __rsub__(self, other):
        return as_interval(other) - self

    def __mul__(self, other):
        other = as_interval(other)
        vals = [self.lo * other.lo, self.lo * other.hi,
                self.hi * other.lo, self.hi * other.hi]
        return Interval(min(vals), max(vals))
    __rmul__ = __mul__

    def __truediv__(self, other):
        other = as_interval(other)
        assert other.lo > 0 or other.hi < 0
        return self * Interval(1 / other.hi, 1 / other.lo)

    def __pow__(self, n):
        out = Interval(1)
        for _ in range(n):
            out = out * self
        return out


def as_interval(x):
    return x if isinstance(x, Interval) else Interval(x)


def cyclotomic_t(t):
    return t**5 + t**4 - 4*t**3 - 3*t**2 + 3*t + 1


# These are rational isolating intervals for the five values 2 cos(2*pi*k/11),
# obtained by sign changes of the displayed integer minimal polynomial.
root_data = [
    (16825070656, 16825070657),
    (8308300260, 8308300261),
    (-2846296766, -2846296765),
    (-13097214679, -13097214678),
    (-19189859473, -19189859472),
]
roots = [Interval(F(lo, 10**10), F(hi, 10**10)) for lo, hi in root_data]
for t in roots:
    assert cyclotomic_t(t.lo) * cyclotomic_t(t.hi) < 0
assert all(roots[i + 1].hi < roots[i].lo for i in range(len(roots) - 1))


R = F(239, 500)


def block_polynomial(a, t):
    e1 = t*t + t - 2
    e2 = t**3 - 2*t + 2
    # det(z I - D_k K_a), in ascending coefficient order.
    return [Interval(a**3),
            -a*a*(3+a)*e1/4,
            a*(1+a)*e2/2,
            -(1+3*a)*e1/4,
            Interval(1)]


def schur_inside(coeffs, radius):
    # Exact interval Schur recursion for p(radius*z).  Every strict margin
    # 1-c_0^2 > 0 certifies that the current factor is Schur stable.
    c = [as_interval(x) * (radius ** j) for j, x in enumerate(coeffs)]
    while len(c) > 1:
        c = [x / c[-1] for x in c]
        margin = Interval(1) - c[0] * c[0]
        if margin.hi < 0:
            return False
        assert margin.lo > 0
        c = [c[j + 1] - c[0] * c[-2-j] for j in range(len(c)-1)]
    return True


for a in (F(1, 5), F(3, 10)):
    for t in roots:
        assert schur_inside(block_polynomial(a, t), R)


# At a=1/4 the middle block has a real root in (-1,-R), by an exact sign
# change.  This is the obstruction to a unimodal three-point gap sequence.
a = F(1, 4)
t = roots[2]
p = block_polynomial(a, t)
def eval_interval(poly, z):
    out = Interval(0)
    for c in reversed(poly):
        out = out * z + c
    return out

left = eval_interval(p, -1)
right = eval_interval(p, -R)
assert left.lo > F(33, 100)
assert right.hi < -F(1, 20000)


# Independent literal determinant check for every nontrivial character of
# Z/11 and S={1,2,9,10}; the congruence is exact in Q[u]/Phi_11.
u, z, a_sym, t_sym = sp.symbols("u z a t")
phi11 = sum(u**j for j in range(11))
for k in range(1, 6):
    speeds = [1, 2, 9, 10]
    diag = [u ** ((k * v) % 11) for v in speeds]
    mat = sp.Matrix(4, 4, lambda i, j:
                    z if i == j else 0)
    for i in range(4):
        for j in range(4):
            mat[i, j] -= diag[i] * (a_sym * int(i == j) + (1-a_sym)/4)
    det = sp.expand(mat.det())
    e1 = (u**k + u**(11-k))**2 + (u**k + u**(11-k)) - 2
    e2 = (u**k + u**(11-k))**3 - 2*(u**k + u**(11-k)) + 2
    expected = (z**4 - (1+3*a_sym)*e1*z**3/4
                + a_sym*(1+a_sym)*e2*z**2/2
                - a_sym**2*(3+a_sym)*e1*z/4 + a_sym**3)
    assert sp.rem(sp.Poly(sp.together(det-expected).as_numer_denom()[0], u),
                  sp.Poly(phi11, u)).is_zero

print("source Fourier-block identity: exact")
print(f"middle block signs: p(-1) > 33/100 and p(-239/500) < -1/20000")
print("endpoint Schur certificates: exact rational interval margins")
print("ALL_CHECKS_OK")
