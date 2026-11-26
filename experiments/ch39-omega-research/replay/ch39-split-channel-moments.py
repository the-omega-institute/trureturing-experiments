"""Exact three-moment separation with an excitation floor and a soft-mode counterexample."""
from fractions import Fraction as Q


NATIVE = Q(624, 625)
NEXT = Q(625, 626)


def corner(radius, floor):
    a, b, c = NATIVE+radius, NATIVE-5*radius, NATIVE+13*radius
    return (a-floor*NEXT)*(c-floor*b)-(b-floor*a)**2


def two_energies(radius, floor):
    a, b = NATIVE+radius, NATIVE-5*radius
    high = (b-floor*a)/(a-floor*NEXT)
    upper_mass = (a-floor*NEXT)/(high-floor)
    return high, upper_mass, NEXT-upper_mass


def tolerance(floor):
    lo, hi = Q(0), Q(1,10**6)
    assert corner(lo,floor) < 0 < corner(hi,floor)
    for _ in range(100):
        mid = (lo+hi)/2
        if corner(mid,floor) < 0:
            lo = mid
        else:
            hi = mid
    # All inequalities needed for attainment hold throughout the root bracket.
    assert NATIVE-5*hi > floor*(NATIVE+hi)
    assert NEXT*(NATIVE-5*hi) > (NATIVE+hi)**2
    assert NATIVE > floor*NEXT
    # A smaller channel count cannot reach the moment cube even in its best corner.
    assert Q(623,624)*(NATIVE+5*hi) < (NATIVE-hi)**2
    print("gap floor",float(floor),"per-reading critical error",float((lo+hi)/4))
    return lo/2, hi/2


def verify():
    soft = Q(1,10**10)
    high, upper, lower = two_energies(0, soft)
    assert 0 < lower < NEXT and 0 < soft < 1 < high
    assert upper*high+lower*soft == upper*high**2+lower*soft**2 == NATIVE
    excess = upper*high**3+lower*soft**3-NATIVE
    assert excess == NATIVE*soft*(1-soft)*(NEXT-NATIVE)/(NATIVE-soft*NEXT)
    lo,hi = tolerance(Q(1,2))
    assert Q("2e-8") < lo < hi < Q("3e-8")
    lo,hi = tolerance(Q("8e-8"))
    assert lo > Q("3e-15")


if __name__ == "__main__":
    verify()
