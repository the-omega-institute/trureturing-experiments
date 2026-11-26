"""Uniform effective flux at positive spectral parameter for both shell clocks."""
from pathlib import Path
import runpy
from flint import arb

model = runpy.run_path(str(Path(__file__).with_name("ch39-shell-schur-flux.py")))
shell = model["shell"]
inside = model["inside"]


def torus_nonzero(poly):
    return abs(poly[0, 0]) > sum(abs(z) for key, z in poly.items() if key != (0, 0))


def all_degree_three(u):
    c = model["schur_data"](u)
    edge = lambda i,j: model["bond"](c, i, j)
    x = edge(7, 6)
    y = edge(11, 12)[0, 0]
    if not (y.real > 0 and 0 < y.imag < y.real):
        return False
    if not all(torus_nonzero(edge(i, j)) for i, j in ((6, 14), (7, 12))):
        return False
    for factor in (1, -y):
        terms = model["imag_terms"](x, factor)
        if not terms[0][0] > sum((a*a+b*b).abs_upper().sqrt() for a, b in terms[1:]):
            return False
    return True


def certify_interval(low, high):
    try:
        if all_degree_three(arb(low).union(arb(high))):
            return 1
    except (AssertionError, ZeroDivisionError):
        pass
    assert high-low > arb(".0001")
    middle = (arb(low)+high)/2
    return certify_interval(low, middle) + certify_interval(middle, high)


if __name__ == "__main__":
    # Bound the high-u path expansion uniformly. The angular triangle inequality and
    # Euler's product for sinh bound any path product by its direct edge factor.
    degrees = [sum(i in (a,b) for a,b,_,_ in shell["edges"]) for i in inside]
    assert min((shell["V"][i] + d).lower() for i,d in zip(inside,degrees)) >= 4
    eta = arb(len(inside))/(6400+4-(len(inside)-1))
    angle_error = 4*eta.asin()
    assert angle_error < arb(".0076")
    print("large-parameter face error", angle_error.str(12))
    count = certify_interval(40, 6400)
    print("all clock phases have degree three for u >= 40; interval cells", count)
