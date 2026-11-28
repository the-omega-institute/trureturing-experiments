"""Generate the exact center jets and interval angle grid, then verify both simplex charts."""

from itertools import combinations
from pathlib import Path
import subprocess
import tempfile

from flint import arb, ctx
import sympy as sp

SCALE = 1024
BINS = 4096
UNITS = 2**40


def integer_code(expr, bounds):
    """C++ integer expression and a bound on every intermediate's absolute value."""
    if expr.is_Integer:
        return str(expr), abs(int(expr))
    if expr.is_Symbol:
        return str(expr), bounds[expr]
    args = [integer_code(arg, bounds) for arg in expr.args]
    values, limits = zip(*args)
    if expr.is_Add:
        code, limit = "(" + " + ".join(values) + ")", sum(limits)
    elif expr.is_Mul:
        code, limit = "(" + " * ".join(values) + ")", 1
        for bound in limits:
            limit *= bound
    elif expr.is_Pow and expr.exp.is_Integer and expr.exp >= 0:
        power = int(expr.exp)
        code, limit = "(" + " * ".join([values[0]] * power) + ")", limits[0]**power
    else:
        raise ValueError(expr)
    assert limit < 2**53, (expr, limit)
    return code, limit


def generate(header):
    edges = list(combinations(range(5), 2))
    q = sp.symbols("q:10")
    distance = sp.zeros(5)
    for (i, j), value in zip(edges, q):
        distance[i,j] = distance[j,i] = value
    twice_gram = sp.Matrix(4, 4, lambda i, j: distance[0,i+1] + distance[0,j+1] - distance[i+1,j+1])
    poly = sp.Poly(twice_gram.det(method="berkowitz"), *q)
    partials = [poly.diff(x) for x in q]
    second = [p.diff(y) for i, p in enumerate(partials) for y in q[i:]]
    numbers = [poly, *partials, *second]
    temps, outputs = sp.cse([p.as_expr() for p in numbers], symbols=sp.numbered_symbols("u"), optimizations="basic")
    bounds = dict.fromkeys(q, SCALE)
    print("void jets(const Vector& center, std::array<double,66>& result) {", file=header)
    for i, x in enumerate(q):
        print(f"const int64_t {x} = static_cast<int64_t>(center[{i}]*{SCALE});", file=header)
        print(f"assert({x} == center[{i}]*{SCALE});", file=header)
    for name, value in temps:
        code, bounds[name] = integer_code(value, bounds)
        print(f"const int64_t {name} = {code};", file=header)
    for i, (expression, polynomial) in enumerate(zip(outputs, numbers)):
        code, _ = integer_code(expression, bounds)
        denominator = 16 * SCALE**polynomial.total_degree()
        print(f"result[{i}] = static_cast<double>({code})/{denominator}.0;", file=header)
    print("}", file=header)

    # The third derivatives are linear. Their extrema over [0,1]^10 occur
    # by choosing opposite endpoints for the positive and negative coefficients.
    print("constexpr double third_bound[10][10][10] = {", file=header)
    for first in partials:
        print("{",file=header)
        for y in q:
            row = []
            for z in q:
                terms = first.diff(y).diff(z).coeffs()
                bound = max(sum(c for c in terms if c>0), -sum(c for c in terms if c<0))
                row.append(str(float(bound)/16))
            print("{"+",".join(row)+"},",file=header)
        print("},",file=header)
    print("};",file=header)

    ctx.prec = 100
    print(f"constexpr int bins = {BINS};",file=header)
    print("constexpr double atan_lower[bins+1] = {",file=header)
    for i in range(BINS+1):
        angle = arb.pi()/2 if i==BINS else (arb(i)/(BINS-i)).atan()
        fixed = (angle*UNITS).lower().floor().unique_fmpz()
        assert fixed is not None and angle >= arb(fixed)/UNITS and angle <= arb(fixed+1)/UNITS
        print(f"{fixed}.0/{UNITS}.0,",file=header)
    print("};",file=header)
    print(f"constexpr double angle_unit = 1.0/{UNITS}.0;",file=header)


def main():
    source = Path(__file__).with_suffix(".cpp")
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        with (work/"regge-generated.h").open("w") as header:
            generate(header)
        executable = work/"regge-upper"
        subprocess.run(["c++","-std=c++17","-O3","-fopenmp","-fno-fast-math","-ffp-contract=off",
                        "-I",str(work),str(source),"-o",str(executable)],check=True)
        subprocess.run([str(executable)],check=True)


if __name__ == "__main__":
    main()
