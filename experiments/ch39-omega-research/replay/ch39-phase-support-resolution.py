"""Finite resolution for identifying the explicit 22-site off-roster perturbation."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
amb = run_path(str(root/"ch39-phase-support-ambiguity.py"))
base, B, multiply, rounded = [amb[k] for k in ("base","B","multiply","rounded")]
mass = Q(1,256)


def tables():
    direction, x, coefficients, total, offsets = amb["direction"]()
    fits = [base["exact_fit"](r) for r in base["data"]("ch39-cutoff-threshold-inverses.json")]
    reference = {i: sum((f.get(i,B(0)) for f in fits), B(0))/B(3) for i in set.union(*(set(f) for f in fits))}
    moved = {i: rounded(reference.get(i,B(0)) + B(mass)*v) for i,v in direction.items()}
    assert all(v.lo > 0 for v in moved.values())
    limits = {i: reference[i]/(B(0)-v) for i,v in direction.items() if v.hi < 0}
    assert all(limits[2315].hi < val.lo for i,val in limits.items() if i != 2315)
    assert Q("0.00615858") < limits[2315].lo <= limits[2315].hi < Q("0.00615859")
    return reference, moved, x, coefficients, total, [j*j for j in offsets if j]


def times_linear(poly, left, right):
    # Homogeneous coefficients in (1-t), t, with outward rounding.
    p = [B(0)]*(len(poly)+1)
    for k,c in enumerate(poly):
        p[k] += c*left
        p[k+1] += c*right
    return [rounded(v) for v in p]


def positive_bernstein(row, depth):
    if min(v.lo for v in row) > 0:
        return 1
    assert depth > 0
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = [rounded((u+v)*B(Q(1,2))) for u,v in zip(row,row[1:])]
        left.append(row[0])
        right.append(row[-1])
    return positive_bernstein(left,depth-1) + positive_bernstein(right[::-1],depth-1)


def main():
    ref, moved, x, coefficients, total, roots = tables()
    assert len(ref) == 50 and len(moved) == 72
    def absolute(y):
        assert y.lo*y.hi > 0
        return y if y.lo > 0 else B(-y.hi,-y.lo)
    def difference(j):
        t = B(Q(-1,j*j))
        h = t/x[-1]-B(1)
        v = B(0)
        for c in reversed(coefficients):
            v = rounded(v*h+c)
        return absolute(multiply(B(mass),v,*(B(1)+B(a)*t for a in roots)) /
                        multiply(total,B(j*j),*(t-y for y in x)))
    for j,bits in [(66,271),(4096,37)]:
        v=difference(j)
        assert Q(2)**(1-bits)<v.lo <= v.hi<Q(2)**(2-bits)
        print("Offset",j,"absolute difference",float(v.lo),"least symmetric precision",bits)
    # Bind the supremum on all positive real a by compactifying a=L*t/(1-t).
    # Polynomial numerator is a^5 Q(-1-1/(a*x_last)) times prod(a-j^2).
    L=1 << 22
    translated=[sum((B((-1)**i*comb(i,k))*coefficients[i] for i in range(k,5)),B(0))
                for k in range(5)]
    numerator=[B(0)]*6
    for k,c in enumerate(translated):
        divisor = x[-1]**k if k else B(1)
        numerator[5-k]=rounded(B(mass*L**(5-k))*c/(total*divisor))
    for a in roots:
        numerator=times_linear(numerator,B(-a),B(L))
    numerator=times_linear(numerator,B(1),B(0))
    denominator=[B(1)]
    for v in x:
        denominator=times_linear(denominator,B(1),B(L)*v)
    assert len(numerator)==len(denominator)==73
    for sign in (-1,1):
        bern=[rounded((d*B(Q(1,2**35))+B(sign)*n)/B(comb(72,i)))
              for i,(d,n) in enumerate(zip(denominator,numerator))]
        print("All-spacing error bound, sign",sign,"positive Bernstein intervals",positive_bernstein(bern,14))
    print("No positive offset separates this pair at 36 bits. Offset 4096 succeeds at 37.")


if __name__=="__main__":
    main()
