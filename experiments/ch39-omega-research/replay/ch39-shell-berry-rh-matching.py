"""Hermitian inertia checks and the common-scale test with unit spectral normalization."""
from itertools import combinations, product
from math import isqrt
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb, arb, ctx


core = run_path(str(Path(__file__).with_name("ch39-shell-rh-arc-laplacian.py")))
Z, target = core["Z"], core["target"]
ctx.prec = 192


def covariant_laplacian(r,cutoff):
    points = [p for p in product(range(isqrt(r)+1),repeat=3) if sum(x*x for x in p)==r]
    matrix = [[acb(0) for p in points] for q in points]
    for i,j in combinations(range(len(points)),2):
        x,y,z = points[i]
        a,b,c = points[j]
        if r-x*a-y*b-z*c>cutoff:
            continue
        overlap = acb(r+x*a+y*b+z*c+arb(r).sqrt()*(z+c), x*b-y*a)
        phase = overlap/abs(overlap)
        matrix[i][i] += 1
        matrix[j][j] += 1
        matrix[i][j] = -phase
        matrix[j][i] = -phase.conjugate()
    return matrix


def count_below(matrix,cut):
    # Hermitian elimination preserves inertia. At the two rational cuts around
    # each suggested cluster every leading pivot must have a certified sign.
    a = [row.copy() for row in matrix]
    for i in range(len(a)):
        a[i][i] -= cut
    negative = 0
    for k in range(len(a)):
        pivot = a[k][k].real
        assert not pivot.contains(0)
        negative += int(pivot<0)
        for i in range(k+1,len(a)):
            term = a[i][k]/pivot
            for j in range(k+1,i+1):
                a[i][j] -= term*a[j][k].conjugate()
    return negative


def eigenvalue_boxes(matrix):
    approximate = np.array([[complex(v) for v in row] for row in matrix])
    seeds = np.linalg.eigvalsh(approximate)
    # The integer endpoints use units of 1e-12. Validity comes from exact counts.
    cuts = []
    for value in seeds:
        mid = round(value*10**12)
        if cuts and mid-8 <= cuts[-1][1]:
            cuts[-1][1] = mid+8
        else:
            cuts.append([mid-8,mid+8])
    out = []
    for lo,hi in cuts:
        a,b = arb(int(lo))/10**12,arb(int(hi))/10**12
        count = count_below(matrix,b)-count_below(matrix,a)
        assert count>0 and a>0
        out.append((a,b,count))
    assert sum(n for _,_,n in out)==len(seeds)
    return out


def stiffness(r,s,boxes):
    def count(k,end):
        return sum(n*Z((s*s+r+k*row[end]).sqrt()) for row in boxes for n in [row[2]])/Z(s)
    # Upper eigenvalue bounds reach the target at smaller k.
    lower = core["increasing_root"](lambda k: target-count(k,1), 10**8, 42)
    upper = core["increasing_root"](lambda k: target-count(k,0), 10**8, 42)
    return lower.union(upper)


if __name__ == "__main__":
    # With g4^{-2} equal to the unscaled trace, alpha_U^{-1}=4*pi*K4*Z(s).
    coeff = 4*arb.pi()*2140/13
    s_lo = core["increasing_root"](lambda s: arb("40.9906")-coeff*Z(s),64)
    s_hi = core["increasing_root"](lambda s: arb("40.6142")-coeff*Z(s),64)
    print("unit-normalized common scale boundaries",s_lo.str(13),s_hi.str(13))
    enclosing = []
    for r,d in ((17,3),(26,3),(29,3),(867,76)):
        spec = eigenvalue_boxes(covariant_laplacian(r,d))
        initial = stiffness(r,arb(0),spec)
        allowed = stiffness(r,s_lo.lower(),spec).union(stiffness(r,s_hi.upper(),spec))
        enclosing.append(allowed)
        print("r",r,"lowest eigenvalue",spec[0][0].union(spec[0][1]).str(8),
              "critical stiffness",initial.str(10),"allowed",
              allowed.lower().str(12),allowed.upper().str(12))
    for a,b in combinations(enclosing,2):
        assert not a.overlaps(b)
