"""The 624-port realization on two minimally extended shell registers."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent


def verify():
    pole = (17,17,17)
    shell = {p for p in product(range(30),repeat=3) if sum(x*x for x in p)==867}
    def rotate(p):
        return p[1],p[2],p[0]
    near = (13,13,23)
    chosen = (near,rotate(near),rotate(rotate(near)))
    register = [(p,k) for p in sorted(shell) for k in range(3 if p in chosen else 1)]
    def lift(tag):
        p,k = tag
        return rotate(p),(k+int(p==chosen[2]))%3
    outer = ((pole,0),(pole,0))
    todo = set(product(register,repeat=2))-{outer}
    cycles = Counter()
    triangle = ((13,13,23),(7,17,23),(11,11,25))
    assert all(sum((x-y)**2 for x,y in zip(triangle[j-1],triangle[j]))<=56
               for j in range(3))
    while todo:
        pair = min(todo)
        start = pair
        tri = triangle
        length = 0
        while pair in todo:
            todo.remove(pair)
            pair = tuple(map(lift,pair))
            tri = tuple(map(rotate,tri))
            length += 1
        assert pair==start and tri==triangle
        cycles[length] += 1
    assert cycles=={3:85,9:41}
    for k in range(9):
        dimension = sum(count for length,count in cycles.items() if k*length%9==0)
        assert dimension == (126 if k%3==0 else 41)
    d = [int(j==0)-int(j==3) for j in range(9)]
    corr = [F(sum(x*d[(j+k)%9] for j,x in enumerate(d)),2) for k in range(9)]
    assert corr==[1,0,0,F(-1,2),0,0,F(-1,2),0,0]
    ctx.prec = 150
    selection = run_path(str(ROOT/"ch39-rounded-phase-selection.py"))
    n,gain,branch,r,mass = next(selection["matches"]())
    assert (n,gain,branch)==(11490,len(register)**2-1,1)
    zero = acb.zeta_zero(n)
    w = (zero-1)/zero
    assert w.real > arb("0.999") and 0 < w.imag < arb("0.001")
    # Constant mode of w times any permutation has the same solution.
    x = 2*w/(625-623*w)
    average = (1+624*x)/625
    b = 2*average-1
    assert (x-w*(2*average-x)).contains(0)
    assert ((b-1)/(b+1)-624*(w-1)/(w+1)).contains(0)
    assert all(not (1+w**k).contains(0) for k in (3,9))
    print("Two 25-state registers give",cycles,"internal cycles; gain",gain,
          "at zero",n,"r",r.str(10),"bare scale",mass.str(12))


if __name__ == "__main__":
    verify()
