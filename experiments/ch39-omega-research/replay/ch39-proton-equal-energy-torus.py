"""荷五八站在任意两钟相位及小辐边亏额下，等能量非负电流的静态振幅下界。"""
from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    clock = run_path(str(root/"ch39-native-light-projector-pump.py"))
    ctx.prec = 96
    base, links = clock["clock_matrices"](shell, 5)
    adjoint = clock["adjoint"]
    soft = shell["soft"]
    hmin = acb.zeta_zero(1).imag
    offset = soft["point"]["study"]["offset"]
    assert offset > 0 and arb(23)/25*arb(2)**-20/(hmin*hmin+offset) < arb("1e-8")
    sites = (0, 2, 3, 5, 6, 9, 11, 13)
    assert base[6,6].real == 18 and 6 not in (*soft["star"]["orbit"], soft["pole"])
    assert max(base[i,i].real for i in sites)<=18 and max(base[i,i].real for i in (3,5))<=10
    pairs = [(i,j, arb("-.0012") if {i,j} in ({6,3},{6,5}) else arb(".0052"))
             for i,j in combinations_with_replacement(sites,2)]
    squared = lambda z: z.real**2 + z.imag**2

    @lru_cache(None)
    def layer(which, n, k):
        phase = acb(0,arb.pi()*(2*k+1)/n).exp()
        value = links[which]*phase
        return value + adjoint(value)

    @lru_cache(None)
    def perturbation(n):
        return 4*arb(2).sqrt()*(arb.pi()/(2*n)).sin() + 6*arb("1e-8")

    n0 = 16
    cells = [(n0,i,j) for i in range(n0) for j in range(n0)]
    counts = Counter()
    while cells:
        n, i, j = cells.pop()
        green = (base+layer(0,n,i)+layer(1,n,j)).inv()
        error = perturbation(n)
        denominator = 1-error*sum(squared(z) for z in green.entries()).sqrt()
        if denominator > 0:
            norms = {k:sum(squared(green[m,k]) for m in range(19)).sqrt() for k in sites}
            if all(green[k,l].real - error*norms[k]*norms[l]/denominator > floor
                   for k,l,floor in pairs):
                counts[n] += 1
                continue
        assert n < 1024
        cells.extend((2*n,2*i+di,2*j+dj) for di in (0,1) for dj in (0,1))
    assert sum(count*(1024//n)**2 for n,count in counts.items()) == 1024**2
    c = (arb(".0052")-arb(".0001"))/arb(".0064")
    pmin = c-(1-c)**2/(2*c)
    assert 19*pmin**2-1>10 and pmin > arb(7)/10
    at_zero = arb(13)/16
    dominant = (at_zero + (2*at_zero-1).sqrt())/2
    energy_gap = 19*dominant**2 - 11
    assert energy_gap > arb("1.20669")
    print("精确相消所要求的两端能量差下界", energy_gap.lower().str(13))
    print("两钟全环面各层终格数、假设小实部时的单向质量下界", sorted(counts.items()), pmin.str(9))


if __name__ == "__main__":
    verify()
