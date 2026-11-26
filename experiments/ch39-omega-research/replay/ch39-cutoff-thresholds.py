"""Integer comparison at the closed-fermion one-loop crossing, for exact positive fits."""
from fractions import Fraction as Q
from functools import cache
import json
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
word = run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
rates = run_path(str(root/"ch39-two-prime-height-law.py"))["rate_numerators"]()
bounds = run_path(str(root/"ch39-rh-clock-integer-replay.py"))
Bounds, read_cell = bounds["Bounds"], bounds["read_cell"]
SCALE = word["DEN"]
HEIGHT = 1 << 400
OFFSETS = [*range(43), 1682]
MODES = (7, 6, 3, 3, 2, 2)
MASS_BRACKETS = ((103800021814678,103800021814680),
                 (103311466292796,103311466292799))


@cache
def data(name):
    return json.loads((root/name).read_text())


def supports():
    case = data("ch39-extra-offset-cutoffs.json")[-1]
    assert case["offsets"] == OFFSETS and case["cut"] == 2323
    original = set(case["support"])
    return [sorted(original-{old}|{new}) for old,new in ((39,2315),(455,29),(39,39))]


def observations(ids):
    prefix = data("ch39-tick-separators.json")["prefix"]
    levels = data("ch39-minimum-tick-finite-inputs.json")["levels"]
    extra = data("ch39-fixed-register-extra-target.json")
    assert extra["j"] == 1682 and extra["bits"] == 288
    heights = [prefix[i-1][0] for i in ids]
    rows = []
    for j in OFFSETS:
        v = (j*HEIGHT)**2
        rows.append([(h*h*SCALE//(h*h+v), -(-((h+1)**2*SCALE)//((h+1)**2+v))) for h in heights])
    angles = [[] for _ in range(4)]
    for h in heights:
        k,l = [-43235114615559591407865398948523*(((h>>112)*r+2**311)//2**312) for r in rates]
        for row,exponent in zip(angles,(k,l,l-k,l+k)):
            row.append(word["power"](exponent)[0])
    targets = [extra["floor"] if j == 1682 else levels[j] for j in OFFSETS]
    targets = [(v,v+int(j!=0)) for j,v in zip(OFFSETS,targets)]
    assert targets[0] == (2**288,2**288)
    return rows+angles, targets+[(0,0)]*4


def exact_fit(record, minimum=Q(4,100000), system=None, size=48, cutoff=2323, target_bits=288, enclose=None):
    ids, bits, v = record["support"], record["bits"], record["inverse"]
    matrix,target = observations(ids) if system is None else system
    n = len(ids)
    assert n == len(set(ids)) == len(matrix) == len(v) == size and max(ids) == cutoff
    den = 2**bits*SCALE
    def linear(row, cells):
        low = high = 0
        for a,(l,h) in zip(row,cells,strict=True):
            low += a*(l if a>=0 else h)
            high += a*(h if a>=0 else l)
        return low,high
    eps = Q(0)
    b = []
    for i,row in enumerate(v):
        rownorm = 0
        for j,col in enumerate(zip(*matrix)):
            l,h = linear(row,col)
            diag = den*int(i==j)
            rownorm += max(abs(l-diag),abs(h-diag))
        eps = max(eps,Q(rownorm,den))
        l,h = linear(row,target)
        b.append(Bounds(Q(l,2**(target_bits+bits)),Q(h,2**(target_bits+bits))))
    assert eps < Q(1,2**100)
    pad = eps*max(max(abs(u.lo),abs(u.hi)) for u in b)/(1-eps)
    enclose = Bounds.positive_scale if enclose is None else enclose
    w = [enclose(Bounds(x.lo-pad,x.hi+pad)) for x in b]
    assert min(x.lo for x in w)>minimum
    return dict(zip(ids,w))


def physical_inputs():
    cells = data("ch39-rh-clock-integer-data.json")["inputs"]
    return (read_cell(cells["pivot"][0]), list(map(read_cell,cells["lambda"])),
            read_cell(cells["q02"][0]), read_cell(cells["target"][0]))


def kinetic(ids):
    prefix = data("ch39-tick-separators.json")["prefix"]
    if max(ids)>len(prefix):
        prefix = prefix + data("ch39-adjacent-offset-certificate.json")["tail"]
    if max(ids)>len(prefix):
        prefix = prefix + data("ch39-paired-products-data.json")["tail"]
    z = [Bounds(Q(HEIGHT**2,(prefix[i-1][0]+1)**2),Q(HEIGHT**2,prefix[i-1][0]**2)).positive_scale()
         for i in ids]
    mean = sum(z,Bounds(0))/Bounds(len(z))
    return [(mean/x).positive_scale() for x in z]


# Fixed-point arithmetic for interval signs of the arrowhead secular function.
# All four quotient corners are rounded out. No point at a pole is evaluated.
UNIT = 1 << 120


def grid(x):
    return (x.lo*UNIT//1,-((-x.hi*UNIT)//1))


def product(a,b):
    corners = [x*y for x in a for y in b]
    return min(corners)//UNIT,-((-max(corners))//UNIT)


def secular(x, poles, rho, pivot, lam, higgs):
    low,high = pivot
    for (l,h),(r,s) in zip(poles,rho,strict=True):
        assert x<l or h<x
        ends = (x-l,x-h)
        a = min(v*UNIT//d for v in (r,s) for d in ends)
        b = max(-((-v*UNIT)//d) for v in (r,s) for d in ends)
        low += x*a//UNIT
        high -= (-x*b)//UNIT
    l,h = product(lam,(low,high))
    return x-h-higgs[1],x-l-higgs[0]


def endpoint(mass, ids, weights, passive_above=False, kinetic_scale=Q(1), closed_weight_cap=Q(0)):
    p,modes,q02,target = physical_inputs()
    assert 0<p.lo<=p.hi<1 and q02.hi < 8*mass**2
    q = q02/Bounds(mass**2)
    poles = [grid(Bounds(kinetic_scale)*c) for c in kinetic(ids)]
    rho = [grid((Bounds(1)-p)*weights[i]) for i in ids]
    assert all(a[1]<b[0] for a,b in zip(poles,poles[1:]))
    for k,(r,s) in enumerate(rho):
        if r<=0:
            # An exactly uncoupled reference above all gaps being tested
            # contributes only closed modes; keep its kinetic coefficient.
            assert passive_above and r==s==0 and k>0 and Q(poles[k-1][0],UNIT)>q.hi
    pure = sum(Q(h,UNIT)<q.lo for l,h in poles)
    assert 0<pure<len(ids) and Q(poles[pure][0],UNIT)>q.hi
    assert Q(poles[0][0]*mass**2,UNIT)>1000**2
    exposure = Bounds(1)
    active = 25*pure
    for l,h in poles[:pure]:
        exposure = (exposure*(q/Bounds(Q(l,UNIT),Q(h,UNIT)))**25).positive_scale()
    secular_p = p
    if closed_weight_cap:
        # All added poles exceed twice the scale and the tested gaps stop below it.
        # Their total contribution to x sum(rho/(x-C)) lies in [-cap,0].
        assert poles[pure][1]<UNIT*kinetic_scale
        secular_p = Bounds(p.lo-closed_weight_cap,p.hi)
    for lam,mult in zip(modes,MODES,strict=True):
        signs = (poles,rho,grid(secular_p),grid(lam),grid(Bounds(Q(125**2,mass**2))))
        for j in range(pure+1):
            lo = poles[j-1][1]+1 if j else 0
            hi = poles[j][0]-1
            # Endpoints belong to the same pole gap. Only confirmed signs
            # move the bound inwards; final checks certify both retained ends.
            for _ in range(90):
                x = (lo+hi)//2
                a,b = secular(x,*signs)
                if b<0:
                    lo=x
                elif a>0:
                    hi=x
                else:
                    break
            assert secular(lo,*signs)[1]<0 and secular(hi,*signs)[0]>0
            rootcell = Bounds(Q(lo,UNIT),Q(hi,UNIT))
            assert rootcell.lo*mass**2>1000**2
            if rootcell.hi<q.lo:
                exposure = (exposure*(q/rootcell)**mult).positive_scale()
                active += mult
            else:
                assert rootcell.lo>q.hi
    ratio = (exposure**11/target).positive_scale()
    sign = 1 if ratio.lo>1 else -1
    assert sign==1 or ratio.hi<1
    return sign,active,ratio


def main():
    records = data("ch39-cutoff-threshold-inverses.json")
    assert [r["support"] for r in records] == supports()
    assert MASS_BRACKETS[0][0]-MASS_BRACKETS[1][1]>488555521878
    for rec,bracket in zip(records[:2],MASS_BRACKETS,strict=True):
        weights = exact_fit(rec)
        for m,expected in zip(bracket,(1,-1)):
            sign,active,ratio = endpoint(m,rec["support"],weights)
            assert (sign,active)==(expected,144)
            print(m,active,sign, float(ratio.lo),float(ratio.hi))


if __name__ == "__main__":
    main()
