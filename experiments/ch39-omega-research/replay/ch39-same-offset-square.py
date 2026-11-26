"""Second moment of the 1682 time kernel on exactly the same five-tick shell clock."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
base = run_path(str(root/"ch39-cutoff-thresholds.py"))
adj = run_path(str(root/"ch39-adjacent-range.py"))
data = base["data"]
SCALE = base["SCALE"]
HEIGHT = base["HEIGHT"]
J = 1682


def square_system(ids, prefix, target, paired=False):
    # All selected 400-bit rows, including the tail, use the identical phase law.
    angles = [[] for _ in range(4)]
    heights = [prefix[i-1][0] for i in ids]
    for h in heights:
        k,l = [-43235114615559591407865398948523*(((h>>112)*r+2**311)//2**312) for r in base["rates"]]
        for row, m in zip(angles,(k,l,l-k,l+k)):
            row.append(base["word"]["power"](m)[0])
    rows = []
    offsets = [*base["OFFSETS"],1683] if paired else base["OFFSETS"]
    for j,power in [*((j,1) for j in offsets),(J,2)]:
        v = (j*HEIGHT)**2
        rows.append([(h**(2*power)*SCALE//(h*h+v)**power,
                     -(-((h+1)**(2*power)*SCALE)//((h+1)**2+v)**power)) for h in heights])
    levels = data("ch39-minimum-tick-finite-inputs.json")["levels"][:43]
    more = [data("ch39-fixed-register-extra-target.json")["floor"]]
    if paired:
        more.append(data("ch39-adjacent-offset-certificate.json")["target1683"])
    more.append(target)
    rhs = [(v,v+int(i!=0)) for i,v in enumerate(levels+more)]
    return rows+angles,rhs+[(0,0)]*4


def verify(paired=False):
    record = data("ch39-paired-products-data.json" if paired else "ch39-same-offset-square-data.json")
    case = record["case"]
    cutoff = 4691 if paired else 3492
    offsets = [*base["OFFSETS"],1683] if paired else base["OFFSETS"]
    size = len(offsets)+5
    assert case["cut"]==cutoff and case["offsets"]==offsets+[J]
    prefix = data("ch39-tick-separators.json")["prefix"]+data("ch39-adjacent-offset-certificate.json")["tail"]
    if paired:
        prefix = prefix+record["tail"]
    assert all(a[0]+1<b[0] for a,b in zip(prefix,prefix[1:]))
    system = square_system(case["support"],prefix,record["square_floor"],paired)
    weights = base["exact_fit"](case,Q(1,100000),system,size=size,cutoff=cutoff)
    assert len(weights)==size
    unit = 1<<40
    coeff = case["a"]
    c = case["c"]
    pd = 1<<case["depth"]
    assert len(coeff)==len(offsets)+1 and len(c)==4 and len(prefix)==cutoff and prefix[-1][0]+1 < 2**415
    rhs = system[1][:-4]
    level = sum(Q(a*(l if a>=0 else r),unit*2**288) for a,(l,r) in zip(coeff,rhs))
    cost = Q(sum(map(abs,coeff)),unit)+Q(sum(map(abs,c)),pd)
    cap = Q(1,1000000)
    assert level>cap+cost/2**272
    grid = adj["GRID"]
    cheb = adj["cheb"]
    phase_pad = Q(sum(map(abs,c)),pd*2**96)
    for h,*angles in prefix[:-1]:
        e = sum(adj["time_term"](j,h,a) for j,a in zip(offsets,coeff[:-1]))
        z = (h+int(coeff[-1]>0))**2
        e -= (-coeff[-1]*z*z*grid)//(z+(J*HEIGHT)**2)**2
        upper = Q(e,unit*grid)+phase_pad
        for a,t in zip(c,angles):
            l,r = cheb(t,5)
            upper += Q(a*(r if a>=0 else l),pd*2**72)
        assert upper<cap
    # Explicit finite vote table, checked for each target and for the phase law itself.
    counts = case["weights"]
    assert min(counts)>0 and sum(counts)==2**288 and len(counts)==size
    for row,(yl,yh) in zip(*system,strict=True):
        lo = sum(w*l for w,(l,r) in zip(counts,row))
        hi = sum(w*r for w,(l,r) in zip(counts,row))
        assert max(yh*SCALE-lo,hi-yl*SCALE)<SCALE*2**16
    if paired:
        A,D = 1682**2,1683**2-1682**2
        first = data("ch39-fixed-register-extra-target.json")["floor"]
        second = data("ch39-adjacent-offset-certificate.json")["target1683"]
        target_lo = Q((A+D)*second-A*(first+1),D)
        target_hi = Q((A+D)*(second+1)-A*first,D)
        low=high=Q(0)
        for i,n in zip(case["support"],counts):
            h=prefix[i-1][0]
            def product(u):
                z=Q(HEIGHT**2,u*u)
                return n/((1+A*z)*(1+(A+D)*z))
            low+=product(h)
            high+=product(h+1)
        assert max(target_hi-low,high-target_lo)<2**16
    print("Squared probe, paired =",paired,"least highest index",cutoff,"positive reference heights",size)
    print("RH second moment",float(Q(record["square_floor"],2**288)))


if __name__ == "__main__":
    run_path(str(root/"ch39-prime-clock-height-horizon.py"))["main"](5,15)
    import sys
    verify("--paired" in sys.argv)
