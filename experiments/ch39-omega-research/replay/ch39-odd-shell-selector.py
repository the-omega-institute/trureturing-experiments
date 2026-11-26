"""Ground energy of the conjugation-odd shell return and its first even perturbation."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

shell = run_path(str(Path(__file__).with_name("ch39-native-flavor-eigenaxes.py")))
first,second = 11431,11363


def energy(data, eps):
    trace, odd, even = data
    return eps*trace-(odd**2+eps**2*even).sqrt()


def slope(data, eps):
    trace, odd, even = data
    return trace-eps*even/(odd**2+eps**2*even).sqrt()


def verify():
    data = {}
    for idx,z in enumerate(acb.zeta_zeros(10815,1278),10815):
        t = shell["plane_response"](z.imag)
        mean = t.trace().real/2
        even = t[0,1].real**2 + (t[0,0].real-mean)**2
        data[idx] = mean,abs(t[0,1].imag),even
    lead = data[first]
    runner = data[10965]
    assert all(lead[1] > d[1] for idx,d in data.items() if idx!=first)
    assert all(runner[1] > d[1] for idx,d in data.items() if idx not in (first,10965))
    assert lead[1]-runner[1] > arb(".0000636")
    assert lead[1]/runner[1] > arb("1.00258")
    other = data[second]
    left,right = arb(".0349009618863"),arb(".0349009618864")
    assert energy(lead,left) < energy(other,left)
    assert energy(lead,right) > energy(other,right)
    stop = arb(".035")
    intervals = [(arb(0),stop,0)]
    while intervals:
        lo,hi,depth = intervals.pop()
        eps = lo.union(hi)
        change = slope(other,eps)-slope(lead,eps)
        if not change < 0:
            assert depth < 16
            mid = ((lo+hi)/2).mid()
            intervals.extend(((lo,mid,depth+1),(mid,hi,depth+1)))
    # No third spectrum can cross either of these two before the stated endpoint.
    queue = [(arb(0),stop,0,[n for n in data if n not in (first,second)])]
    while queue:
        lo,hi,depth,candidates = queue.pop()
        eps,mid = lo.union(hi),((lo+hi)/2).mid()
        threshold = energy(lead,mid) if mid < left else energy(other,mid)
        chosen = lead if mid < left else other
        # It suffices to beat one of the two for every eps in this interval.
        keep = [n for n in candidates if not (energy(data[n],mid)-threshold+
                (eps-mid)*(slope(data[n],eps)-slope(chosen,eps))) > 0]
        if keep:
            assert depth < 24
            queue.extend(((lo,mid,depth+1,keep),(mid,hi,depth+1,keep)))
    assert energy(lead,arb(0))<0 and lead[0]<0 and other[0]<0
    full = {n:energy(d,arb(1)) for n,d in data.items()}
    assert full[11474]<0 and all(full[11474]<v for n,v in full.items() if n!=11474)
    print("pure odd winner 11431, gap", (lead[1]-runner[1]).str(12),
          "first crossover to 11363 in",left.str(15),right.str(15),"full return 11474")


if __name__ == "__main__":
    verify()
