"""Sharp total native control cost for two preparations that keep the third pole closed."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
model=run_path(str(root/"ch39-clock-schur-mediators.py"))
base,B,root_cell=[model[k] for k in ("base","B","root_cell")]
cross=run_path(str(root/"ch39-kinetic-reading-crossing.py"))
rounded,mul=cross["rounded"],cross["multiply"]


def main():
    old=cross["ReadingBounds"]()
    record=base["data"]("ch39-cutoff-threshold-inverses.json")[1]
    w=base["exact_fit"](record)
    z=old.z[0]
    tail=old.z[1:]
    factor=mul(*(B(1)+B(a)*z for a in cross["A"]))/mul(*(z-v for v in tail))
    leakage=mul(factor,*tail).hi
    sensitivity=factor.hi*sum(t/(Q(1,a)+z.lo) for a,t in zip(cross["A"],old.tails))
    minimum=w[1].lo-2*leakage-cross["TOL"]*sensitivity
    assert minimum>Q(".2166")
    # Third-pole boundary: this archived sign applies to arbitrary positive
    # spectra satisfying the 44 time rows, including their stated tolerances.
    row=next(r for r in base["data"]("ch39-kinetic-reading-roots.json") if r["mass"] is None and r["sign"]==1)
    assert Q(row["c"])==Q(".00082657012849025")
    old.verify_sample(row)
    H=B(125**2)*old.poles[2]/old.q2
    assert H.hi<old.poles[0].lo and ((B(1)-old.p)/(B(1)-H*old.z[0])).hi<1
    word=base["word"];unit=word["DEN"]
    cosine=[B(Q(l,unit),Q(h,unit)) for l,h in [word["power"](n)[0] for n in range(44)]]
    floor=B(minimum)*z*(B(1)-cosine[41])/B(2)
    assert floor.lo>Q(".000834")>Q(row["c"])
    ids,z,p,*_=model["calibration"]()
    weights=model["base_weights"]()
    s=list(map(root_cell,p))
    assert s[0].lo>max(v.hi for v in s[1:]) and s[1].lo>max(v.hi for v in s[2:])
    sine=[B(0)]+[root_cell((B(1)-v)/B(2)) for v in cosine[1:]]
    # Each one of the first 42 increments on height 1 dominates every first
    # increment on another height. Check also the polygon feasibility margin.
    for j in range(1,43):
        assert (s[0]*(sine[j]-sine[j-1])).lo>(s[1]*sine[1]).hi
    assert (s[0]*(sine[43]-sine[42])).hi<(s[1]*sine[1]).lo
    # At most one index has |n_i|>21. The rest all retain at least cos(21 beta/2).
    low_half=root_cell((B(1)+cosine[42])/B(2))
    other_half=root_cell((B(1)+cosine[21])/B(2))
    sum_min=(sum(s,B(0))-s[0])*other_half+s[0]*low_half
    assert sum_min.lo>2*s[0].hi and (sum(s,B(0))-s[0]).lo>2*s[0].hi
    c40=weights[1]*z[0]*(B(1)-cosine[40])/B(2)
    assert Q(".0008107")<c40.lo<=c40.hi<Q(".0008108")
    # Interpolation from a low witness to the diagonal action can hit this
    # definite 144-doublet model on the same heights and weights.
    target=B(Q(".00082249054717514"))
    mean=sum(z,B(0))/B(len(z))
    parent=run_path(str(root/"ch39-phase-support-parent.py"))["parent_endpoint"]
    assert sum((weights[i]*zi for i,zi in zip(ids,z)),B(0)).lo>target.hi>c40.hi
    for m,sg in ((69725800000000,1),(69725820000000,-1)):
        assert parent(m,ids,[weights],kinetic_scale=target/mean)==(sg,144)
    print("First-height weight exceeds",float(minimum),"for every positive 44-row fit.")
    print("Two preparations, 41 net shell turns: universal closure at calibrated diagonal kinetics.")
    print("Every scheme of at most 40 net turns, at arbitrary locations and mixture probability, admits 144.")


if __name__=="__main__":main()
