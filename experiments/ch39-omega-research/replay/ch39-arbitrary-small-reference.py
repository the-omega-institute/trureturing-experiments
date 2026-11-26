"""Any missing height in the optimal prefix can carry arbitrarily little positive mass."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
m=run_path(str(root/"ch39-cutoff-thresholds.py"))
B=m["Bounds"]
data=m["data"]


def setup():
    rec=data("ch39-cutoff-threshold-inverses.json")[1]
    w=m["exact_fit"](rec)
    ids=rec["support"]
    prefix=data("ch39-tick-separators.json")["prefix"]
    missing=[i for i in range(1,2324) if i not in w]
    assert len(missing)==2275 and min(missing)==27
    mean_lower=sum(Q(2**800,(prefix[i-1][0]+1)**2) for i in ids)/len(ids)
    rmax=Q(61,200)
    assert Q(2**800,prefix[26][0]**2)/mean_lower<rmax
    lower=Q(48,49)
    upper=(48+rmax)/49
    assert lower/rmax>2
    # The archived integer preconditioner V has ||VA-I|| < 2^-100.
    norm=Q(max(sum(map(abs,row)) for row in rec["inverse"]),2**rec["bits"])
    assert norm/(1-Q(1,2**100)) < 2**219
    # For any b>=0 take epsilon=2^(-b-300) and decrease the old vector
    # by epsilon * A^{-1} column. This wider bound handles every b at once.
    cap=Q(1,2**81)
    assert (48*2**219+1)<2**300 and min(v.lo for v in w.values())>cap
    allowed={i:B(v.lo-cap,v.hi+cap) for i,v in w.items()}
    return ids,allowed,(lower,upper)


def main():
    ids,allowed,(small,large)=setup()
    # At the low mass use the largest kinetic scale (smallest exposure);
    # at the high mass use the smallest scale (largest exposure).
    for mass,k,sign in ((104000000000000,large,1),(104500000000000,small,-1)):
        s,count,_=m["endpoint"](mass,ids,allowed,kinetic_scale=k,closed_weight_cap=Q(1,2**300))
        assert s==sign and count==144
    assert 104000000000000-m["MASS_BRACKETS"][1][1]>688533707200
    print("At all 2275 missing heights, the 49-site fit lies between 1.040e14 and 1.045e14 GeV.")
    print("An inserted vote of 2^(-b-300) changes every bounded linear reading by less than 2^-b.")


if __name__=="__main__":
    main()
