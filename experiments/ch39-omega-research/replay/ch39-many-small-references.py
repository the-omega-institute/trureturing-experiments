"""Largest threshold over all nonempty subsets of the 2275 additional heights."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
s=run_path(str(root/"ch39-arbitrary-small-reference.py"))
m=s["m"]
B=m["Bounds"]


def extrema(ids):
    prefix=s["data"]("ch39-tick-separators.json")["prefix"][:2323]
    scale=1<<156
    height=1<<400
    cells=[(height**2*scale//((h+1)**2),-(-height**2*scale//(h*h))) for h,*_ in prefix]
    old=(sum(cells[i-1][0] for i in ids),sum(cells[i-1][1] for i in ids))
    remaining=sorted(set(range(1,2324))-set(ids),reverse=True)
    def ratio(total,count):
        return B(Q(48*total[0],old[1]*count),Q(48*total[1],old[0]*count))
    sums=list(old)
    best_by_count=[]
    replacements=[]
    for n,i in enumerate(remaining,1):
        c=cells[i-1]
        sums=[x+y for x,y in zip(sums,c)]
        best_by_count.append(ratio(sums,48+n))
        if n<len(remaining):
            other=cells[remaining[n]-1]
            changed=(sums[0]+other[0]-c[1],sums[1]+other[1]-c[0])
            replacements.append(ratio(changed,48+n))
    at=min(range(len(remaining)),key=lambda i:best_by_count[i].lo)
    best=best_by_count[at]
    second=min([v.lo for i,v in enumerate(best_by_count) if i!=at]+[replacements[at].lo])
    assert at+1==2157 and second-best.hi>Q(1,10**9)
    return best,second,best_by_count[-1],remaining[at],remaining[at+1]


def main():
    ids,allowed,_=s["setup"]()
    best,next_scale,full,first,last=extrema(ids)
    print("Optimal added tail starts at zero",first,"previous missing zero",last)
    print("Best scale",float(best.lo),"runner-up scale lower bound",float(next_scale))
    print("All 2275 extra heights",float(full.lo))
    # For any nonempty addition distribute at most 2^-300 of total weight.
    # The correction of old weights stays in the same boxes irrespective of
    # the number of new sites. Their poles in baseline units all exceed 2.
    neg = Q(1,2**300)
    for mass,k,wanted in (
        (658554352740000,best.hi,1),(658554352750000,best.lo,-1),
        (658553353760000,next_scale,-1),
        (627439358420000,full.hi,1),(627439358440000,full.lo,-1)):
        sign,n,_ = m["endpoint"](mass,ids,allowed,kinetic_scale=k,closed_weight_cap=neg)
        assert (sign,n)==(wanted,144)
    assert (first,last)==(161,160)
    print("Unique largest-mass roster adds all missing indices >=161, total 2157 new sites.")
    print("Mass lies in (658554352740000,658554352750000) GeV; every other added roster is below 658553353760000.")
    print("With the whole prefix filled, mass lies in (627439358420000,627439358440000) GeV.")
    return best,next_scale,full


if __name__=="__main__":
    main()
