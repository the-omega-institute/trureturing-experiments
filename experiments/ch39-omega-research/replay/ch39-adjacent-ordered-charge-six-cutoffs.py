"""相邻零点荷六首枢轴的保序截止及最低电弱尺度的左右差异。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def solve_three(n,abc,bracket,high_pair,low_weak,high_weak,mz,p0,p1):
    a,b,c=abc
    reference=arb("-.074256577194663365545")
    lo,hi=map(arb,bracket)
    assert a>1 and b<0<c
    # 中点迭代只寻找小区间。最终两端以高能两根的严格外包复核。
    for _ in range(23):
        mid=(lo+hi)/2
        hint=a*low_weak(mid)+(1-a)*reference
        if hint>0:
            lo=mid
        else:
            assert hint<0
            hi=mid
    low = p0*lo.union(hi)
    middle,highest = high_pair(abc,low)
    rest=b*high_weak(middle,p0,p1)+c*high_weak(highest,p0,p1)
    dl,dr = [a*low_weak(t)+rest for t in (lo,hi)]
    assert dl>0>dr and 0<low<middle<highest<1
    minimum=mz*lo.union(hi).sqrt()
    print(n,"最小 Q/GeV、另两平方截止、弱差端点",minimum.str(12),
          middle.str(12),highest.str(12),dl.str(6),dr.str(6),flush=True)
    return low,middle,highest,minimum


def verify(root):
    gravity = run_path(str(root / "ch39-three-index-free-atoms.py"))
    native = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    pair = run_path(str(root / "ch39-ordered-cross-gravity-roots.py"))["high_pair"]
    high_weak = run_path(str(root / "ch39-full-neutral-free-cutoff.py"))["weak_column"]
    center = 2_306_888_100_085
    numbers = (center-1,center,center+1)
    moments,heat,fields = gravity["broken_fields"](root,numbers)
    ctx.prec = 256
    mp,mz = arb("1.22089e19"),arb("91.1876")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    low_weak = native["weak_series"]()
    result = []
    for n,rows,bracket in zip(numbers, fields, ((".991","1"),("1","2"),("1","2")),strict=True):
        _,levels,pivots,_ = next(t for t in native["shell_weights"](root,n) if t[0]==6)
        weights = [v.real for v in pivots[0]]
        assert sum(weights).contains(1)
        abc = (sum(weights[:17]),weights[17],weights[18])
        a,b,c = abc
        if n==center-1:
            depth = max((-sum(weights[i:])).upper() for i in range(1,19))
            bound = low_weak(arb(1))+depth*(low_weak(arb(1))-high_weak(arb(1),p0,p1))
            print("左邻所有电弱以上保序截止的弱差上界",bound.str(13),flush=True)
            assert bound<arb("-.00027260")
        def high_pair(abc,low):
            return pair(abc,low,rows,heat,moments)
        if n==center-1:
            promoted=(a-weights[6],b+weights[6],c)
            lo_pos,mid_pos,hi_pos,minimum=solve_three(n,promoted,("1","1.001"),high_pair,low_weak,high_weak,mz,p0,p1)
            assert mz<minimum and lo_pos<mid_pos<hi_pos
            print("左邻将第七层移入中截止后，最低 Q/GeV",minimum.str(12),flush=True)
        low,middle,highest,minimum=solve_three(n,abc,bracket,high_pair,low_weak,high_weak,mz,p0,p1)
        assert (minimum<mz if n==center-1 else minimum>mz)
        result.append((levels,[low]*17+[middle,highest]))
    return result


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
