"""前两个算术时间读数迫使一个正概率的高度带偏离参考味相位。"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from runpy import run_path

from flint import acb,arb,arb_poly,ctx


root=Path(__file__).parent
source=run_path(str(root/"ch39-flavor-scale-time-prefix.py"))
flavor=source["determinant"]["invariant"]
lower_height,upper_height=10**6,4*10**6


def frequencies():
    geometry=flavor["flavor"]["geom"]
    speed=geometry["a"]-geometry["b"]
    source_height=acb.zeta_zero(8359).imag
    origins=[source_height+angle/speed for angle in flavor["reference_rotations"]()]
    return speed,2*arb.pi()/abs(speed),origins


def exclude_block(block):
    ctx.prec=128
    origin_index,first,last=block
    _,period,origins=frequencies()
    origin=origins[origin_index]
    width=arb("1.4e-7")
    for j in range(first,last):
        center=origin+j*period
        band=center+arb(0,width)
        if abs(acb(arb(1)/2,band).zeta())>0:
            continue
        # 两个外侧端点间的完整计数相等即排除整个小窗。
        left=(center-2*width).zeta_nzeros().unique_fmpz()
        right=(center+2*width).zeta_nzeros().unique_fmpz()
        assert left is not None and left==right,(origin_index,j)
    return last-first


def verify():
    ctx.prec,ctx.cap=1536,3
    times=source["clock"]["clock"]["times"]()
    moments=source["clock"]["clock"]["cleared_moments"](times,2)
    a,b=source["determinant"]["scale_interval"]()
    lower=(200*(b/upper_height)**2).upper()
    upper=(200*(a/lower_height)**2).lower()
    polynomial=-arb_poly([-lower,1])*arb_poly([-upper,1])
    score=sum(polynomial[i]*moments[i] for i in range(3))
    norm=sum(abs(v) for v in source["clock"]["residues_of"](polynomial,source["clock"]["partial_fractions"](2)))
    positive_max=(upper-lower)**2/4
    forced=(score-arb("1e-7")*norm)/positive_max
    assert forced>arb(".20006") and score/norm>arb("1.0614e-5")

    ctx.prec=192
    speed,period,origins=frequencies()
    fc=flavor["flavor"]
    c0=fc["singlet"](11661)[1]
    axes=[fc["singlet_at"]({724:z,804:acb(1)})[1] for z in (acb(1),acb(0,1))]
    amp2=sum(v*v for v in axes)
    tolerance=arb("1e-9")
    slope=abs(speed)*(amp2-(abs(c0)+tolerance)**2).sqrt()
    assert slope*arb("1.4e-7")>tolerance
    blocks=[]
    for i,o in enumerate(origins):
        first=int(((lower_height-o)/period).floor().unique_fmpz())
        last=int(((upper_height-o)/period).ceil().unique_fmpz())
        blocks.extend((i,j,min(j+1024,last+1)) for j in range(first,last+1,1024))
    with ProcessPoolExecutor(max_workers=4) as pool:
        count=sum(pool.map(exclude_block,blocks))
    assert count==53170
    trace,pair,_,_=flavor["invariants"]()
    assert pair[2]<0 and trace[1]<0
    drift=abs(pair[2])*forced*tolerance**2
    assert drift>arb("4.00e-35")
    print("高度带最低概率、覆盖的小窗数",forced.str(13),count)
    print("精确参考味支撑的时间误差下界，两个质量不变量组合的亏额", (score/norm).str(13),drift.str(13))


if __name__=="__main__":
    verify()
