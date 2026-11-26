"""三个物理判据的首枚联合候选对算术权与共同尺度的全谱限制。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


local = run_path(str(Path(__file__).with_name("ch39-three-mode-time-localization.py")))


def verify():
    ctx.prec,ctx.cap=1536,3
    times=local["clock"]["times"]()
    (first,window),_=local["required_windows"](times)
    total=local["mask"]["time_data"]()[0]
    ctx.prec=192
    hit,next_height=[acb.zeta_zero(n).imag for n in (8359,8360)]
    target_weight=1/(hit**2*total)
    # 该有理高度严格位于首个实例与下一枚零点之间。
    cutoff=((hit+next_height)/2).mid()
    assert hit<cutoff<next_height and cutoff.zeta_nzeros().unique_fmpz()==8359
    main=cutoff/(2*arb.pi())
    tail=main.log()/(arb.pi()*cutoff)+(arb(7)/8-8359)/cutoff**2
    error=arb(".28")*(cutoff.log()+arb(1)/2)/cutoff**2
    probability=target_weight+(tail+error)/total
    assert target_weight>arb("6.0442e-7") and probability<arb(".006682")
    print("联合实例的标准单点概率、计入全尾的集合上界",target_weight.str(13),probability.str(13))
    ctx.prec=1536
    cap=(200*window/hit**2).upper()
    max_probability=local["weight_cap"](times,cap,41,[cap],arb(2)**-176)
    assert max_probability<arb(".038541")
    floor=hit/first*window.lower().sqrt()
    assert floor>arb("598.6724")
    print("小尺度任意权重的联合概率上界、纯联合支撑的共同尺度下界",max_probability.str(13),floor.str(13))


if __name__=="__main__":
    verify()
