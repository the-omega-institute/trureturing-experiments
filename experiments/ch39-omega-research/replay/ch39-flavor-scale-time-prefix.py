"""在参考味行列式唯一尺度上找到算术时间前缀的严格分界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx


root=Path(__file__).parent
determinant=run_path(str(root/"ch39-rh-flavor-determinant.py"))
clock=run_path(str(root/"ch39-three-mode-time-localization.py"))


def certify_weights(times,scale,degree,heights):
    assert len(heights)==degree+1
    moments=clock["clock"]["cleared_moments"](times,degree)
    center=scale.mid()
    nodes=[200*(center/h)**2 for h in heights]
    ratio=(center/scale)**2
    probabilities=[]
    for i,node in enumerate(nodes):
        cardinal=arb_poly([1])
        denominator=arb(1)
        for j,other in enumerate(nodes):
            if i!=j:
                cardinal*=arb_poly([-other,1])
                denominator*=node-other
        # 倾斜矩乘以共同尺度的幂。先换到 ratio=1 附近以保留抵消。
        weight_poly=arb_poly([cardinal[j]*moments[j] for j in range(degree+1)])
        weight=weight_poly(arb_poly([1,1]))(ratio-1)/denominator
        assert weight>0
        actual=node/ratio
        for j in range(1,degree+1):
            weight*=1+j*j*actual/200
        probabilities.append(weight)
    assert sum(probabilities).contains(1)
    for j in range(1,degree+1):
        value=sum(p/(1+j*j*x/(200*ratio)) for x,p in zip(nodes,probabilities,strict=True))
        assert value.overlaps(times[j])
    return probabilities


def separation_margin(times,scale,degree,extra_power,lower_height,upper_height):
    left=(200*(scale/upper_height)**2).upper()
    right=(200*(scale/lower_height)**2).lower()
    assert left<right
    nonnegative=arb_poly([0,1])**extra_power*arb_poly([-left,1])*arb_poly([-right,1])
    moments=clock["clock"]["cleared_moments"](times,degree)
    order=(degree-nonnegative.degree())//2
    tilted=[sum(nonnegative[j]*moments[i+j] for j in range(nonnegative.degree()+1))
            for i in range(2*order+1)]
    gram=arb_mat([[tilted[i+j] for i in range(order)] for j in range(order)])
    coefficients=gram.solve(arb_mat([[-tilted[i+order]] for i in range(order)]))
    polynomial=arb_poly([coefficients[i,0].mid() for i in range(order)]+[1])
    numerator=-nonnegative*polynomial**2
    nominal=sum(numerator[i]*moments[i] for i in range(numerator.degree()+1))
    residues=clock["residues_of"](numerator,clock["partial_fractions"](degree))
    return nominal/sum(abs(c) for c in residues)


def verify():
    ctx.prec,ctx.cap=1536,3
    times=clock["clock"]["times"]()
    low,high=determinant["scale_interval"]()
    scale=low.union(high)
    first_indices=(2556808,3935900,4752291,5932061,6885101,8724973,12448034,21829052,78456164)
    ctx.prec=384
    heights=[acb.zeta_zero(1).imag]
    for index in first_indices:
        heights.extend(z.imag for z in acb.zeta_zeros(index,2))
    lower_height,upper_height=heights[1:3]
    assert lower_height<scale*heights[0]<upper_height
    ctx.prec=1536
    probabilities=certify_weights(times,scale,18,heights)
    margin=separation_margin(times,scale,19,1,lower_height,upper_height)
    assert margin>arb("6.62e-30")
    print("十九枚的最小原概率下界",min(v.lower() for v in probabilities).str(12))
    print("第十九项前缀的逐项最大误差下界",margin.str(13))


if __name__=="__main__":
    verify()
