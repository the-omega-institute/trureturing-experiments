"""Exclude both lower top orbits with joint current errors and neutral Yukawa running."""
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
fit = run_path(str(root/"ch39-finite-current-width-inverse.py"))
flavour = run_path(str(root/"ch39-neutral-yukawa-window.py"))
neutral = run_path(str(root/"ch39-neutral-current-threshold-inverse.py"))
study, q = fit["study"], fit["q"]
adjoint, comp = study["adjoint"], study["comp"]


def trajectory(index,upper,onset,h3,h8,casimir,time):
    scalar = (h3.upper(),h8.upper()) if upper else (h3.lower(),h8.lower())
    window = (onset.upper(),h3.lower()) if upper else (onset.lower(),h3.upper())
    breaks = [q(23),window[0],h3.lower(),h3.upper(),scalar[1],time]
    state = flavour["before"](index,upper)
    for a,b in zip(breaks,breaks[1:]):
        assert a < b
        mid = (a+b)/2
        flow = adjoint["flow"](mid>scalar[0],mid>scalar[1])
        if window[0]<mid<window[1]:
            flow = {**flow,"NEUTRAL_YUKAWA":casimir}
        n = int(float((b-a).upper())*12)+1
        for _ in range(n):
            state,tube = comp["advance"](state,(b-a)/n,upper,flow=flow,need_tube=True)
            assert all(abs(v)<2 for v in tube)
    return [1/(r*v) for r,v in zip(adjoint["ratios"],state[:3],strict=True)]


def verify():
    _,_,L,z,rho,eta,c,(triplet,octet) = fit["inferred_parameters"]()
    time = study["nominal"][2]
    cal = time+L.log()/2
    e,charges,weak,pair = neutral["first_window"](triplet,octet,c,z/rho,eta,1/rho)
    onset,h3,h8 = cal+e.log()/2,cal+triplet.log(),cal+octet.log()
    minimum,maximum,_ = flavour["charge_bounds"](*charges)
    assert weak/(h3-onset)<q(".001")
    exclude(onset,h3,h8,minimum,maximum,weak,pair)


def exclude(onset,h3,h8,minimum,maximum,weak,pair):
    time = study["nominal"][2]
    for index,offset in enumerate((".00000266",".00000212")):
        trial = time+q(offset)
        lo = trajectory(index,False,onset,h3,h8,maximum,trial)
        hi = trajectory(index,True,onset,h3,h8,minimum,trial)
        gap = study["target_at"](trial)-hi[1]+weak
        separation = lo[1]-hi[0]+pair
        assert gap > arb(2)**-21 and separation > arb(2)**-20
        for upper,bound in ((False,maximum),(True,minimum)):
            trajectory(index,upper,onset,h3,h8,bound,q(33))
        print("joint six-point and Yukawa exclusion",index,gap.lower().str(10),separation.lower().str(10))


if __name__ == "__main__":
    verify()
