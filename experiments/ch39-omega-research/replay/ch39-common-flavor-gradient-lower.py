"""跨高度固定系数与相位时，参考味双钟驻点的四分量梯度分离。"""
import json
from pathlib import Path
from runpy import run_path
import sys
from flint import acb, arb, ctx


def verify(root):
    model = run_path(str(root / "ch39-mincut-three-window-angle.py"))
    signs = run_path(str(root / "ch39-twentyseven-coefficient-window.py"))["signs"]
    inertia = run_path(str(root / "ch39-neutral-mean-phase-inertia.py"))
    winding = run_path(str(root / "ch39-mincut-global-winding.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    flavor = run_path(str(root / "ch39-chiral-neutral-flavor-phase.py"))
    ctx.prec = 192
    faces = winding["terms"](shell)
    reference = shell["geo"]["phases"](acb.zeta_zero(11661).imag)
    d = (reference[724]/reference[804]).arg()
    a,b = [flavor["singlet_at"]({724:z,804:acb(1)})[1] for z in (acb(1),acb(0,1))]
    assert a>0 and abs(b)<arb("1e-45")
    distances = {1:d, -1:d-2*(b/a).atan()}
    levels = (arb("18.94088212"), arb("18.94414547"))
    band = model["sign"]["density"]["eigen_band"]()
    assert levels[0]<band[0]<band[1]<levels[1]
    groups = json.loads(Path(__file__).with_name("ch39-common-flavor-gradient-separators.json").read_text())
    assert [(r["sign"],r["order"]) for r in groups] == [(1,56),(1,71),(-1,78),(-1,92),(-1,99)]
    keys = ((1,0),(0,1),(1,-1))
    co = []
    for q in range(1,100):
        bare = winding["coeff"](faces,q,arb(0))
        soft = [sum((z**q for z,n in faces[k] if n==1), acb(0)) for k in keys]
        co.append((bare,soft))

    def separated(t, group):
        s,order = group["sign"],group["order"]
        x,y = s*t,s*(t-distances[s])
        columns = []
        for q,(bare,soft) in enumerate(co[:order],1):
            phases = [acb(0,q*z).exp() for z in (x,y,x-y)]
            column = []
            for elements in (bare,soft):
                f,g,h = [(z*phase).imag for z,phase in zip(elements,phases,strict=True)]
                column.extend((f+h,g-h))
            columns.append(column)
        vectors = group["directions"]
        position = int(round(float(t)*1e9))
        before = [row for row in vectors if row[0]<=position]
        after = [row for row in vectors if row[0]>position]
        if before and after:
            left, u = before[-1]
            right, v = after[0]
            vector = [(right-position)*a+(position-left)*b for a,b in zip(u,v,strict=True)]
            magnitude = right-left
        else:
            vector = (before[-1] if before else after[0])[1]
            magnitude = 1
        return all(sum(v*z for v,z in zip(vector,column,strict=True)) > magnitude*arb(".001")
                   for column in columns)

    cut = arb(".075")
    for group in groups[:4]:
        s,order = group["sign"],group["order"]
        endpoints = (-arb.pi()/16,cut) if order in (56,78) else (cut,arb.pi()/16)
        pending = [endpoints]
        screened = certified = 0
        while pending:
            left,right = pending.pop()
            t = left.union(right)
            x,y = s*t,s*(t-distances[s])
            if (abs(y)>arb.pi()/16 or signs(inertia,2,levels[0],x,y)[0]>=16
                    or signs(inertia,2,levels[1],x,y)[1]<16):
                screened += 1
            elif separated(t,group):
                certified += 1
            else:
                assert right-left>arb("1e-9"),(s,order,t)
                midpoint=(left+right)/2
                pending.extend(((left,midpoint),(midpoint,right)))
        print("角距支、排除至绕数、曲率排除段与正投影段",s,order,screened,certified)
    assert separated(arb(".01046"),groups[-1])
    print("原三十六圈全局见证同一相位跨高度固定时，99 圈仍不能共同驻定")


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
