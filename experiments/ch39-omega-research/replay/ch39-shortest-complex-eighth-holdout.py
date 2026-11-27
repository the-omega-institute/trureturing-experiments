"""保持第七壳认证复门，对同一个一次共轭的正算子作第八壳引力留出。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


NUMBER = 2_306_888_100_093


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-"+name+".py")))

    x,y = lib("shortest-common-complex-gate")["experiment"](root,fifth,columns,levels)
    ctx.prec = 256
    denominator = 1+x*x+y*y
    c = 1-2/denominator
    z = 2*acb(x,y)/denominator
    station = fifth*columns
    moved = acb_mat(station)
    shell = lib("three-ray-chiral-boundary")
    _, pairs = lib("behrend-rational-chord-pairs")["inspect"](shell["geo"]["sites"],867)
    assert pairs == [(6,11,12),(7,12,12),(8,18,12)]
    for a,b,_ in pairs:
        for j in range(2):
            moved[a,j] = c*station[a,j]-z.conjugate()*station[b,j]
            moved[b,j] = z*station[a,j]+c*station[b,j]
    values,basis = lib("adjacent-pvm-transport")["spectrum"](shell,6,NUMBER)
    _,nodes,pivots,_ = next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root,NUMBER) if row[0]==6)
    assert all(a.overlaps(b) for a,b in zip(values,nodes,strict=True))
    weights = [v.real for v in pivots[2]]
    assert sum(weights).contains(1)
    _,heat,(field,) = lib("three-index-free-atoms")["broken_fields"](root,[NUMBER])
    overlap = basis.conjugate().transpose()*moved
    low,heavy,high = levels
    amplitudes = [low+(heavy-low)*abs(overlap[i,0])**2+(high-low)*abs(overlap[i,1])**2 for i in range(19)]
    assert (sum(amplitudes)-(17*low+heavy+high)).contains(0)
    newton,vacuum = arb(0),arb(0)
    for w,a in zip(weights,amplitudes,strict=True):
        l = max(low.lower(),a.lower())**2
        r = min(high.upper(),a.upper())**2
        assert 0<l<=r
        midpoint = (l.union(r)).mid()
        ak = midpoint.sqrt()
        moment_k,moment_w,der_k,der_w = [arb(0) for _ in range(4)]
        for g,alpha,mass in field:
            moment_k += alpha*midpoint*heat(2,mass/midpoint)/(2*arb.pi())
            moment_w -= g*midpoint**2*heat(3,mass/midpoint)/(4*arb.pi())
            kernel = (-mass/l).exp().union((-mass/r).exp())
            der_k += alpha*kernel
            der_w -= g*kernel
        interval = l.sqrt().union(r.sqrt())
        shift = interval-ak
        newton += w*(moment_k+shift*interval*der_k/arb.pi())
        vacuum += w*(moment_w+shift*interval**3*der_w/(2*arb.pi()))
    print("第八壳完整根盒的牛顿项、真空项",newton.str(16),vacuum.str(16),flush=True)
    assert arb("-6.379")<newton<arb("-6.361")
    assert arb("1.267")<vacuum<arb("1.273")
    return newton,vacuum


def verify(root):
    op = run_path(str(root/"ch39-fifth-full-positive-joint.py"))["verify"](root,return_operator="spectral")
    return experiment(root,*op)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
