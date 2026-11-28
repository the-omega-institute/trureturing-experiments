"""实际两部门各前进一枚零点时，冻结四个协方差并比较重新校准的目标。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


PAIRS = ((2345666, 964890), (2345667, 964891))
CENTERS = (
    (("1.657055709135", "4.507081302329"), ("2.130442742842", ".906317020444"),
     ("5.058382425325", ".640595752598"), ("5.450512188429", "4.752972164129")),
    (("1.755250855303", "3.447397832320"), ("4.600967189839", "2.335319400971")),
)
UP_SCALE = "7.06539940386823251"


def certify_centers(equations, trig, centers):
    radius = arb(".02")
    pad = arb(0, radius.upper())
    points = [(arb(x), arb(y)) for x, y in centers]
    assert all(radius < t < 2*arb.pi()-radius for point in points for t in point)
    assert all(abs(x-a)>2*radius or abs(y-b)>2*radius for (x,y),(a,b) in combinations(points,2))
    identity = arb_mat([[1,0],[0,1]])
    for x,y in points:
        values, derivative = trig(x,y)
        jac = equations*derivative
        pre = jac.mid().inv()
        shift = pre*equations*values
        error = identity-pre*equations*trig(x+pad,y+pad)[1]
        assert all(abs(shift[i,0])/radius + sum(abs(error[i,j]) for j in range(2)) < 1
                   for i in range(2))
    cells = [(0,0,0)]
    leaves = 0
    while cells:
        depth,i,j = cells.pop()
        half = arb.pi()/2**depth
        pad = arb(0,half.upper())
        x,y = (2*i+1)*half+pad, (2*j+1)*half+pad
        inside = any(abs(x-a)<radius and abs(y-b)<radius for a,b in points)
        if inside or any(not v.contains(0) for v in (equations*trig(x,y)[0]).entries()):
            leaves += 1
        else:
            assert depth < 20
            cells.extend((depth+1,2*i+a,2*j+b) for a,b in product((0,1), repeat=2))
    return points,leaves


def shell_axes(model):
    source,adj = model["sources"],model["adj"]
    sites = model["shell"]["geo"]["sites"]
    cycle = [sites.index((y,z,x)) for x,y,z in sites]
    rotation = acb_mat([[int(i==j) for j in cycle] for i in range(len(sites))])
    h = acb(0,1)*adj(source)*(rotation-adj(rotation))*source
    pole = model["shell"]["soft"]["pole"]
    v = adj(acb_mat([[source[pole,j] for j in range(3)]]))
    p = v*adj(v)
    c = acb(0,1)*(h*p-p*h)
    krylov = (v,h*v,h*h*v)
    basis = acb_mat([[w[i,0] for w in krylov] for i in range(3)])
    assert abs(basis.det()) > arb("3.08e-7")
    axes = []
    for i,j in combinations(range(3),2):
        outer = krylov[j]*adj(krylov[i])
        axes.append((outer-adj(outer))/acb(0,2))
    weights = [(axis*c).trace().real for axis in axes]
    assert abs(weights[0]) > arb("3e-7")
    return [axes[i]-(weights[i]/weights[0])*axes[0] for i in (1,2)], c


def readings(model, up, down, change):
    pair = [model["eig"](model["adj"](a)*a) for a in (up*change*arb(UP_SCALE),down*change)]
    return ([s.sqrt() for levels,_ in pair for s in levels]
            + model["mixing_angles"](pair[0][1],pair[1][1]))


def verify(replay):
    def load(name):
        return run_path(str(replay/("ch39-"+name+".py")))
    orbit = load("two-zero-orbit-condition")
    model,_,_ = orbit["calibration"]()
    gram = load("quadratic-gram-source")
    words = load("two-zero-quadratic-source")
    trig = words["trigonometric_columns"]
    refine = load("four-quadratic-source-neighbours")["refine"]
    load("moving-preparation-quark-shape")["target_gaps"]()
    ctx.prec = 224
    normals, commutator = shell_axes(model)
    adj = model["adj"]
    frozen = []
    for index,(up_zero,down_zero) in enumerate(PAIRS):
        for n in (up_zero,down_zero):
            height = acb.zeta_zero(n).imag
            left,right = height.mid()-arb("1e-12"),height.mid()+arb("1e-12")
            assert left<height<right
            assert [t.zeta_nzeros().unique_fmpz() for t in (left,right)] == [n-1,n]
            assert model["hardy"](left)*model["hardy"](right)<0
        up = gram["at_zero"](model,-4,up_zero)
        down = gram["at_zero"](model,2,down_zero)
        a,b,c = model["singular"](up*down.inv())
        assert 0<a<b<c
        reference = readings(model,up,down,model["K"])
        eigenpairs = [model["eig"](adj(y)*y) for y in (up*model["K"],down*model["K"])]
        assert all(0 < vals[0] < vals[1] < vals[2] for vals,_ in eigenpairs)
        mixing = adj(eigenpairs[0][1])*eigenpairs[1][1]
        assert all(abs(z)>0 for z in mixing.entries())
        if index == 0:
            assert all(abs(v-m)<r for v,m,r in zip(reference,model["ref"],model["tol"],strict=True))
        else:
            assert arb("3.06771") < b/a < arb("3.06772")
            assert arb("1.72654") < c/b < arb("1.72655")
            assert max(words["norm"](x.inv()) for x in (up,down)) < 54
            for m,lo,hi in zip(frozen,("41.69","9.92","42.05","15.00"),
                               ("41.70","9.93","42.06","15.01"),strict=True):
                vals = readings(model,up,down,m)
                errors = [abs(v-r)/t for v,r,t in zip(vals,model["ref"],model["tol"],strict=True)]
                assert all(z>arb(lo) for z in errors) and any(z<arb(hi) for z in errors)
                print("冻结源后十项中最近的偏差，以参考半宽计",lo,hi)
        parts,_ = orbit["orbit"](model,up,down)
        equations = words["orbit_coefficients"](parts,adj,normals)
        centers,leaves = certify_centers(equations,trig,CENTERS[index])
        signs = []
        for x,y in centers:
            x,y = refine((x,y,0),arb(".02"),equations,trig)
            m = orbit["phase_sum"](parts,orbit["phases_at"](x,y))
            s = m*adj(m)
            ell = (commutator*s).trace().real
            assert abs(m.det())>0 and not ell.contains(0)
            assert all(d.overlaps(r) for d,r in zip(readings(model,up,down,m),reference,strict=True))
            signs.append(1 if ell>0 else -1)
            if index==0:
                frozen.append(m)
        assert signs == ([-1,1,-1,1],[1,-1])[index]
        print("上下零点编号、覆盖终格、各根交换子迹符号",up_zero,down_zero,leaves,signs)


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
