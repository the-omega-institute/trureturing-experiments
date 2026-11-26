"""冻结三壳 171 点后，把等权平方倒距离函数外推到 N−2、N+2。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx


CENTER = 2_306_888_100_085
OFFSETS = (-2,2)
WINDOWS = (
    ((".00000601",".00000629"), (".00306",".00316"), (".8573",".8659")),
    ((".0000133",".0000137"), (".000265",".000271"), (".01392",".01410")),
)


def response(energy, training):
    gaps = [energy-node for node,_ in training]
    assert all(not gap.contains(0) for gap in gaps)
    weights = [1/gap**2 for gap in gaps]
    assert all(w>0 for w in weights)
    value = sum(w*y for w,(_,y) in zip(weights,training,strict=True))/sum(weights)
    assert arb("91.1876")/arb("1.22089e19") < value < arb(".999")
    return value


def interlacing(raw):
    # 排序各端点，使真实重根或过近的特征值也得到外包。
    p,d = raw[0,0].real,raw[1,1].real
    rad = abs(acb((p-d)/2,abs(raw[0,1])))
    mid = (p+d)/2
    iso = raw[2,2].real
    helicity = raw[2,3].imag
    four = [mid-rad,mid+rad,iso+helicity,iso-helicity]
    lower = sorted(z.lower() for z in four)
    upper = sorted(z.upper() for z in four)
    return [(lower[j],upper[j+1]) for j in range(3)]


def report_quarks(root,n,training,frames,shell,ordered,native):
    height = acb.zeta_zero(n).imag
    projs = native["native"](root,n)
    adj = native["adj"]
    energy = run_path(str(root/"ch39-three-gram-source-flavor.py"))["projections"]
    states = [[],[]]
    for which,(q,rows) in enumerate(zip((-4,2),projs,strict=True)):
        weights = [response(t,training) for t in ordered(shell,q,height)]
        raw = sum((P*w for P,w in zip(rows,weights,strict=True)),acb_mat(4,4))
        bounds = interlacing(raw)
        fits = [not (hi<arb(a) or lo>arb(b)) for (lo,hi),(a,b) in zip(bounds,WINDOWS[which],strict=True)]
        print("任意复三平面",n,q,[(lo.str(12),hi.str(12)) for lo,hi in bounds],fits,flush=True)
        for label,frame in enumerate(frames):
            effective=adj(frame)*raw*frame
            levels,proj=energy(effective*effective)
            states[label].append(proj)
            print("固定空间 Yukawa",n,label,q,*[z.sqrt().str(12) for z in levels],flush=True)
    for label,(up,down) in enumerate(states):
        mix=[[(pu*pd).trace().real for pd in down] for pu in up]
        j=(up[0]*down[0]*up[1]*down[1]).trace().imag
        print("固定空间 CKM 模平方和有符号 J",n,label,[[v.str(12) for v in row] for row in mix],j.str(12),flush=True)


def report_gravity(root,n,training,field,moments,heat):
    pivots = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    weak = run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"]
    _,levels,c,_ = next(x for x in pivots["shell_weights"](root,n) if x[0]==6)
    coeffs=[z.real for z in c[0]]
    assert sum(coeffs).contains(1)
    y=[response(t,training) for t in levels]
    x=[z*z for z in y]
    mp=arb("1.22089e19")
    p0,p1=(arb("91.1876")/mp)**2,(arb("1.176397e14")/mp)**2
    assert all(z>p0 for z in x)
    columns=[(*moments(field,heat,z),weak(z,p0,p1)) for z in x]
    total=[sum(a*v for a,v in zip(coeffs,col,strict=True)) for col in zip(*columns,strict=True)]
    print("荷六截止范围与 K,W,D",n,min(z.lower() for z in y).str(12),
          max(z.upper() for z in y).str(12),*[z.str(13) for z in total],
          "含目标",[z.contains(a) for z,a in zip(total,(arb(1),arb("1e-122"),arb(0)),strict=True)],flush=True)


def verify(root):
    records,frames=run_path(str(root/"ch39-three-shell-two-space-readout.py"))["verify"](root,return_data=True)
    assert len(records)==171
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ordered=run_path(str(root/"ch39-bounded-common-rational-amplitude.py"))["ordered"]
    native=run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    gravity=run_path(str(root/"ch39-three-index-free-atoms.py"))
    numbers=[CENTER+offset for offset in OFFSETS]
    moments,heat,fields=gravity["broken_fields"](root,numbers)
    for n,field in zip(numbers,fields,strict=True):
        ctx.prec=320
        report_quarks(root,n,records,frames,shell,ordered,native)
        report_gravity(root,n,records,field,moments,heat)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
