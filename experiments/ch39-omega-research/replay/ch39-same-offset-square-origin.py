"""Evaluate the independent RH target and generate a finite witness for the second moment."""
import json
from pathlib import Path
from runpy import run_path
import sys
from flint import acb, arb, arb_mat, arb_series, ctx

root = Path(__file__).resolve().parent
model = run_path(str(root/"ch39-same-offset-square.py"))
base = model["base"]


def calibrate(case, paired=False):
    ctx.prec,ctx.cap = 1500,4
    def xi(j):
        t = arb_series([arb(1)/2+j,1])
        return t*(t-1)/2*((-t/2)*arb.pi().log()).exp()*(t/2).gamma()*t.zeta()
    zero = xi(0)
    mass = zero[2]/zero[0]
    J = model["J"]
    v = xi(J)
    mean = v[1]/(2*J*v[0]*mass)
    square = mean/2+(2*v[2]/v[0]-(v[1]/v[0])**2)/(4*mass)
    target = int((square*2**288).floor().unique_fmpz())
    ids = case["support"]
    prefix = base["data"]("ch39-tick-separators.json")["prefix"]+base["data"]("ch39-adjacent-offset-certificate.json")["tail"]
    extra = []
    if paired:
        ctx.prec = 544
        floor = run_path(str(root/"ch39-tick-separators-origin.py"))["floor_at"]
        heights = [v.imag for v in acb.zeta_zeros(len(prefix)+1,case["cut"]-len(prefix))]
        logs = [arb(z).log() for z in (2,3,arb(3)/2,6)]
        extra = [[floor(h,400),*[floor((h*z).cos(),72) for z in logs]] for h in heights]
        prefix += extra
    observe = lambda ids:model["square_system"](ids,prefix,target,paired)
    inverse = run_path(str(root/"ch39-cutoff-threshold-origin.py"))["build"]([ids],observe)[0]
    ctx.prec,ctx.cap = 1800,4
    z = [1/acb.zeta_zero(i).imag**2 for i in ids]
    offsets = [*base["OFFSETS"],1683] if paired else base["OFFSETS"]
    matrix = [[1/(1+j*j*t) for t in z] for j in offsets]+[[1/(1+J*J*t)**2 for t in z]]
    for row in observe(ids)[0][-4:]:
        matrix.append([(arb(l)/base["SCALE"]).union(arb(r)/base["SCALE"]) for l,r in row])
    ys = [arb(1)]
    for j in offsets[1:]:
        v = xi(j)
        ys.append(v[1]/(2*j*v[0]*mass))
    inv = arb_mat([[v.mid() for v in row] for row in matrix]).inv()
    trial = (inv*arb_mat([[v] for v in ys+[square]+[0]*4])).entries()
    w = [int((x.mid()*2**288).floor().unique_fmpz()) for x in trial[:-1]]
    w.append(2**288-sum(w))
    result = {"case":{**case,**inverse,"weights":w},"square_floor":target}
    if paired:
        result["tail"]=extra
    return result


if __name__=="__main__":
    paired = "--paired" in sys.argv
    args = [v for v in sys.argv[1:] if v!="--paired"]
    path = root/("ch39-paired-products-data.json" if paired else "ch39-same-offset-square-data.json")
    if args:
        path.write_text(json.dumps(calibrate(json.loads(Path(args[0]).read_text()),paired))+"\n")
    else:
        record = json.loads(path.read_text())
        assert record==calibrate(record["case"],paired)
    print("RH squared-kernel target and finite witnesses calibrated.")
