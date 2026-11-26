"""在轻模分界 9 处用共同有理函数抑制重费米子顶角，重算质量、CKM 与质量四次迹。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


root = Path(__file__).parent
products = run_path(str(root / "ch39-minimal-product-transport.py"))
model = products["study"]
adjoint = model["adjoint"]
shell = model["shell"]
sectors = ((1,-4), (1,2), (-3,6))


def filter_value(level, power):
    return 1/(1 + (level/9)**power)


def spectrum(lap, center, q):
    eigenvalues, vectors = lap.eig(right=True, algorithm="rump")
    assert all(not z.overlaps(prev) for j,z in enumerate(eigenvalues) for prev in eigenvalues[:j])
    levels = [abs(z.real*(z.real-center)) for z in eigenvalues]
    order = sorted(range(19),key=lambda i: float(levels[i]))
    levels = [levels[i] for i in order]
    assert levels[2]<9<levels[3] and all(a<b for a,b in zip(levels, levels[1:]))
    basis = acb_mat(19,19)
    for j, index in enumerate(order):
        norm = sum(abs(vectors[s,index])**2 for s in range(19)).sqrt()
        for s in range(19):
            basis[s,j] = vectors[s,index]/norm
    transport = acb_mat([[-lap[i,j]/8 if i!=j else 0 for j in range(19)] for i in range(19)])
    unit = arb("2.435e18")/724**2 * (32 if q%2==0 else 1)
    masses = [arb(0)]*3 + [unit*d for d in levels[3:]]
    return levels, masses, (basis, transport*basis)


def source(number, weights, higgs):
    height = acb.zeta_zero(number).imag
    neutral = shell["laplacian"](0,height)
    center = sum(neutral[i,i].real for i in range(19))/18
    parts = {q:spectrum(shell["laplacian"](q,height),center,q) for q in (1,-4,2,-3,6)}
    h_modes = model["edge"]["mode_columns"](model["edge"]["filters"]["spectra"](3,height,center))
    for j in range(3):
        anchor = next(h_modes[i,j] for i in range(19) if abs(h_modes[i,j])>arb(".01"))
        for i in range(19):
            h_modes[i,j] *= anchor.conjugate()/abs(anchor)
    h = h_modes*higgs/(adjoint(higgs)*higgs)[0,0].real.sqrt()
    lap = shell["laplacian"](3,height)
    transported = acb_mat([[-lap[i,j]/8 if i!=j else 0 for j in range(19)] for i in range(19)])*h
    matrices = []
    for left,right in sectors:
        values = (h,transported) if right==-4 else (h.conjugate(),transported.conjugate())
        y = acb_mat(19,19)
        for ix,(a,b,c) in enumerate(model["masks"]):
            r = parts[right][2][b]
            hr = acb_mat([[values[c][i,0]*r[i,j] for j in range(19)] for i in range(19)])
            y += weights[ix,0]*(parts[left][2][a].transpose()*hr)
        matrices.append(y)
    return parts, matrices


def observables(parts, matrices, values):
    light, heavy_norms = [], []
    heat = arb(0)
    for (left,right), y, colors in zip(sectors,matrices,(3,3,1)):
        ml, mr = parts[left][1], parts[right][1]
        fl, fr = values[left], values[right]
        out = acb_mat([[fl[i]*y[i,j]*fr[j] for j in range(19)] for i in range(19)])
        light.append(acb_mat([[out[i,j] for j in range(3)] for i in range(3)]))
        heavy_entries = [(i,j,abs(out[i,j])**2) for i in range(19) for j in range(19)
                         if i>=3 or j>=3]
        heavy_norms.append(sum(x for i,j,x in heavy_entries).sqrt())
        heat += colors*sum((ml[i]**2+mr[j]**2)*x for i,j,x in heavy_entries)
    return light, heavy_norms, heat


def check():
    refs = run_path(str(root / "ch39-shell-flavor-matching-windows.py"))
    mixing = run_path(str(root / "ch39-finite-orthogonal-ckm.py"))
    ctx.prec = 448
    inputs = ((22027, model["WITNESS"], model["HIGGS"]),
              (28325, products["SECOND_WITNESS"], products["SECOND_HIGGS"]))
    for number in (11661,11934,12207,22027,28325):
        height = acb.zeta_zero(number).imag
        neutral = shell["laplacian"](0,height)
        center = sum(neutral[i,i].real for i in range(19))/18
        for q in (1,-4,2,-3,6):
            roots = shell["laplacian"](q,height).eig(algorithm="rump")
            levels = [abs(r.real*(r.real-center)) for r in roots]
            low = [d for d in levels if d < 9]
            high = [d for d in levels if d > 9]
            assert len(low)==3 and len(high)==16
            assert all(d < arb("7.66") for d in low) and all(d > arb("10.56") for d in high)
    assert filter_value(arb("7.66"),256)**2 > arb(".987")
    for number,w,z in inputs:
        parts, matrices = source(number,model["column"](w),model["column"](z))
        entries = []
        for (left,right), y, colors in zip(sectors,matrices,(3,3,1)):
            ml,mr = parts[left][1],parts[right][1]
            entries.append((left,right,[(i,j,colors*(ml[i]**2+mr[j]**2)*abs(y[i,j])**2)
                                         for i in range(19) for j in range(19) if i>=3 or j>=3]))
        for power in range(2,153,2):
            sq = {q:[filter_value(d,power)**2 for d in p[0]] for q,p in parts.items()}
            heat = sum(sq[q][i]*sq[r][j]*x for q,r,rows in entries for i,j,x in rows)
            if power < 152:
                assert heat > arb("125.25")**2
            else:
                assert heat < arb("125.25")**2
        for power in (0,150,152,256):
            values = {q:[filter_value(d,power) if power else arb(1) for d in p[0]]
                      for q,p in parts.items()}
            light, norms, heat = observables(parts,matrices,values)
            masses = [y for mat in light for y in sorted(
                (v.real.sqrt() for v in (mat*adjoint(mat)).eig(algorithm="rump")),key=float)]
            target = [y for g in ("up","down","leptons") for y in refs[g]]
            error = max(abs(y/ref-1).upper() for y,ref in zip(masses,target))
            angles = mixing["observables"](*light[:2])[4:8]
            goals = [arb(t).sin() for t in (".22704",".04275",".00377")] + [arb("1.139")]
            angle_error = max(abs(x-ref).upper() for x,ref in zip(angles,goals))
            if power == 256:
                assert heat < arb("4.22e-11")
                assert error < arb("3.3e-8") and angle_error < arb("1e-8")
            print(number, power, "质量误差、混合误差、最大含重模范数、重质量平方矩 GeV²",
                  *(v.str(8) for v in (error,angle_error,max(n.upper() for n in norms),heat)))


if __name__ == "__main__":
    check()
