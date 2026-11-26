"""三边界上的共同复 Yukawa 系数、联合秩与物理质量分裂的幅度下界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

light = run_path(str(Path(__file__).with_name("ch39-boundary-light-yukawa.py")))
model = light["boundary_model"]
adjoint = light["adjoint"]


def canonical(profile):
    root = acb(0,2*arb.pi()/3).exp()
    waves = profile*acb_mat([[root**(i*j) for j in range(3)] for i in range(3)])
    for k in range(3):
        norm = sum(abs(waves[i,k])**2 for i in range(19)).sqrt()
        for i in range(19):
            waves[i,k] /= norm
    return waves


def split_bound(q, up, down, higgs):
    # Align each right-handed Fourier channel by a unit-modulus phase.
    diff = acb_mat(19,3)
    for k in range(3):
        a = [higgs[i,0]*up[i,k] for i in range(19)]
        b = [higgs[i,0].conjugate()*down[i,k] for i in range(19)]
        pair = sum(x.conjugate()*y for x,y in zip(b,a))
        align = pair/abs(pair)
        for i in range(19):
            diff[i,k] = a[i]-align*b[i]
    # Split the diagonal coefficient using a positive site weight and bound
    # both rectangular operator norms by their three-by-three Gram matrices.
    left, right = acb_mat(19,3), acb_mat(19,3)
    for i in range(19):
        weight = (sum(abs(diff[i,k])**2 for k in range(3))/sum(abs(q[i,k])**2 for k in range(3)))**(arb(1)/4)
        for k in range(3):
            left[i,k] = q[i,k]*weight
            right[i,k] = diff[i,k]/weight

    def norm_upper(matrix):
        gram = adjoint(matrix)*matrix
        return max(sum(abs(gram[i,j]).upper() for j in range(3)) for i in range(3)).sqrt()

    return norm_upper(left)*norm_upper(right)


def local_map(left, right, scalar):
    return [[left[i,j]*scalar[i,0]*right[i,k] for i in range(19)]
            for j in range(3) for k in range(3)]


def phase_spread(q, up, down, higgs):
    # One third of a trace with a fixed right-unitary permutation has
    # modulus at most the largest down singular value. The numerator
    # gives its individual site contributions in the up rank-one norm.
    entries = []
    for i in range(19):
        norm = abs(higgs[i,0])*(sum(abs(q[i,k])**2 for k in range(3))
                                *sum(abs(up[i,k])**2 for k in range(3))).sqrt()
        tr = sum(q[i,k]*down[i,(-k)%3] for k in range(3))*higgs[i,0].conjugate()
        entries.append(tr/(3*norm))
    positive = min(z.real.lower() for z in entries)
    ratio = (arb(".01401")+arb(".00009"))/(arb(".8616")-arb(".0043"))
    def lower_bound(degrees):
        angle = degrees*arb.pi()/180
        return min((z.real*angle.cos()-abs(z.imag)*angle.sin()).lower() for z in entries)
    lo,hi = arb(0),arb(90)
    for _ in range(32):
        mid = (lo+hi)/2
        if lower_bound(mid)>ratio:
            lo = mid
        else:
            hi = mid
    assert positive > arb(".32") and lower_bound(lo)>ratio
    print("同相站点的底/顶比下界、容许目标层级所需相位弧长下界（度）",positive.str(10),(2*lo).str(10))


def verify():
    ctx.prec = 256
    height = acb.zeta_zero(11431).imag
    boundary = model["boundary_sites"]()
    pole = model["soft"]["pole"]
    p = {q:light["extension"](q,height,boundary) for q in (1,-4,2,-3,6)}
    h = light["extension"](3,height,[pole])
    quark = local_map(p[1],p[-4],h)+local_map(p[1],p[2],h.conjugate())
    keep = [i for i in range(19) if i!=pole]
    square = acb_mat([[row[i] for i in keep] for row in quark])
    determinant = abs(square.det())
    assert arb("1.82619e-39") < determinant < arb("1.82620e-39")
    kernel = acb_mat(19,1)
    kernel[pole,0] = 1
    solved = square.solve(acb_mat([[-row[pole]] for row in quark]))
    for k,i in enumerate(keep):
        kernel[i,0] = solved[k,0]
    assert all(z.contains(0) for row in (acb_mat(quark)*kernel).tolist() for z in row)
    lepton = acb_mat(local_map(p[-3],p[6],h.conjugate()))*kernel
    lepton_change = acb_mat([[lepton[3*i+j,0] for j in range(3)] for i in range(3)])
    assert abs(lepton_change.det()) > 0
    print("十八个非极点系数的夸克子式及极点取一时的轻子变化行列式",
          determinant.str(13), abs(lepton_change.det()).str(13))

    h = h/sum(abs(h[i,0])**2 for i in range(19)).sqrt()
    q, up, down = [canonical(p[charge]) for charge in (1,-4,2)]
    coefficient = split_bound(q,up,down,h)
    assert coefficient < arb(".037376")
    # The conservative top/bottom edge of the separately quoted 1 TeV intervals.
    target_gap = arb(".8616")-arb(".0043")-arb(".01401")-arb(".00009")
    assert target_gap/coefficient > arb("22.55")
    print("共同站点系数对最大奇异值差的界，以及 1 TeV 所需的最大系数下界",
          coefficient.str(12), (target_gap/coefficient).str(12))
    phase_spread(q,up,down,h)


if __name__ == "__main__":
    verify()
