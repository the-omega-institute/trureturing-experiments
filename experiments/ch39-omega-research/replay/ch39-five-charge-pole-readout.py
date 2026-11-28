"""十六站盲点随软化参数的延拓及补回轨道后的局部秩。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


# Each chart fixes one nonzero eigenvector coordinate per charge.
CHARTS = (
    ("3.8031598886818507", "2.57243365038201",   (1,1,0,0,1)),
    ("2.480025418497736",  "3.710751656797576",  (0,0,1,1,0)),
    ("1.6980941908310234", "4.648816431388861",  (0,0,1,0,1)),
    ("2.4405211865612157", "3.7111059405035363", (0,0,1,1,0)),
    ("3.8426641206183736", "2.5720793666760513", (1,1,0,0,1)),
    ("2.4845918065736123", "3.6950357645185226", (0,0,1,1,0)),
)


def verify(root):
    ctx.prec = 200
    common = run_path(str(Path(__file__).with_name("ch39-common-charge-amplitudes.py")))
    pump = run_path(str(root / "ch39-native-light-projector-pump.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    charges = common["CHARGES"]
    adj = pump["adjoint"]
    first_line = pump["first_line"]
    sites = shell["geo"]["sites"]
    rotation = [sites.index((b,c,a)) for a,b,c in sites]
    orbits = sorted({tuple(sorted({j,rotation[j],rotation[rotation[j]]})) for j in range(19)})
    pole = shell["soft"]["pole"]
    assert [o for o in orbits if len(o) == 1] == [(pole,)]
    triples = [o for o in orbits if len(o) == 3]
    soft = acb_mat(19,19)
    for site in shell["soft"]["star"]["orbit"]:
        soft[pole,pole] -= 1
        soft[site,site] -= 1
        soft[pole,site] += 1
        soft[site,pole] += 1

    blocks = []
    for frame in pump["rotation_bases"](shell):
        rows = [next(j for j in range(frame.ncols()) if abs(frame[o[0],j]) > 0) for o in triples]
        operators = []
        for charge in charges:
            base, links = pump["clock_matrices"](shell,charge)
            operators.append(tuple(adj(frame)*m*frame for m in (base,*links)))
        blocks.append((rows,operators,adj(frame)*soft*frame))
        if len(blocks) == 2:
            assert all(frame[pole,j].is_zero() for j in range(frame.ncols()))

    def matrices(angles, deficit, block):
        _, operators, perturbation = block
        for q,(base,link1,link2) in zip(charges,operators,strict=True):
            factor = acb(0,q)
            a = link1*(factor*angles[0]).exp()
            b = link2*(factor*angles[1]).exp()
            matrix = base+a+adj(a)+b+adj(b)+deficit*perturbation
            yield matrix, ((a-adj(a))*factor,(b-adj(b))*factor,perturbation)

    def system(angles, deficit, selected, pivots):
        lines = [common["eigenline_derivatives"](matrix,slopes,first_line,pivot)
                 for (matrix,slopes),pivot in zip(matrices(angles,deficit,blocks[1]),pivots,strict=True)]
        value, (dx,dy,ds) = common["determinant_derivatives"](lines,selected)
        return (arb_mat([[value.real],[value.imag]]),
                arb_mat([[dx.real,dy.real],[dx.imag,dy.imag]]), arb_mat([[ds.real],[ds.imag]]))

    radius = arb("1e-8")
    deficit = arb(0).union(radius)
    eye = arb_mat([[1,0],[0,1]])
    for omitted, (x,y,pivots) in enumerate(CHARTS):
        center = [arb(x),arb(y)]
        box = [c+arb(0,radius.upper()) for c in center]
        selected = [j for i,j in enumerate(blocks[1][0]) if i != omitted]
        residual, jacobian, _ = system(center,arb(0),selected,pivots)
        pre = jacobian.mid().inv()
        assert not pre.det().contains(0)
        _, whole, drift = system(box,deficit,selected,pivots)
        move, rate = pre*residual, pre*drift
        defect = eye-pre*whole
        ratio = max(move[i,0].abs_upper()/radius + rate[i,0].abs_upper()
                    + sum(z.abs_upper() for z in row) for i,row in enumerate(defect.tolist()))
        assert ratio < arb(".05")
        floors = []
        for block in blocks:
            rows = block[0]
            lines = [first_line(matrix)[0] for matrix,_ in matrices(box,deficit,block)]
            lifted = acb_mat([[line[j] for line in lines] for j in rows])
            floors.append(common["left_inverse_bound"](lifted))
        floor = min(z.lower() for z in floors)
        assert floor > arb(".00467")
        print(triples[omitted], "参数盒自映比及十八站奇异值下界", ratio.str(11), floor.str(11))


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve())
