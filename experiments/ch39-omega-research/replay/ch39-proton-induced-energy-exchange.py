"""八站诱导子图能量条件对全壳传播与子图传播的不同作用。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx

SITES = (0,2,3,5,6,9,11,13)
LOW_SUPPORT = {120:(3,6,9),163:(3,9,11),177:(3,5,9)}


def receiver(propagator, support):
    row = SITES.index(5)
    x,y,z = [propagator[row,SITES.index(i)] for i in support]
    areas = [(u.conjugate()*v).imag for u,v in ((y,z),(z,x),(x,y))]
    assert all(a>0 for a in areas) or all(a<0 for a in areas)
    weights = {site:a/sum(areas) for site,a in zip(support,areas,strict=True)}
    return [weights.get(i,arb(0)) for i in SITES]


def energy_form(matrix, left, right):
    return sum(a*matrix[i,j].real*b for i,a in enumerate(left) for j,b in enumerate(right))


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=192
    edges=[(i,j) for i,j,_,_ in shell["geo"]["edges"] if i in SITES and j in SITES]
    assert len(edges)==19 and all(shell["soft"]["pole"] not in e for e in edges)
    degrees=[sum(i in e for e in edges) for i in SITES]
    assert degrees==[5,5,2,2,7,6,5,6]
    for index,low_support in LOW_SUPPORT.items():
        full = shell["laplacian"](5,acb.zeta_zero(index).imag)
        green = full.inv()
        induced = acb_mat([[degrees[r] if i==j else full[i,j] for j in SITES]
                           for r,i in enumerate(SITES)])
        assert all(abs(induced[i,j].real)<=1 for i in range(8) for j in range(8) if i!=j)
        assert full[6,6].real==18
        assert all(full[j,j].real <= (10 if j in (3,5) else 18) for j in SITES)
        for i in SITES:
            for j in SITES:
                lower = "-.00085" if {i,j} in ({3,6},{5,6}) else ".006"
                assert green[i,j].real>arb(lower)
        small_amplitude = arb(110)/137
        assert (small_amplitude+(2*small_amplitude-1).sqrt())/2 > arb(7)/9
        cancellation = arb(120)/137
        dominant = (cancellation+(2*cancellation-1).sqrt())/2
        assert 8*dominant**2-3>arb("3.07603")
        assert 19*dominant**2-11>arb("3.43057")

        local_green = induced.inv()
        low = receiver(local_green,low_support)
        high = receiver(local_green,(5,9,13))
        energy_low = energy_form(induced,low,low)
        energy_high = energy_form(induced,high,high)
        assert energy_low<2<energy_high
        direction=[v-u for u,v in zip(low,high,strict=True)]
        linear=2*energy_form(induced,low,direction)
        quadratic=energy_form(induced,direction,direction)
        assert quadratic>0
        t=(-linear+(linear**2+4*quadratic*(2-energy_low)).sqrt())/(2*quadratic)
        assert 0<t<1
        print(index,"接收支撑",low_support,(5,9,13),"能量",energy_low.str(12),energy_high.str(12),
              "精确相消且等能的插值参数",t.str(12))


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
