"""十九站原壳在所有钟相位与软边亏额下的谱宽界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def spectral_bounds(root):
    """逐荷返回完整原壳的 (最低值下界，最高值上界，最高/最低值比上界)。"""
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 160
    edges = shell["geo"]["edges"]
    pole = shell["soft"]["pole"]
    tuned = set(shell["soft"]["star"]["orbit"])
    vertices = len(shell["geo"]["sites"])
    spokes = {j if i==pole else i for i,j,_,_ in edges if pole in (i,j)}
    assert pole==14 and spokes==set(range(vertices))-{pole} and vertices==19
    assert len(edges)==135 and sum(d in (724,804) for _,_,d,_ in edges)==12
    degree = [sum(i==v or j==v for i,j,_,_ in edges) for v in range(vertices)]
    maximum = arb(2*max(degree))
    assert maximum==36
    middle = [(d,z) for i,j,d,z in edges if (i,j)==(0,16)]
    assert len(middle)==1 and middle[0][0]==632
    assert all(d not in (724,804) for i,j,d,_ in edges if (i,j) in ((0,14),(14,16)))
    least_weight = 1-arb("1e-8")
    result = {}
    for q in (1,-4,2,-3,6):
        mismatch = abs(1-middle[0][1]**q)
        assert mismatch>arb(".0868")
        path_floor = least_weight/(1+arb(3*vertices).sqrt()/mismatch)**2
        assert maximum/path_floor<279000
        fixed = acb_mat(vertices,vertices)
        for i,j,d,z in edges:
            if d in (724,804):
                continue
            spoke = pole in (i,j)
            phase = acb(1) if spoke else z**q
            weight = least_weight if spoke and (i in tuned or j in tuned) else arb(1)
            fixed[i,i] += weight
            fixed[j,j] += weight
            fixed[i,j] -= weight*phase
            fixed[j,i] -= weight*phase.conjugate()
        eigenvalues = fixed.eig(multiple=True)
        assert all(v.imag.contains(0) for v in eigenvalues)
        floor = min(v.real.lower() for v in eigenvalues)
        ratio = (maximum/floor).upper()
        assert floor>arb(".03403") and ratio<1058
        result[q] = floor,maximum,ratio
    return result


if __name__=="__main__":
    for q,(floor,upper,ratio) in spectral_bounds(Path(__file__).resolve().parent).items():
        print("荷、磁谱底下界、谱顶上界、谱宽比上界",q,floor.str(12),upper,ratio.str(12))
