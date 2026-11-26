"""从贝伦德壳现有两张图的边差作用求矢量谱隙及单边软化锐界。"""
from pathlib import Path
from runpy import run_path

from flint import arb_mat, fmpq, fmpq_mat, fmpz_mat

study = run_path(str(Path(__file__).with_name("ch39-noisy-direct-vector-gap.py")))
poles, q, window = (study[k] for k in ("poles","q","window"))
neutral = poles["neutral"]


def capture(vertices, lengths, trees):
    edges = [(i,j) for i,p in enumerate(vertices) for j,other in enumerate(vertices[:i])
             if neutral["graph"]["squared"](p,other) in lengths]
    n = len(vertices)
    lap = fmpz_mat(n,n)
    for i,j in edges:
        lap[i,i] += 1
        lap[j,j] += 1
        lap[i,j] -= 1
        lap[j,i] -= 1
    minor = fmpz_mat([[lap[i,j] for j in range(n-1)] for i in range(n-1)])
    assert minor.det() == trees
    return lap, edges


def least_positive(lap):
    polynomial = lap.charpoly()
    roots = sorted(polynomial.complex_roots(), key=lambda pair:float(pair[0].real.mid()))
    assert roots[0][0].real == 0 and roots[0][1] == 1
    assert all(r.imag == 0 and r.real > 0 for r,_ in roots[1:])
    value,multiplicity = roots[1]
    assert multiplicity == 2
    return value.real


def retained_weight(lap, edges, ratio):
    n = lap.nrows()
    projector = fmpq_mat([[fmpq(1,n) for _ in range(n)] for _ in range(n)])
    shifted = arb_mat(lap)+arb_mat(projector)
    for i in range(n):
        shifted[i,i] -= ratio
    inverse = shifted.inv()
    # 对固定的差源 b=e_i-e_j 作秩一降阶；退到等号时得到恰在阈值的模。
    limits = [1-1/(inverse[i,i]+inverse[j,j]-2*inverse[i,j]) for i,j in edges]
    return max(v.lower() for v in limits), max(v.upper() for v in limits)


def verify():
    unit,_,colour,right = neutral["graph_constants"](port_builder=capture)
    cases = []
    time = poles["study"]["nominal"][2]
    for logs in poles["retune"]["log_sets"]:
        (base, _, _, _, _, masses), _, _ = poles["retune"]["predictions"](logs,time)
        m1,m2 = poles["reference_readings"](logs)[:2]
        error = window(q("1e-11"))
        length,_,_ = poles["charged_geometry"](m1*error,m2*error)
        triplet = masses[0]*(base/length).sqrt()*window(q("1e-7"))
        cases.append(4*triplet**2/unit)

    for name,(lap,edges) in (("colour",colour),("right",right)):
        value = least_positive(lap)
        print(name, "first positive eigenvalue", value.str(12), flush=True)
        for ratio in cases:
            assert 0 < ratio < min(q(1), value)
            ordinary_mass = 2*(value/ratio).sqrt()
            uniform = ratio/value
            weak = retained_weight(lap, edges, ratio)
            if name == "colour":
                assert ordinary_mass > q("4.395") and uniform < q(".208") and weak[1] < 0
            else:
                assert ordinary_mass > q("3.986") and uniform < q(".252")
                assert q(".1240") < weak[0] <= weak[1] < q(".1245")
            print("mass/triplet",ordinary_mass.str(10),"uniform threshold",uniform.str(10),
                  "single-link threshold", weak[0].str(9), weak[1].str(9), flush=True)


if __name__ == "__main__":
    verify()
