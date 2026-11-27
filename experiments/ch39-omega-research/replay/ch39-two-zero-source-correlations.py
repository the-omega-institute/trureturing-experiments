"""在两个已选实际零点的整盒上，六质量迫使三源系数出现非对角相关。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def shape(masses):
    product = masses[0]*masses[1]*masses[2]
    return sum(x*x for x in masses)/product**(arb(2)/3)


def certified_inputs():
    model = run_path(str(Path(__file__).with_name("ch39-two-zero-quark-window.py")))
    shell, load, sources = model["shell"], model["load"], model["sources"]
    adjoint = model["adj"]
    root = Path(__file__).resolve().parent
    ports = [shell["soft"]["pole"], *shell["boundary_sites"]()]
    gram = load("native-rh-ray-pvm")["gram_cholesky"](root)
    raw = acb_mat([[gram[j][ports.index(i)] if i in ports else 0
                    for j in (1,2,3)] for i in range(19)])
    matrices = [model["prepare"](q,arb(x),arb(y)) for q,x,y in (
        (-4,"-0.680061861483031","0.873223965353125"),
        (2,"0.565247837594779","-0.630343671648397"),
    )]
    factor = abs(matrices[1].det()/matrices[0].det())**(arb(2)/3)
    masses = [m+arb(0,r.upper()) for m,r in zip(model["ref"][:6],model["tol"][:6],strict=True)]
    required = shape(masses[:3])/shape(masses[3:])
    assert required > arb("74.25")
    coordinates = (adjoint(sources)*sources, adjoint(sources)*raw)
    return matrices, coordinates, factor, adjoint


def verify():
    matrices, coordinates, factor, adjoint = certified_inputs()
    for transform in coordinates:
        outputs = [a*transform for a in matrices]
        up, down = [adjoint(a)*a for a in outputs]
        ratios = [factor*up[i,i].real/down[i,i].real for i in range(3)]
        assert all(v < arb("1.083") for v in ratios)
        print("三条源分别给出的归一化迹之比", [x.str(9) for x in ratios])
    condition = arb("74.25")/arb("1.083")
    correlation = (condition-1)/(2*(condition+1))
    assert condition > arb("68.55") and correlation > arb(".4856")


if __name__ == "__main__":
    verify()
