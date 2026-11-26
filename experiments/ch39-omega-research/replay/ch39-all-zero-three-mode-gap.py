"""三个内站点的 Schur 消元给两只弦钟整张环面上的三模共同谱隙。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify():
    directory = Path(__file__).parent
    shell = run_path(str(directory/"ch39-three-ray-chiral-boundary.py"))
    algebra = run_path(str(directory/"ch39-three-mode-higgs-hierarchy.py"))
    adjoint,identity = algebra["adjoint"],algebra["identity"]
    ctx.prec = 320
    inside = (6,7,8)
    outside = [i for i in range(19) if i not in inside]
    omega = acb(-arb(1)/2,arb(3).sqrt()/2)
    fourier = acb_mat([[omega**(i*k)/arb(3).sqrt() for k in range(3)] for i in range(3)])
    cut = arb("9.6001")
    for q in (1,2,3,4,6):
        flat_clock = acb_mat(19,19)
        links = [acb_mat(16,3) for _ in range(3)]
        for i,j,d,z in shell["geo"]["edges"]:
            value = -acb(1) if shell["soft"]["pole"] in (i,j) else -z**q
            flat_clock[i,i] += 1
            flat_clock[j,j] += 1
            flat_clock[i,j] = value
            flat_clock[j,i] = value.conjugate()
            clock = {724:1,804:2}.get(d,0)
            if clock:
                assert i in outside and j in inside
            if i in inside and j in outside:
                links[clock][outside.index(j),inside.index(i)] = value.conjugate()
            if i in outside and j in inside:
                links[clock][outside.index(i),inside.index(j)] = value
        principal = acb_mat([[flat_clock[i,j] for j in outside] for i in outside])
        levels = sorted((x.real for x in principal.eig()),key=float)
        assert levels[2] < arb("8.685") < cut < levels[3]
        inverse = (principal-cut*identity(16)).inv()
        center = acb_mat([[flat_clock[i,j] for j in inside] for i in inside])-cut*identity(3)
        inner = adjoint(fourier)*center*fourier
        transformed = [b*fourier for b in links]
        pieces = [[adjoint(a)*inverse*b for b in transformed] for a in transformed]
        # Cyclic covariance makes all these matrices diagonal in the same
        # basis; the interval tests independently check that decomposition.
        for m in [inner,*[item for row in pieces for item in row]]:
            assert all(m[i,j].contains(0) for i in range(3) for j in range(3) if i != j)
        margins = []
        for k in range(3):
            center_value = (inner[k,k]-sum(pieces[i][i][k,k] for i in range(3))).real
            radius = 2*sum(abs(pieces[i][j][k,k]) for i,j in ((0,1),(0,2),(1,2)))
            margins.append(center_value-radius)
        assert min(margins) > arb(".04")
        print("荷、外部主块第三值、Schur 最小相位一致下界",q,levels[2].str(10),min(margins).str(10))
        if q == 2:
            negative = acb_mat(flat_clock)
            for i,j,d,_ in shell["geo"]["edges"]:
                if d in (724,804):
                    negative[i,j],negative[j,i] = -negative[i,j],-negative[j,i]
            first = min(v.real for v in flat_clock.eig())
            last = min(v.real for v in negative.eig())
            assert first < arb(".042") and last > arb("2.035")
    assert cut-6*arb("1e-8") > arb("9.6")


if __name__ == "__main__":
    verify()
