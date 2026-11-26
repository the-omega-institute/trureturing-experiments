"""三射线边界的轻模延拓、局部 Yukawa 投影与代际质量比。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

boundary_model = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))


def adjoint(a):
    return a.conjugate().transpose()


def extension(charge, height, boundary):
    lap = boundary_model["laplacian"](charge, height)
    interior = [i for i in range(19) if i not in boundary]
    principal = acb_mat([[lap[i,j] for j in interior] for i in interior])
    coupling = acb_mat([[lap[i,j] for j in boundary] for i in interior])
    harmonic = -principal.solve(coupling)
    result = acb_mat(19,len(boundary))
    for j,site in enumerate(boundary):
        result[site,j] = 1
    for k,site in enumerate(interior):
        for j in range(len(boundary)):
            result[site,j] = harmonic[k,j]
    equations = acb_mat([[lap[i,j] for j in range(19)] for i in interior])
    assert all(z.contains(0) for row in (equations*result).tolist() for z in row)
    return result


def squared_yukawa(left, right, scalar):
    local = acb_mat(19,19)
    for i in range(19):
        local[i,i] = scalar[i,0]
    y = left.transpose()*local*right
    left_metric = (adjoint(left)*left).transpose()
    right_metric = adjoint(right)*right
    # This form is similar to the positive Hermitian matrix in canonical
    # coordinates; it avoids an arbitrary orthonormalization of the three modes.
    form = left_metric.solve(y*right_metric.solve(adjoint(y)))
    roots = form.eig(multiple=True)
    assert all(z.imag.contains(0) and z.real>0 for z in roots)
    masses = sorted((z.real.sqrt() for z in roots), key=float)
    assert masses[0] < masses[1] < masses[2]
    return form, masses


def projectors(form, masses):
    identity = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    squares = [v*v for v in masses]
    result = []
    for i, eigen in enumerate(squares):
        p = identity
        for j, other in enumerate(squares):
            if i != j:
                p = p*(form-other*identity)/(eigen-other)
        result.append(p)
    return result


def verify():
    ctx.prec = 224
    height = acb.zeta_zero(11431).imag
    boundary = boundary_model["boundary_sites"]()
    geo = boundary_model["geo"]
    sites = geo["sites"]
    cycle = [sites.index((y,z,x)) for x,y,z in sites]
    assert set(cycle[i] for i in boundary) == set(boundary) and cycle[boundary[0]] != boundary[0]
    tuned = boundary_model["soft"]["star"]["orbit"]
    assert set(cycle[i] for i in tuned) == set(tuned)
    assert cycle[boundary_model["soft"]["pole"]] == boundary_model["soft"]["pole"]
    # A proper coordinate rotation preserves the dot and determinant defining
    # geometric transport. The clock-labelled edges also preserve orientation.
    edge_lengths = {(i,j):d for i,j,d,_ in geo["edges"]}
    for (i,j), d in edge_lengths.items():
        a,b = cycle[i],cycle[j]
        assert edge_lengths[tuple(sorted((a,b)))] == d
        if d in (724,804):
            assert a < b
    profiles = {q: extension(q, height, boundary) for q in (1,-4,2,-3,6)}
    higgs = extension(3, height, [boundary_model["soft"]["pole"]])
    norm = sum(abs(higgs[i,0])**2 for i in range(19)).sqrt()
    higgs = higgs/norm
    up, yu = squared_yukawa(profiles[1], profiles[-4], higgs)
    down, yd = squared_yukawa(profiles[1], profiles[2], higgs.conjugate())
    _, ye = squared_yukawa(profiles[-3], profiles[6], higgs.conjugate())
    for name, values in (("上型",yu),("下型",yd),("带电轻子",ye)):
        assert arb(".80") < values[0]/values[-1] < arb(".87")
        print(name, "归一轻 Yukawa、最小/最大", [x.str(9) for x in values],
              (values[0]/values[-1]).str(9))
    up_projectors, down_projectors = projectors(up,yu), projectors(down,yd)
    for i,p in enumerate(up_projectors):
        row = []
        for j,q in enumerate(down_projectors):
            mixed = p*q
            value = sum(mixed[k,k] for k in range(3))
            assert value.imag.contains(0)
            assert abs(value.real-int(i == j)) < arb("1e-25")
            row.append(value.real)
        assert sum(row).contains(1)
        print("CKM 模平方", [x.str(9) for x in row])


if __name__ == "__main__":
    verify()
