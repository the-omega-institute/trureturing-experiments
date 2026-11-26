"""三低模夸克共用一个局域 Higgs 时，顶底层级的最低 Higgs 谱截断。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adjoint(a):
    return a.conjugate().transpose()


def cols(a, end):
    return acb_mat([[a[i,j] for j in range(end)] for i in range(a.nrows())])


def norm(a):
    return sum((z.conjugate()*z).real for row in a.tolist() for z in row).sqrt()


def positive(a):
    for k in range(1, a.nrows()+1):
        minor = acb_mat([[a[i,j] for j in range(k)] for i in range(k)]).det()
        assert minor.imag.contains(0) and minor.real > 0


def identity(n):
    m = acb_mat(n,n)
    for i in range(n):
        m[i,i] = 1
    return m


def verify():
    directory = Path(__file__).parent
    shell = run_path(str(directory/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 384
    height = acb.zeta_zero(11431).imag
    bases, energies = {}, {}
    for q in (1,4,-2,3,-3,-6):
        levels, modes = shell["laplacian"](q,height).eig(right=True)
        order = sorted(range(19),key=lambda j:float(levels[j].real))
        energies[q] = [levels[j].real for j in order]
        assert all(a < b for a,b in zip(energies[q],energies[q][1:]))
        if q != 3:
            order = order[:3]
        vectors = []
        for j in order:
            v = acb_mat([[modes[i,j]] for i in range(19)])
            vectors.append(v/norm(v))
        basis = acb_mat([[v[i,0] for v in vectors] for i in range(19)])
        assert all(z.contains(0) for row in (adjoint(basis)*basis-identity(len(order))).tolist()
                   for z in row)
        bases[q] = basis
    higgs = bases[3]
    maps = []
    for q,lq,conjugated in ((4,1,True),(-2,1,False),(-6,-3,False)):
        entries = acb_mat(9,19)
        for k in range(19):
            diagonal = acb_mat(19,19)
            for i in range(19):
                diagonal[i,i] = higgs[i,k]
            if conjugated:
                diagonal = diagonal.conjugate()
            matrix = adjoint(bases[lq])*diagonal*bases[q]
            if conjugated:
                matrix = matrix.conjugate()
            for i in range(3):
                for j in range(3):
                    entries[3*i+j,k] = matrix[i,j]
        maps.append(entries)
    up, down, lepton = maps
    for a in maps:
        ground = acb_mat([[a[3*i+j,0] for j in range(3)] for i in range(3)])
        values = sorted([v.real.sqrt() for v in (ground*adjoint(ground)).eig()],key=float)
        assert arb(".229") < values[0] and values[-1] < arb(".242")
        print("Higgs 基态的三个奇异值",*[v.str(10) for v in values])
    for count, ceiling in ((3,"1.728"),(9,"4.104")):
        u,d = cols(up,count),cols(down,count)
        positive(arb(ceiling)**2*(adjoint(d)*d)/3-adjoint(u)*u)
    d9 = cols(down,9)
    assert abs(d9.det()) > 0
    coordinate = -d9.solve(acb_mat([[down[i,9]] for i in range(9)]))
    witness = acb_mat(coordinate.tolist()+[[1]])
    witness = witness/norm(witness)
    assert all(z[0].contains(0) for z in (cols(down,10)*witness).tolist())
    surviving = cols(up,10)*witness
    surviving = acb_mat([[surviving[3*i+j,0] for j in range(3)] for i in range(3)])
    top = max(v.real.sqrt() for v in (surviving*adjoint(surviving)).eig())
    assert arb(".0078833") < top < arb(".0078835")
    assert abs(witness[9,0]) > arb(".9982")
    print("十模下型精确为零时的上型最大奇异值",top.str(12))

    right_inverse = adjoint(down)*(down*adjoint(down)).inv()
    transfer = up*right_inverse
    rest = up-transfer*down
    positive(arb("1.591")**2*identity(9)-transfer*adjoint(transfer))
    positive(arb(".03046")**2*identity(9)-rest*adjoint(rest))
    lower = (arb(".8573")-arb("1.591")*arb(3).sqrt()*arb(".01410"))/arb(".03046")
    assert lower > arb("26.869")
    print("完整十九个 Higgs 模式所需共同局域系数的下界",lower.str(12))

    ground = cols(higgs,1)
    equal = acb_mat([[ground[i,0]/(abs(ground[i,0])*arb(19).sqrt())] for i in range(19)])
    lap = shell["laplacian"](3,height)
    trial_energy = (adjoint(equal)*lap*equal)[0,0].real
    assert arb("1.57366") < trial_energy < arb("1.57368")
    levels = energies[3]
    excited = (trial_energy-levels[0])/(levels[1]-levels[0])
    assert excited < arb(".02211")
    relative = (excited/(1-excited)).sqrt()
    for a, ceiling in ((up,".45656"),(down,".45887")):
        higher = acb_mat([[a[i,j] for j in range(1,19)] for i in range(9)])
        positive(arb(ceiling)**2*identity(9)-higher*adjoint(higher))
    def ground_max(a):
        m = acb_mat([[a[3*i+j,0] for j in range(3)] for i in range(3)])
        return max(v.real.sqrt() for v in (m*adjoint(m)).eig())
    upper_ratio = ((ground_max(up)+arb(".45656")*relative) /
                   (ground_max(down)-arb(".45887")*relative))
    assert upper_ratio < arb("1.796")
    print("共用径向凸势的激发份额上界与顶底比上界",excited.str(10),upper_ratio.str(10))

    # The hierarchy implies x*(D*D - 3 U*U/R²)x <= 0. A single positive
    # multiplier suffices to certify a lower bound for all such profiles.
    required = arb(".8573")/arb(".01410")
    kinetic = acb_mat(19,19)
    for i,e in enumerate(levels):
        kinetic[i,i] = e
    constrained = kinetic + 30000*(adjoint(down)*down - 3*adjoint(up)*up/required**2)
    positive(constrained-arb("14.8147")*identity(19))
    assert arb("14.8147")-trial_energy > arb("13.241")
    print("层级所需动能下界、相对等模长比较态的能量代价下界",14.8147,13.241)


if __name__ == "__main__":
    verify()
