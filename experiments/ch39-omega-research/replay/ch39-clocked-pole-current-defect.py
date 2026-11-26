"""同一双时钟电荷五传播核的轮换保护、单点破缺与有限动量检验。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

joint = run_path(str(Path(__file__).with_name("ch39-joint-path-minimum.py")))
study = joint["spectrum"]["study"]
geo = study["geo"]
sites = geo["sites"]
pole = sites.index(joint["pole"])
size = len(sites)
rotation = [sites.index((y, z, x)) for x, y, z in sites]
orbits = [(i, rotation[i], rotation[rotation[i]]) for i in range(size)
          if i != pole and i == min(i, rotation[i], rotation[rotation[i]])]
identity = acb_mat([[int(i == j) for j in range(size)] for i in range(size)])


def geometry():
    # Integer edge signature in the pole frame, before normalization and clocks.
    signed = {}
    for i, j, d, _ in geo["edges"]:
        x, y = sites[i], sites[j]
        real = geo["R"]+geo["dot"](joint["pole"], x)+geo["dot"](x, y)+geo["dot"](y, joint["pole"])
        imag = geo["det"](joint["pole"], x, y)
        clock = d if d in (724, 804) else 0
        signed[i, j] = (d, real, imag, clock)
        signed[j, i] = (d, real, -imag, -clock)
    assert len(signed) == 270
    assert all(signed[rotation[i], rotation[j]] == sig for (i, j), sig in signed.items())
    for j in range(size):
        if j != pole:
            d, real, imag, clock = signed[pole, j]
            assert real > 0 and imag == clock == 0 and d != 200


def charge_blocks(height):
    phases = geo["phases"](height)
    base = acb_mat(size, size)
    extras = acb_mat(size, size)
    for i, j, d, z in geo["edges"]:
        block = extras if d in phases else base
        block[i, i] += 1
        block[j, j] += 1
        block[i, j] = -(z*phases.get(d, 1))**5
        block[j, i] = block[i, j].conjugate()
    # Each omitted clocked edge adds a positive semidefinite rank-one form.
    lower = base-arb("1.2908")*identity
    for k in range(1, size+1):
        minor = acb_mat([[lower[i, j] for j in range(k)] for i in range(k)])
        assert minor.det().real > 0
    return base+extras


def second_current():
    a, b = joint["a"], joint["b"]
    assert joint["word"](joint["baryon_loops"]) == (-35, 13)
    phase = (acb(-83, 13*arb(3).sqrt())/86)**5
    vec = acb_mat(size, 1)
    i, j, k = orbits[0]
    vec[i, 0] = arb(a)/(a+b)
    vec[j, 0] = b*phase/(2*(a+b))
    vec[k, 0] = b*phase.conjugate()/(2*(a+b))
    return vec


def finite_samples(laplacian):
    basis = [(i, last) for first, second, last in orbits for i in (first, second)]
    embedding = acb_mat([[int(v == i)-int(v == j) for i, j in basis] for v in range(size)])
    momenta = [0, arb(1)/4, 1, 2, 4, 7, 11, 16, 24, 40, 80, 160]
    greens = [(laplacian+s*identity).inv() for s in momenta]
    assert all(g[i, i].real < arb(".126") for g in greens for i in range(size))
    least = 1
    for mark, _, _ in orbits:
        # The rows omit only the scalar -delta/(1+delta*G_xx) of the exact defect.
        samples = acb_mat([[g[pole, mark]*(g[mark, i]-g[mark, j]) for i, j in basis]
                           for g in greens])
        recover = embedding*samples.inv()
        upper = sum(sum(abs(recover[i, j]) for j in range(len(basis)))**2
                    for i in range(size)).sqrt()
        least = min(least, (1/upper).lower())
    assert least > arb("4.75e-23")
    print("Twelve-point uniform inverse bound over six mark orbits", least.str(12))
    return greens[0]


def verify():
    ctx.prec = 384
    assert len(orbits) == 6 and size == 19
    geometry()
    laplacian = charge_blocks(acb.zeta_zero(11431).imag)
    green = finite_samples(laplacian)
    current = second_current()
    delta = -arb(2)**-20/study["center"]
    for perturbation in (arb(1), delta):
        changed = acb_mat(laplacian)
        changed[0, 0] += perturbation
        amplitude = (changed.solve(current))[pole, 0]
        rank_one = -perturbation*green[pole, 0]*(green*current)[0, 0]/(1+perturbation*green[0, 0])
        assert abs(amplitude-rank_one) < arb("1e-45")
        print("Local diagonal, static absolute amplitude", perturbation.str(11), abs(amplitude).str(16))
        if perturbation > 0:
            assert arb(".00051481424107") < abs(amplitude) < arb(".00051481424109")
        else:
            assert arb("4.49205544049e-18") < abs(amplitude) < arb("4.49205544050e-18")
    # Distribute the same total mass over all three sites.
    restored = acb_mat(laplacian)
    for i in orbits[0]:
        restored[i, i] += delta/3
    assert abs(restored.solve(current)[pole, 0]) < arb("1e-46")


if __name__ == "__main__":
    verify()
