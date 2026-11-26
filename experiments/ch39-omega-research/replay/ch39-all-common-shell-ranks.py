"""分类一个固定实截断在全部正零点、五种荷上能保留的非平凡内部重数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def character_bases(shell):
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    orbits = [[i] if i == j else [i,j,rotation[j]] for i,j in enumerate(rotation)
              if i == min(i,j,rotation[j])]
    omega = acb(-arb(1)/2, arb(3).sqrt()/2)
    for k in range(3):
        active = [o for o in orbits if len(o) == 3 or k == 0]
        basis = acb_mat(19,len(active))
        for j,orbit in enumerate(active):
            for r,i in enumerate(orbit):
                basis[i,j] = omega**(k*r)/arb(len(orbit)).sqrt()
        yield basis


def build_family(shell, charge):
    base = acb_mat(19,19)
    clock = [acb_mat(19,19),acb_mat(19,19)]
    for i,j,d,z in shell["geo"]["edges"]:
        value = -acb(1) if shell["soft"]["pole"] in (i,j) else -z**charge
        base[i,i] += 1
        base[j,j] += 1
        if d in (724,804):
            assert i < 6 and j in (6,7,8)
            clock[(724,804).index(d)][i,j] = value
        else:
            base[i,j] = value
            base[j,i] = value.conjugate()
    return base,clock


def character_cuts(base, clock, basis, cuts):
    adj = lambda a: a.conjugate().transpose()
    block = adj(basis)*base*basis
    terms = [adj(basis)*t*basis for t in clock]
    center, = [j for j in range(basis.ncols()) if not basis[6,j].is_zero()]
    outer = [j for j in range(basis.ncols()) if j != center]
    fixed = acb_mat([[block[i,j] for j in outer] for i in outer])
    levels = [v.real for v in fixed.eig(multiple=True)]
    links = [acb_mat([[t[i,center]] for i in outer]) for t in [block,*terms]]
    result = []
    for t in cuts:
        shifted = acb_mat(fixed)
        for i in range(shifted.nrows()):
            shifted[i,i] -= t
        inverse = shifted.inv()
        g = [[(adj(a)*inverse*b)[0,0] for b in links] for a in links]
        scalar = (block[center,center]-t-sum(g[i][i] for i in range(3))).real
        assert not scalar.contains(0)
        sign = 1 if scalar > 0 else -1
        # For u on the first circle, optimize the second phase exactly.
        if not abs(scalar) > 2*sum(abs(g[i][j]) for i,j in ((0,1),(0,2),(1,2))):
            verify_circle_sign(scalar, sign, g)
        assert all(v < t or v > t for v in levels)
        result.append(sum(v < t for v in levels) + int(sign < 0))
    return tuple(result)


def verify_circle_sign(scalar, sign, g):
    count = 2048
    angular_radius = arb.pi()/count
    for i in range(count):
        u = acb(0,2*arb.pi()*i/count).exp()
        # u=e^{ix}, v=e^{iy}; coefficients of v are g02+conj(u)g12.
        constant = scalar-2*(u*g[0][1]).real
        amplitude = 2*abs(g[0][2]+u.conjugate()*g[1][2])
        variation = 2*(abs(g[0][1])+abs(g[1][2]))*angular_radius
        assert sign*constant-amplitude > variation, (i,sign,constant,amplitude,variation)


def verify():
    ctx.prec = 192
    shell = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))
    cuts = [arb(t) for t in ("3.1","7.1","8.69","9.59","12","14.8","18.32","18.35")]
    expected = [(1,1,1,1,3,3,5,5),(0,0,1,1,2,2,5,5),(0,0,1,1,2,2,5,5)]
    # The reduced matrices are exact cyclic restrictions. Testing each sign
    # includes the intervening circle arcs, not just the listed centers.
    bases = list(character_bases(shell))
    charges = (1,2,3,4,6)
    for q in charges:
        base,clock = build_family(shell,q)
        assert [character_cuts(base,clock,basis,cuts) for basis in bases] == expected
    heights = acb.zeta_zeros(11416,34)
    ordered = []
    for q in charges:
        for z in heights:
            levels = sorted((v.real for v in shell["laplacian"](q,z.imag).eig()), key=float)
            assert all(a < b for a,b in zip(levels,levels[1:]))
            ordered.append(levels)
    possible = {1,3,7,15}
    for rank in range(1,19):
        last_low = max((row[rank-1] for row in ordered), key=float)
        first_high = min((row[rank] for row in ordered), key=float)
        if rank not in possible:
            assert last_low > first_high
        elif rank != 1:
            assert first_high-last_low < arb("3.1")
    print("整张环面的共同截断端点", *[v.str(6) for v in cuts])
    print("全部非平凡重数与字符",[(1,(1,0,0)),(3,(1,1,1)),(7,(3,2,2)),(15,(5,5,5))])
    print("排除其余十四种重数的实际零点",11416,11449)


if __name__ == "__main__":
    verify()
