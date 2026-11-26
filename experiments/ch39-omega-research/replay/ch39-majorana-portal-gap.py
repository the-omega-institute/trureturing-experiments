"""Full-energy low quasiparticle gap of a particle-hole-paired local shell."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

native = run_path(str(Path(__file__).with_name("ch39-native-flavor-eigenaxes.py")))
inside,origin,portals = (native[k] for k in ("inside","origin","portals"))
geometry,columns,source = (native[k] for k in ("geo","columns","source"))
unit = acb_mat([[int(i==j) for j in inside] for i in inside])


def parent_blocks(phases):
    m = geometry["matrix"](0, phases)
    block = acb_mat([[m[i,j] for j in inside] for i in inside])
    f = columns*acb_mat([[m[origin,portals[0]],0],[0,m[origin,portals[1]]]])
    return block, f


def parent_data(height):
    block, f = parent_blocks(geometry["phases"](height))
    inv = block.inv()
    assert sum(abs(v)**2 for v in inv.entries()) < 1
    t0 = -f*(source.transpose()*inv*source)*f.conjugate().transpose()
    assert 0 < abs(t0[0,1].imag) < arb(".025")
    return block,f


def sign_at(data,energy):
    m,f = data
    def ret(shift):
        compressed = source.transpose()*(m-shift*unit).solve(source)
        return -f*compressed*f.conjugate().transpose()
    k = (ret(energy)-ret(-energy).conjugate())/2
    mean = k.trace().real/2
    top = mean+(((k[0,0]-k[1,1]).real/2)**2+abs(k[0,1])**2).sqrt()
    diff = top-energy
    assert 0 not in diff
    return diff > 0


def verify():
    best = {}
    for idx,z in enumerate(acb.zeta_zeros(10815,1278),10815):
        data = parent_data(z.imag)
        if sign_at(data,arb(".02350")):
            best[idx] = data
    assert set(best)=={11431,10965}
    brackets = {}
    for idx,data in best.items():
        lo,hi = arb(".0235"),arb(".025")
        for _ in range(32):
            m = ((lo+hi)/2).mid()
            if sign_at(data,m):
                lo = m
            else:
                hi = m
        brackets[idx] = lo.union(hi)
    assert brackets[11431]-brackets[10965] > arb(".0000928")
    print("full 38-mode gap at 11431",brackets[11431].str(13),
          "second at 10965",brackets[10965].str(13),"others below .02350")


if __name__ == "__main__":
    verify()
