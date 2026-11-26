"""在十九点贝伦德球壳上用六次几何三角输运制备唯一基态的端口。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

port = run_path(str(Path(__file__).with_name("ch39-single-port-ground-selection.py")))
geo = port["study"]["geo"]


def verify():
    ctx.prec = 192
    edges = {(i,j):d for i,j,d,_ in geo["edges"]}
    def unit(i,j,k):
        return geo["unit"](i,j)*geo["unit"](j,k)*geo["unit"](i,k).conjugate()
    count = unclocked = 0
    for i,j,k in combinations(range(len(geo["sites"])),3):
        triangle = ((i,j),(j,k),(i,k))
        if all(edge in edges for edge in triangle):
            angle = abs(unit(i,j,k).arg())
            assert angle < arb(".156547")
            count += 1
            if all(edges[edge] not in (724,804) for edge in triangle):
                assert angle < arb(".108207") or angle > arb(".156545")
                unclocked += 1
    assert 5*arb(".156547") < arb(".892")
    assert 5*arb(".156547")+arb(".108207") < arb(".892")
    assert 6*arb(".156545") > arb(".920")
    assert 5*arb(".156547")-arb(".156545") < arb(".892")
    assert arb.pi()-6*arb(".156547") > arb(".920")
    print("Admissible triangles, all and without clocked edges",count,unclocked)
    for loops,counts,clocks_off in [
            ([(5,13,15),(5,7,13),(5,6,7)],(3,2,1),False),
            ([(4,5,14),(5,13,15),(5,14,15)],(1,5,1),True)]:
        run_word(loops,counts,clocks_off,edges,unit)


def run_word(loops,counts,clocks_off,edges,unit):
    factor = acb(1)
    for vertices, times in zip(loops,counts):
        i,j,k = vertices
        triangle = ((i,j),(j,k),(i,k))
        assert all(edge in edges for edge in triangle)
        if clocks_off:
            assert all(edges[edge] not in (200,724,804) for edge in triangle)
        z = unit(i,j,k)
        if z.imag < 0:
            z = z.conjugate()
        assert z.imag > 0
        factor *= z**times
        print("Oriented triangle (positive phase)",[geo["sites"][a] for a in vertices],
              "angle",z.arg().str(13),"repetitions",times)
    ratio = (arb(25)/6).sqrt()
    vector = acb_mat([[(arb(6)/31).sqrt()],[factor*5/arb(31).sqrt()]])
    assert arb(".90572") < factor.arg() < arb(".90577")
    assert abs(abs(factor)-1)<arb("1e-34")
    print("Preparation phase",factor.arg().str(14),"amplitude ratio",ratio.str(14))
    gap,full_gap,energy = port["verify"](vector,allowance=arb("3e-9"))
    assert gap > arb("5.12e-9") and full_gap > arb("2.92e-15")
    assert energy > -arb(2)**-20


if __name__=="__main__":
    verify()
