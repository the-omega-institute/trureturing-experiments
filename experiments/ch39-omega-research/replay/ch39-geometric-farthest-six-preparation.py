"""沿最短边把原三源送到全部六个最远站的五荷统一制备界。"""
from pathlib import Path
from runpy import run_path
from flint import arb, acb, acb_mat, ctx

ROOT = Path(__file__).resolve().parent
def load(name):
    return run_path(str(ROOT/name))

sh = load('ch39-three-ray-chiral-boundary.py')
pump = load('ch39-native-light-projector-pump.py')
gram = load('ch39-native-rh-ray-pvm.py')
least = load('ch39-terminal-port-light-bands.py')["least_singular"]
ctx.prec = 192
sites = sh["geo"]["sites"]
pole = sh["soft"]["pole"]
d = lambda i,j: sum((x-y)**2 for x,y in zip(sites[i],sites[j]))
B = sh["boundary_sites"]()
assert set(B) == {i for i in range(19) if d(i,pole)==min(d(j,pole) for j in range(19) if j!=pole)}
S = [i for i in range(19) if d(i,pole)==max(d(j,pole) for j in range(19))]
assert S == list(range(6)) and d(S[0],pole) == 544
shortest = min(d(i,j) for i in S for j in B)
assert shortest == 244
pairs = [(i,j) for i in S for j in B if d(i,j)==shortest]
assert [sum(a==i for a,b in pairs) for i in S] == [1]*6
assert [sum(b==j for a,b in pairs) for j in B] == [2]*3
rotate = [sites.index((y,z,x)) for x,y,z in sites]
reflect = [sites.index((x,z,y)) for x,y,z in sites]
T = {0,rotate[0],rotate[rotate[0]]}
assert T == {0,3,4}
assert {reflect[i] for i in T} == set(S)-T
adj = pump["adjoint"]
C = acb_mat([row[1:] for row in gram["gram_cholesky"](ROOT)[1:]]).transpose()
s = least(C)
assert s > arb('.0000899007')
omega = acb(-1,arb(3).sqrt())/2
turns = {B[0]:0,rotate[B[0]]:1,rotate[rotate[B[0]]]:2}
assert set(turns)==set(B)
for q,t in zip((1,-4,2,-3,6),(".5753",".5193",".5210",".5744",".5165")):
    base, clocks = pump["clock_matrices"](sh,q)
    L0 = base
    W = acb_mat([[base[i,j]/arb(2).sqrt() if (i,j) in pairs else acb(0)
                  for j in B] for i in range(19)])
    overlaps = []
    for k,F in enumerate(pump["rotation_bases"](sh)):
        v,_ = pump["first_line"](adj(F)*L0*F)
        v = F*acb_mat([[z] for z in v])
        b = acb_mat([[omega**(-k*turns[j])/arb(3).sqrt()] for j in B])
        overlaps.append(abs((adj(v)*W*b)[0,0]))
    print(q,[v.str(15) for v in overlaps],flush=True)
    assert overlaps[0]>arb(t) and min(overlaps[1:])>arb(".9092")

def F(m,s):
    m,s = arb(m),arb(s)
    return m*(1-s*s).sqrt()-s*(1-m*m).sqrt()
k0=F(".5165",".17")
k12=F(".9092",arb(1)/3)
assert k0>arb(".3634") and k12>arb(".71")
assert k0*arb(".0000899007")>arb(".0000326")
print("lower bounds",k0.str(18),k12.str(18),(k0*arb(".0000899007")).str(18))
