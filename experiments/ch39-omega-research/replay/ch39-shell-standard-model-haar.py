"""Twenty-two independent shell angles, two Abelian holonomies and two Haar color matrices."""
from pathlib import Path
import runpy

from flint import arb, acb, acb_mat, fmpz

root = Path(__file__).resolve().parent
sampler = runpy.run_path(str(root/"ch39-zero-haar-clocks.py"))
geo = sampler["geo"]
theta = geo["theta"]
pi = arb.pi()
distances = sorted({d for _,_,d,_ in geo["edges"]})
classes = {}
for d in distances:
    signature = tuple(int(p) for p,e in fmpz(d*(3468-d)).factor() if e%2)
    classes.setdefault(signature,[]).append(d)
    assert 2*(1734-d)%1734 != 0
assert len(classes)==21 and classes[(2,)]==[12,68,136]
assert all(len(v)==1 for k,v in classes.items() if k!=(2,))
# In Q(sqrt(-2)), alpha_12,alpha_68,alpha_136 have (real,imaginary,denominator)
# (287,24,289),(49,10,51),(47,14,51).
assert 287*49-2*24*10 == 47*289 and 287*10+24*49 == 14*289
assert (49-10)%3==0 and (49+10)%3!=0 and 289%3 and (287-24)%3
assert len(distances)-1==22
color = [d for d in distances if d not in (12,32,52,136,200,724,804)]
assert len(color)==16
def unit_interval(phase):
    turns=phase/(2*pi)
    return turns-turns.floor()

def su3(phases):
    u,v = map(unit_interval,phases[:2])
    remainder = (1-u).sqrt()
    lengths = [(1-remainder).sqrt(),(remainder*v).sqrt(),(remainder*(1-v)).sqrt()]
    column = [s*acb(0,a).exp() for s,a in zip(lengths,phases[2:5])]
    x,y,z = column
    norm = (lengths[0]**2+lengths[1]**2).sqrt()
    assert norm>0
    vcol = [-y.conjugate()/norm,x.conjugate()/norm,acb(0)]
    # Complex conjugate of the bilinear cross product.
    wcol = [(column[(i+1)%3]*vcol[(i+2)%3]-column[(i+2)%3]*vcol[(i+1)%3]).conjugate()
            for i in range(3)]
    frame = acb_mat([[column[i],vcol[i],wcol[i]] for i in range(3)])
    h = unit_interval(phases[5])
    a = h.sqrt()*acb(0,phases[6]).exp()
    b = (1-h).sqrt()*acb(0,phases[7]).exp()
    return frame*acb_mat([[1,0,0],[0,a,-b.conjugate()],[0,b,a.conjugate()]])

if __name__=="__main__":
    print("independent from 23 lengths",22,"color clock representatives",color)
    for index in sampler["first_hits"]():
        gamma = acb.zeta_zero(index).imag
        matrices = [su3([gamma*theta(d) for d in color[start:start+8]]) for start in (0,8)]
        for g in matrices:
            product = g.conjugate().transpose()*g
            assert abs(g.det()-1)<arb("1e-40")
            assert all(abs(product[i,j]-int(i==j))<arb("1e-40") for i in range(3) for j in range(3))
        if index==8:
            u,v = matrices
            comm = u*v*u.conjugate().transpose()*v.conjugate().transpose()
            print("first accepted color commutator trace",sum(comm[j,j] for j in range(3)).str(12))
