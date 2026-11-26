"""Mass-ratio obstruction for every angle and every positive single-orbit source."""
from functools import lru_cache
from itertools import product
from pathlib import Path
import runpy
from flint import acb, acb_poly, arb

ray=runpy.run_path(str(Path(__file__).with_name("ch39-orbit-flux-windows.py")))
fit, wall = ray["fit"], ray["wall"]
modes=((0,0),(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1))


def coefficients(orbit):
    # M_JJ and H_J are clock independent. B and H_W are affine in the two
    # conjugated phases. Their Hermitian products therefore have exactly these harmonics.
    samples={}
    for i,j in product(range(3),repeat=2):
        model,orbits=fit["model"](2*arb.pi()*i/3,2*arb.pi()*j/3)
        _,mats=ray["wall_samples"](model,orbits[orbit])
        def sample(mat):
            syms=fit["symbols"](fit["block"](mat,fit["ii"],fit["oo"]))
            return syms+[mat[wall.index(7),wall.index(6)],mat[wall.index(11),wall.index(12)]]
        samples[i,j]=list(map(sample,mats))
    rows=[]
    for v in range(5):
        row={}
        for a,b in modes:
            vals=[sum(data[k][v]*acb(0,-2*arb.pi()*(i*a+j*b)/3).exp()
                      for (i,j),data in samples.items())/9 for k in range(4)]
            poly=acb_poly.interpolate(ray["nodes"],vals)
            if (v<3 and (a,b) not in ((0,0),(-1,0),(0,-1))) or (v==4 and (a,b)!=(0,0)):
                assert all(z.contains(0) for z in poly.coeffs())
            else:
                row[a,b]=poly
        rows.append(row)
    return rows


@lru_cache(maxsize=2048)
def clock(level,index):
    turns=(arb(index)+arb(1)/2+arb(0,arb(1)/2))/(1<<level)
    return acb(0,2*arb.pi()*turns).exp()


def cover(rows,threshold):
    @lru_cache(maxsize=4096)
    def stiffness(level,index):
        box=arb(index).union(arb(index+1))/(1<<level)
        return [{k:poly(acb(box)) for k,poly in row.items()} for row in rows]
    for j in range(256):
        base=stiffness(8,j)[4][0,0]
        assert arb(".034")<base.arg()<arb(".038")
    todo=[((0,0,0),(0,0,0))]
    leaves=0
    deepest=0
    while todo:
        levels,idx=todo.pop()
        p,q=(clock(levels[i],idx[i]) for i in (0,1))
        phases={(0,0):1,(1,0):p,(-1,0):p.conjugate(),(0,1):q,(0,-1):q.conjugate(),
                (1,-1):p*q.conjugate(),(-1,1):q*p.conjugate()}
        values=[sum(c*phases[k] for k,c in row.items()) for row in stiffness(levels[2],idx[2])]
        syms,top,base=values[:3],values[3],values[4]
        assert base.real>0 and base.imag>0 and base.real>2*base.imag
        # With nonzero spokes and degree three, 0<arg(X)<pi and arg(X)+arg(Y)>pi.
        if top.imag<0 or (top*base).imag>0:
            leaves+=1
            continue
        sizes=[abs(z) for z in syms]
        cap=max(v.upper() for v in sizes)
        if sum(v.lower()>threshold*cap for v in sizes)>=2:
            leaves+=1
            continue
        # Split all three compact dimensions fairly, with initially finer clock widths.
        dim=min(range(3),key=lambda j:levels[j]+(2 if j==2 else 0))
        level=levels[dim]
        assert level<19,(levels,idx)
        deeper=list(levels)
        deeper[dim]+=1
        deepest=max(deepest,deeper[dim])
        for bit in (0,1):
            nextidx=list(idx)
            nextidx[dim]=2*idx[dim]+bit
            todo.append((tuple(deeper),tuple(nextidx)))
    return leaves,deepest


if __name__=="__main__":
    for orbit in range(4):
        print("orbit, boxes, maximum level",orbit,*cover(coefficients(orbit),arb(".05")))
