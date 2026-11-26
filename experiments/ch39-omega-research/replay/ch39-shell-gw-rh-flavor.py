"""Certified chiral index changes along the two R=26 clocks in noncommuting flavor links."""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx


folder=Path(__file__).resolve().parent
model=run_path(str(folder/"ch39-shell-gw-chiral-pair.py"))
ctx.prec=320
gaps={1:(3,-1,2),4:(3,-3,0),13:(2,-1,2)}
edges=(((0,3),(0,1)),((1,4),(1,2)))


def family_matrix(base,p,q):
    # The two modified underlying edges are disjoint.
    links={}
    for ((i,j),(a,b)),t in zip(edges,(p,q)):
        u=[[acb(int(x==y)) for y in range(3)] for x in range(3)]
        u[a][a]=u[b][b]=acb(t.cos())
        u[a][b]=u[b][a]=acb(0,t.sin())
        links[i,j]=u
        links[j,i]=[[u[y][x].conjugate() for y in range(3)] for x in range(3)]
    out=[]
    for i in range(36):
        a,v=divmod(i,3)
        row=[]
        for j in range(36):
            b,w=divmod(j,3)
            link=links.get((a%6,b%6))
            block=link[v][w] if link is not None else acb(int(v==w))
            value=base[a][b]*block
            if i==j:
                value+=arb(3)/4 if a<6 else -arb(3)/4
            row.append(value)
        out.append(row)
    return out


def opposite_modes(k,index,gap):
    n=len(k)
    unit=acb_mat([[int(i==j) for j in range(n)] for i in range(n)])
    scale=arb(8)
    assert max(sum(abs(v) for v in row).upper() for row in k)<scale
    x=acb_mat(k)/scale
    b=gap/scale
    while 1-b>arb("1e-12"):
        x=x*(3*unit-x*x)/2
        b=b*(3-b*b)/2
    start,sign=(n//2,1) if index>=0 else (0,-1)
    comp=[[(unit[i,j]+sign*x[i,j])/2 for j in range(start,start+n//2)] for i in range(start,start+n//2)]
    assert model["count"](comp,arb("1e-9"))==0 and (1-b)/2<arb("1e-9")


if __name__=="__main__":
    angles=[(1-arb(d)/26).acos() for d in (3,5,7)]
    f,g=angles[1]-angles[0],angles[2]-angles[0]
    blocks=[model["make"](s,True) for s in (1,-1,0)]
    tau=(arb(13)/72+arb(17).sqrt()/24).sqrt()
    assert tau<arb("0.595") and tau*arb(".42")<arb(".251")
    for n,indices in gaps.items():
        h=acb.zeta_zero(n).imag
        gap=arb(".0006") if n==13 else arb(".005")
        for base,index in zip(blocks,indices):
            k=family_matrix(base,h*f,h*g)
            assert model["count"](k,-gap)==model["count"](k,gap)==18-index
            opposite_modes(k,index,gap)
        print("positive zero",n,"three net indices",indices,"no opposite-chirality partners; gap >",gap)
    print("central index window density lower bound",((arb(".42")/arb.pi())**2).str(12))
