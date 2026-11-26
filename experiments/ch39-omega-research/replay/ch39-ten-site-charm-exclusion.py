"""十站支撑的一、二维锥约束与粲顶比。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb,acb_mat,arb,ctx


def adjoint(matrix):
    return matrix.conjugate().transpose()


def singulars(column):
    matrix = acb_mat([[column[3*i+j,0] for j in range(3)] for i in range(3)])
    return sorted((v.real.max(0).sqrt() for v in (adjoint(matrix)*matrix).eig(multiple=True)),
                  key=float)


def error_norm(images, energies, start):
    scaled = acb_mat([[images[i,j]/energies[j].sqrt() for j in range(start,images.ncols())]
                      for i in range(9)])
    return max(z.real for z in (adjoint(scaled)*scaled).eig(multiple=True)).sqrt()


def column(matrix,j):
    return acb_mat([[matrix[i,j]] for i in range(matrix.nrows())])


def chart(a,b,ea,eb,tail,target,x,y,half):
    # 方框和 |z|≤1 的交内，真正允许的 z 还须满足 ea+eb|z|²≤0。
    center = acb(x,y)
    r = arb(2).sqrt()*half
    lo = (abs(center)-r).max(0)
    if lo>1:
        return
    hi = (abs(center)+r).min(1)
    weight = -ea-eb*(lo**2 if eb>0 else hi**2)
    if weight<0:
        return
    error = tail*weight.max(0).sqrt()
    sv = singulars(a+center*b)
    # 对中心到方框任一点的上型扰动用 Frobenius 范数估计。
    error += r*sum(abs(row[0])**2 for row in b.tolist()).sqrt()
    if sv[1]-error > target*(sv[2]+error):
        return
    assert half>arb(1)/2**13
    h=half/2
    for dx in (-h,h):
        for dy in (-h,h):
            chart(a,b,ea,eb,tail,target,x+dx,y+dy,h)


def check_group(k,u,target):
    values,vectors=k.eig(right=True)
    order=sorted(range(k.nrows()),key=lambda j:float(values[j].real))
    es=[values[j].real for j in order]
    assert all(a<b for a,b in zip(es,es[1:])) and es[0]<0<es[2]
    cols=[]
    for j in order:
        v=column(vectors,j)
        cols.append(v/sum(abs(row[0])**2 for row in v.tolist()).sqrt())
    images=u*acb_mat([[v[i,0] for v in cols] for i in range(k.nrows())])
    if es[1]>0:
        # Frobenius 范数先给廉价的充分界，必要时再求算子范数。
        sv=singulars(column(images,0))
        radius=(-es[0]*sum(sum(abs(images[i,j])**2 for i in range(9))/es[j]
                            for j in range(1,len(es)))).sqrt()
        if sv[1]-radius>target*(sv[2]+radius):
            return False
        radius=(-es[0]).sqrt()*error_norm(images,es,1)
        if sv[1]-radius>target*(sv[2]+radius):
            return False
    tail=error_norm(images,es,2)
    first,second=column(images,0),column(images,1)
    for a,b,ea,eb in [(first,second,es[0],es[1]),(second,first,es[1],es[0])]:
        chart(a,b,ea,eb,tail,target,arb(0),arb(0),arb(1))
    return True


def verify():
    root=Path(__file__).parent
    flavor=run_path(str(root/"ch39-minimum-flavor-support.py"))
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=192
    up,down=flavor["flavor_maps"](shell,11431)
    sites=shell["geo"]["sites"]
    rot=[sites.index((y,z,x)) for x,y,z in sites]
    ref=arb(".8573")/arb(".01410")
    k=adjoint(down)*down-3*adjoint(up)*up/ref**2
    total=charts=0
    for group in combinations(range(19),10):
        orbit=flavor["orbit"](group,rot)
        if group!=min(orbit):
            continue
        us=acb_mat([[up[i,j] for j in group] for i in range(9)])
        ks=flavor["submatrix"](k,group)
        charts+=check_group(ks,us,arb(".01"))
        total+=1
    assert total==30806 and charts==19
    print("十站支撑循环类、需二维锥的循环类、严格粲顶比下界",total,charts,"0.01")


if __name__=="__main__":
    verify()
