"""三条非负边权共用的固定低能二平面，按第三原边及五种荷认证。"""
from functools import cache
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from flint import acb, acb_mat, arb, ctx


FIXED = ((0, 2), (5, 14))
INNER = (6, 7, 8)
OUTER = tuple(i for i in range(19) if i not in INNER)
CHARGES = (1, -4, 2, -3, 6)


def adj(a):
    return a.conjugate().transpose()


def sub(a, rows, cols):
    return acb_mat([[a[i,j] for j in cols] for i in rows])


def numeric(a):
    return np.array(a.tolist(), complex)


def outside_constraints(base, edges):
    vectors = []
    for i, j in edges:
        if i in INNER and j in INNER:
            continue
        v = [acb(0)]*len(OUTER)
        if i in OUTER and j in OUTER:
            v[OUTER.index(i)], v[OUTER.index(j)] = 1, base[j,i]
        else:
            v[OUTER.index(i if i in OUTER else j)] = 1
        vectors.append(v)
    return acb_mat(vectors).transpose()


def propose(base, affine, edges):
    # 浮点只提出一个有理阈值。两个接受判据在下文独立用球算术核对。
    c = numeric(outside_constraints(base, edges))
    _, _, vh = np.linalg.svd(c.conj().T, full_matrices=True)
    basis = vh.conj().T[:, np.linalg.matrix_rank(c):]
    outside = numeric(sub(base, OUTER, OUTER))
    alpha = np.linalg.eigvalsh(basis.conj().T @ outside @ basis)[1]
    cuts = np.zeros((len(affine), 19, 3), complex)
    for k, (i,j) in enumerate(edges):
        cuts[:,i,k] = 1
        cuts[:,j,k] = affine[:,j,i]
    beta = np.linalg.eigvalsh(affine-cuts @ cuts.conj().transpose(0,2,1))[:,3].min()
    if beta <= alpha + .01:
        return None
    return int(np.ceil(10000*alpha)) + 2


def plane_below(base, edges, threshold):
    c = outside_constraints(base, edges)
    outside = sub(base, OUTER, OUTER)
    n, r = outside.nrows(), c.ncols()
    bordered = acb_mat(n+r, n+r)
    for i in range(n):
        for j in range(n):
            bordered[i,j] = outside[i,j] - (threshold if i == j else 0)
        for j in range(r):
            bordered[i,n+j] = c[i,j]
            bordered[n+j,i] = c[i,j].conjugate()
    eig = bordered.eig(multiple=True)
    return sum(z.real < 0 for z in eig) >= r+2


@cache
def phase(n, i):
    return acb(0, arb.pi()*(2*i+1)/n).exp()


@cache
def change(n):
    return 4*arb(2).sqrt()*(arb.pi()/(2*n)).sin()+6*arb("1e-8")


def positive_count(a, needed):
    if needed <= 0:
        return True
    if needed > 3:
        return False
    # 首先用小主子式证明正定子空间；余下情形由认证的三阶谱决定。
    for i,j in ((0,1),(0,2),(1,2)):
        if a[i,i].real > 0:
            if needed == 1:
                return True
            if a[i,i].real*a[j,j].real - abs(a[i,j])**2 > 0:
                if needed == 2 or a.det().real > 0:
                    return True
    return sum(z.real > 0 for z in a.eig(multiple=True)) >= needed


def deleted_above(base, links, edges, threshold, clock_excess=0):
    deleted = [acb_mat(m) for m in (base, *links)]
    for i,j in edges:
        deleted[0][i,i] -= 1
        deleted[0][j,j] -= 1
        for m in deleted:
            m[i,j] = m[j,i] = 0
    outside = sub(deleted[0], OUTER, OUTER)
    inner = sub(deleted[0], INNER, INNER)
    terms = [sub(m, OUTER, INNER) for m in deleted]
    outer_eig = [v.real for v in outside.eig(multiple=True)]
    @cache
    def polynomial(n):
        cut = threshold + change(n) + 2*clock_excess*(arb.pi()/(2*n)).sin()
        large = sum(v > cut for v in outer_eig)
        assert all(v < cut or v > cut for v in outer_eig)
        matrix = acb_mat(outside)
        for i in range(matrix.nrows()):
            matrix[i,i] -= cut
        r = matrix.inv()
        center = acb_mat(inner)
        for i in range(3):
            center[i,i] -= cut
        for m in terms:
            center -= adj(m)*r*m
        coeff = [-adj(terms[i])*r*terms[j] for i,j in ((0,1),(0,2),(1,2))]
        return 16-large, center, coeff

    todo = [(8,i,j) for i in range(8) for j in range(8)]
    leaves = 0
    finest = 8
    while todo:
        n,i,j = todo.pop()
        needed, center, (x,y,z) = polynomial(n)
        u,v = phase(n,i),phase(n,j)
        extra = u*x + v*y + u.conjugate()*v*z
        if positive_count(center+extra+adj(extra), needed):
            leaves += 1
            finest = max(finest,n)
        elif n >= 2048:
            return None
        else:
            todo.extend((2*n,2*i+a,2*j+b) for a in range(2) for b in range(2))
    return leaves, finest


WEIGHTS = np.r_[np.arange(0,2.01,.125),np.arange(2.25,4.01,.25),
                np.arange(4.5,8.01,.5),np.arange(9,17),24,32,48,64,128,256,1024]


def set_weight(base, links, edge, weight):
    i,j = edge
    updated = [acb_mat(m) for m in (base,*links)]
    updated[0][i,i] += weight-1
    updated[0][j,j] += weight-1
    for m in updated:
        m[i,j] *= weight
        m[j,i] *= weight
    return updated[0], updated[1:]


def proposed_bands(base, affine, triple):
    columns = np.zeros((len(affine),19,3), complex)
    for k,(i,j) in enumerate(triple):
        columns[:,i,k] = 1
        columns[:,j,k] = affine[:,j,i]
    outside = numeric(sub(base, OUTER, OUTER))
    plans = []
    for k, moving in enumerate(triple):
        fixed = tuple(e for e in triple if e!=moving)
        fixed_cols = columns[:,:,[l for l in range(3) if l!=k]]
        bottom = affine - fixed_cols @ fixed_cols.conj().transpose(0,2,1)
        v = columns[:,:,[k]]
        ray = v@v.conj().transpose(0,2,1)
        c = numeric(outside_constraints(base,fixed))
        _,_,vh = np.linalg.svd(c.conj().T,full_matrices=True)
        basis = vh.conj().T[:,np.linalg.matrix_rank(c):]
        a = basis.conj().T @ outside @ basis
        h = basis.conj().T @ v[0,OUTER,:]
        rank_one = h@h.conj().T
        alpha = [np.linalg.eigvalsh(a+(w-1)*rank_one)[1] for w in WEIGHTS]
        if np.linalg.norm(h)>1e-10:
            tail = np.linalg.qr(h,mode="complete")[0][:,1:]
            alpha.append(np.linalg.eigvalsh(tail.conj().T@a@tail)[1])
        else:
            alpha.append(alpha[0])
        beta = [np.linalg.eigvalsh(bottom+(w-1)*ray)[:,3].min() for w in WEIGHTS]
        quality = min(np.array(beta)-alpha[:-1])
        if quality < .003:
            continue
        bands = []
        left = 0
        while left<len(WEIGHTS):
            valid = [j for j in range(left+1,len(alpha)) if beta[left]-alpha[j]>.45*quality]
            if not valid:
                break
            right = max(valid)
            bands.append((left,right,int(np.ceil(10000*alpha[right]))+2))
            left = right
        if left==len(WEIGHTS):
            clocked = abs(np.ptp(v[:,moving[1],0]))>1e-8
            cost = sum((1+(max(0,WEIGHTS[i]-1) if clocked else 0))**2
                       /(beta[i]-alpha[j])**2 for i,j,_ in bands)
            plans.append((cost,fixed,moving,bands))
    return sorted(plans)


def verify_bands(base,links,affine,triple):
    plans = proposed_bands(base,affine,triple)
    for _,fixed,moving,bands in plans:
        output=[]
        for i,j,numerator in bands:
            lo=float(WEIGHTS[i])
            high=float(WEIGHTS[j]) if j<len(WEIGHTS) else None
            threshold=arb(numerator)/10000
            if high is None:
                good_plane=plane_below(base,triple,threshold)
            else:
                top,_=set_weight(base,links,moving,arb(high))
                good_plane=plane_below(top,fixed,threshold)
            if not good_plane:
                break
            low,mats=set_weight(base,links,moving,arb(lo))
            extra=max(0,lo-1) if any(m[moving]!=0 for m in links) else 0
            certificate=deleted_above(low,mats,fixed,threshold,extra)
            if certificate is None:
                break
            output.append((lo,high,numerator,*certificate))
        if len(output)==len(bands):
            return moving,output
    return None


def verify(root):
    model = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    factory = run_path(str(root / "ch39-native-light-projector-pump.py"))
    ctx.prec = 112
    assert set(INNER) == {i for i,p in enumerate(model["geo"]["sites"]) if sorted(p)==[13,13,23]}
    edges = [(i,j) for i,j,_,_ in model["geo"]["edges"] if (i,j) not in FIXED]
    assert len(edges) == len(set(edges)) == 133
    records = {e: [] for e in edges}
    for q in CHARGES:
        base, links = factory["clock_matrices"](model,q)
        assert all(m[i,j] == 0 for m in links for i in range(19) for j in range(19)
                   if i not in OUTER or j not in INNER)
        nodes = np.exp(2j*np.pi*np.arange(8)/8)
        d,x,y = (numeric(m) for m in (base,*links))
        affine = []
        for u in nodes:
            for v in nodes:
                off = u*x+v*y
                affine.append(d+off+off.conj().T)
        affine = np.stack(affine)
        assert plane_below(base,FIXED,arb("8.178"))
        axis = deleted_above(base,links,FIXED,arb("8.55"))
        assert axis is not None
        print("第三增量零轴",q,*axis,flush=True)
        for e in edges:
            triple = (*FIXED,e)
            candidate = propose(base,affine,triple)
            if candidate is None:
                continue
            threshold = arb(candidate)/10000
            if not plane_below(base,triple,threshold):
                continue
            certificate = deleted_above(base,links,triple,threshold)
            if certificate:
                leaves,finest=certificate
                records[e].append(q)
                print("认证",q,*e,candidate,leaves,finest,flush=True)
        print("逐荷已认证",q,sum(q in v for v in records.values()),flush=True)
        for e in edges:
            if q in records[e]:
                continue
            result=verify_bands(base,links,affine,(*FIXED,e))
            if result is not None:
                records[e].append(q)
                print("分段认证",q,*e,*result,flush=True)
            else:
                print("未覆盖",q,*e,flush=True)
    complete = [e for e,c in records.items() if len(c)==len(CHARGES)]
    print("五荷均认证",len(complete),*complete)
    print("其余原边",len(edges)-len(complete),*(e for e in edges if e not in complete))
    assert len(complete) == 133, "仍有第三边或荷尚未通过认证"


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
