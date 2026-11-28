"""由相对奇异标架以及 S、S^-1 的实性检验实壳源能否命中夸克观测盒。"""
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from flint import acb, acb_mat, arb


SCALE = 2**24


def encode(value):
    low = int((value.lower()*SCALE).floor().unique_fmpz())
    high = int((value.upper()*SCALE).ceil().unique_fmpz())
    mid = (low+high)//2
    return mid, max(mid-low,high-mid)


def phase_table(count):
    angles = [(2*arb.pi()*i/count).union(2*arb.pi()*(i+1)/count) for i in range(count)]
    cosine = np.array([encode(v.cos()) for v in angles],dtype=np.int64)
    sine = np.array([encode(v.sin()) for v in angles],dtype=np.int64)
    differences = [x-y for x in angles for y in angles]
    dc = np.array([encode(v.cos()) for v in differences],dtype=np.int64)
    ds = np.array([encode(v.sin()) for v in differences],dtype=np.int64)
    table = [np.array([np.full(count*count,(SCALE,0)[k]), dc[:,k], ds[:,k],
                       np.repeat(cosine[:,k],count),np.repeat(sine[:,k],count),
                       np.tile(cosine[:,k],count),np.tile(sine[:,k],count)],dtype=np.int64)
             for k in (0,1)]
    assert all(np.max(np.abs(part)) < 2*SCALE for part in table)
    return table


def surviving_cells(rows, table, indices):
    encoded = np.array([[encode(v) for v in row] for row in rows],dtype=np.int64)
    assert np.max(np.abs(encoded)) < 4*SCALE
    assert 7*4*6*SCALE**2 < 2**63
    coefficients, errors = encoded[:,:,0],encoded[:,:,1]
    phase, radius = [part[:,indices] for part in table]
    lower = (np.abs(coefficients@phase)-np.abs(coefficients)@radius
             -errors@(np.abs(phase)+radius))
    return indices[np.max(lower,axis=0)<=0]


def gram(rows):
    return acb_mat([[sum(row[i]*row[j].conjugate() for row in rows)
                     for j in range(3)] for i in range(3)])


def unit(vector):
    length = sum(abs(z)**2 for z in vector).sqrt()
    return [z/length for z in vector]


def target_frame(values, eigenvalues, standard):
    up, down = values[:3],values[3:6]
    mixing = acb_mat(standard(*(x.sin() for x in values[6:9]),values[9]))
    direct, inverse = [], []
    for i in range(3):
        # 分别以直接矩阵元 (2,0) 和逆矩阵的主列元 (1,0) 归一化。
        direct_row, inverse_row = [], []
        for j in range(3):
            entry = mixing[i,j].conjugate()
            q = acb(1) if (i,j)==(2,0) else entry/mixing[2,0].conjugate()
            direct_row.append(q*(up[i]/up[2] if i!=2 else 1)*(down[0]/down[j] if j!=0 else 1))
            q = acb(1) if (i,j)==(0,1) else entry/mixing[0,1].conjugate()
            inverse_row.append(q*(up[0]/up[i] if i!=0 else 1)*(down[j]/down[1] if j!=1 else 1))
        direct.append(direct_row)
        inverse.append(inverse_row)
    h, k = gram(direct), gram(inverse)
    nh, nk = [sum(sum(abs(z)**2 for z in row) for row in group) for group in (direct,inverse)]

    def largest_column(matrix, opposite, norm, opposite_norm, levels, target, pivot):
        b,c = [v for i,v in enumerate(levels) if i!=target]
        trace = sum(levels)
        opposite_trace = sum(1/v for v in levels)
        u = [matrix[i,pivot]-(b+c)*norm/trace*int(i==pivot)
             +(b*c/trace*opposite_trace)*(norm/opposite_norm)*opposite[i,pivot] for i in range(3)]
        assert not u[pivot].real.contains(0)
        return unit([u[i]/u[pivot].real if i!=pivot else acb(1) for i in range(3)])

    low = largest_column(k,h,nk,nh,[1/x for x in eigenvalues],0,1)
    top = largest_column(h,k,nh,nk,eigenvalues,2,0)
    middle = [(top[(i+1)%3]*low[(i+2)%3]-top[(i+2)%3]*low[(i+1)%3]).conjugate()
              for i in range(3)]
    return acb_mat([low,middle,top]).transpose()


def axes(matrix):
    # 实左乘只组合实性约束。浮点只用于选择随后固定为二进有理数的实矩阵。
    mid = np.array([[complex(z.mid()) for z in row] for row in matrix.tolist()])
    eigen, basis = np.linalg.eigh((mid@mid.conj().T).real)
    transform = acb_mat(((basis/np.sqrt(eigen)).T).tolist())
    matrix = transform*matrix
    result = []
    for i,j in ((0,1),(0,2),(1,2)):
        result.append([[(matrix[i,a]*matrix[j,b].conjugate()-matrix[j,a]*matrix[i,b].conjugate())
                        /acb(0,2) for b in range(3)] for a in range(3)])
    return result


def coefficients(values, eigenvalues, standard, normals, adjoint):
    frame = target_frame(values,eigenvalues,standard)
    down = values[3:6]
    rows = []
    for inverse in (False,True):
        weights = [(down[0]/v)**2 if inverse else (v/down[2])**2 for v in down]
        weights[0 if inverse else 2] = arb(1)
        diagonal = acb_mat([[weights[i] if i==j else 0 for j in range(3)] for i in range(3)])
        moment = adjoint(frame)*diagonal*frame
        for normal in normals[int(inverse)]:
            row = [sum(normal[a][a].real*moment[a,a].real for a in range(3))]
            for a,b in ((0,1),(0,2),(1,2)):
                z = normal[a][b]*moment[a,b]
                row.extend((2*z.real,-2*z.imag))
            rows.append(row)
    return rows


def certify_phase_box(rows, depth, i, j):
    half = arb.pi()/2**depth
    bump = arb(0,half.upper())
    x, y = (2*i+1)*half+bump, (2*j+1)*half+bump
    phases = (1,(x-y).cos(),(x-y).sin(),x.cos(),x.sin(),y.cos(),y.sin())
    if any(not sum(a*b for a,b in zip(row,phases,strict=True)).contains(0) for row in rows):
        return True
    if depth>=11:
        return False
    if depth==7:
        center_x, center_y = x.mid(), y.mid()
        difference = center_x-center_y
        center = (1,difference.cos(),difference.sin(),center_x.cos(),center_x.sin(),
                  center_y.cos(),center_y.sin())
        if all(sum(a*b for a,b in zip(row,center,strict=True)).contains(0) for row in rows):
            return False
    return all(certify_phase_box(rows,depth+1,2*i+a,2*j+b) for a in (0,1) for b in (0,1))


def verify(replay):
    orbit = run_path(str(replay/"ch39-two-zero-orbit-condition.py"))
    model, up, down = orbit["calibration"]()
    shell = run_path(str(replay/"ch39-shell-pencil-six-mass-window.py"))
    real_frame, _ = shell["shell_frame"](model)
    adj = model["adj"]
    relative = up*down.inv()
    eigenvalues, basis = model["eig"](adj(relative)*relative)
    # 行列式公式收紧最小本征值区间，保留全制备盒。
    eigenvalues[0] = abs(relative.det())**2/(eigenvalues[1]*eigenvalues[2])
    assert 0 < eigenvalues[0] < eigenvalues[1] < eigenvalues[2]
    normals = [axes(adj(real_frame)*v*basis) for v in (down.inv(),adj(down))]
    standard = model["load"]("shell-clock-ckm-blocks")["standard_matrix"]
    phase_count = 128
    table = phase_table(phase_count)
    priority = [1,10,5,1100,550,560,260,430,990,1010]
    stack = [([0.]*10,[1.]*10,0,np.arange(phase_count**2))]
    leaves = deepest = 0
    while stack:
        center, size, depth, cells = stack.pop()
        values = [ref+arb(c)*w+arb(0,(arb(h)*w).upper()) for ref,w,c,h in
                  zip(model["ref"],model["tol"],center,size,strict=True)]
        rows = coefficients(values,eigenvalues,standard,normals,adj)
        cells = surviving_cells(rows,table,cells)
        if 0 < len(cells) <= 30 and depth>=7:
            cells = np.array([ij for ij in cells if not certify_phase_box(rows,7,int(ij)//phase_count,int(ij)%phase_count)], dtype=np.int64)
        if not len(cells):
            leaves += 1
            deepest = max(deepest,depth)
            continue
        assert depth<36,(center,size,len(cells))
        j = max(range(10),key=lambda i:priority[i]*size[i])
        smaller = size.copy()
        smaller[j] /= 2
        for sign in (-1,1):
            shifted = center.copy()
            shifted[j] += sign*smaller[j]
            stack.append((shifted,smaller,depth+1,cells))
    print("完整十窗的终结观测盒数与最大二分深度",leaves,deepest)


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
