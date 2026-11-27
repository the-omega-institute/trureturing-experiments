"""在两个实际制备盒内检验固定见证读数的完整换基轨道及最小条件数。"""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb


def calibration():
    model = run_path(str(Path(__file__).with_name("ch39-two-zero-quark-window.py")))
    up, down = [model["prepare"](q, arb(x), arb(y)) for q,x,y in (
        (-4,"-0.680061861483031","0.873223965353125"),
        (2,"0.565247837594779","-0.630343671648397"),
    )]
    return model, up, down


def frobenius(matrix):
    return sum(abs(z)**2 for z in matrix.entries()).sqrt()


def orbit(model, up, down):
    adj = model["adj"]
    inverse_down = down.inv()
    relative = up*inverse_down
    _, frame = model["eig"](adj(relative)*relative)
    k0 = model["K"]
    inverse_k0 = k0.inv()
    reduced = adj(frame)*down*k0
    gram = reduced*adj(reduced)
    assert all(abs(gram[i,j]) > arb("1.3e-5") for i,j in ((0,1),(0,2),(1,2)))

    direct, inverse = [], []
    for j in (0,1):
        column = acb_mat([[frame[i,j]] for i in range(3)])
        direct.append((inverse_down*column)*(adj(column)*(down*k0)))
        inverse.append(((inverse_k0*inverse_down)*column)*(adj(column)*down))
    direct.append(k0-direct[0]-direct[1])
    inverse.append(inverse_k0-inverse[0]-inverse[1])
    return direct, inverse


def phases_at(x,y):
    return (acb(0,x).exp(),acb(0,y).exp(),acb(1))


def phase_sum(parts, phases):
    return sum((z*part for part,z in zip(parts,phases,strict=True)),acb_mat(3,3))


def norm_witness(matrix, adj):
    midpoint = np.array([[complex(matrix[i,j].mid()) for j in range(3)] for i in range(3)])
    left,_,right_adjoint = np.linalg.svd(midpoint)
    def unit(values):
        column = acb_mat([[acb(complex(v))] for v in values])
        return column/frobenius(column)
    return adj(unit(left[:,0])), unit(right_adjoint[0,:].conjugate())


def certify_orbit(direct, inverse, adj):
    direct_norm = [frobenius(part) for part in direct]
    inverse_norm = [frobenius(part) for part in inverse]
    remainder = sum(direct_norm[i]*inverse_norm[j]*(int(i<2)+int(j<2))**2
                    for i in range(3) for j in range(3) if i!=j)/2
    pi = arb.pi()
    cells = [(0,0,0)]
    accepted = 0
    deepest = 0
    while cells:
        depth,i,j = cells.pop()
        half = pi/2**depth
        phases = phases_at((2*i+1)*half,(2*j+1)*half)
        conjugates = [z.conjugate() for z in phases]
        left,right = norm_witness(phase_sum(direct,phases), adj)
        f = [(left*part*right)[0,0]*z for part,z in zip(direct,phases,strict=True)]
        left,right = norm_witness(phase_sum(inverse,conjugates), adj)
        g = [(left*part*right)[0,0]*z for part,z in zip(inverse,conjugates,strict=True)]
        fs,gs = sum(f),sum(g)
        slope = sum(abs((f[k]*gs-fs*g[k]).imag) for k in (0,1))
        lower = (fs*gs).real-half*slope-half**2*remainder
        if lower > 12680:
            accepted += 1
            deepest = max(deepest,depth)
        else:
            assert depth < 19,(depth,i,j,lower)
            cells.extend((depth+1,2*i+di,2*j+dj) for di in (0,1) for dj in (0,1))
    print("条件数下界、终格数、最大深度",12680,accepted,deepest)


def verify():
    model,up,down = calibration()
    adj = model["adj"]
    val = model["singular"](up*down.inv())
    assert val[0] < val[1] < val[2]
    assert arb("78.373") < val[2]/val[1] < arb("78.375")
    ref = model["ref"]
    target_ckm = acb_mat(model["load"]("shell-clock-ckm-blocks")["standard_matrix"](
        *(x.sin() for x in ref[6:9]),ref[9]))
    target = acb_mat([[ref[i]*target_ckm[i,j]/ref[j+3] for j in range(3)] for i in range(3)])
    val_target = model["singular"](target)
    assert val_target[2]/val_target[1] < arb("78.269")
    direct,inverse = orbit(model,up,down)
    phases = phases_at(arb("5.22992"),arb("4.55861"))
    upper = arb(1)
    # 先利用三个投影之和为 I，再做区间相加。
    changes = [(phases[0]-1),(phases[1]-1)]
    k0 = model["K"]
    point = k0 + phase_sum(direct[:2],changes)
    reciprocal = k0.inv() + phase_sum(inverse[:2],[z.conjugate() for z in changes])
    for matrix in (point,reciprocal):
        upper *= frobenius(adj(matrix)*matrix).sqrt()
    assert upper < arb("12681.771")
    print("指定相位上的条件数上界",upper.upper().str(16))
    certify_orbit(direct,inverse,adj)


if __name__=="__main__":
    verify()
