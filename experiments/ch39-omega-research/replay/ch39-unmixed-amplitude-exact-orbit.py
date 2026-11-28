"""在两个实际零点的完整十项读数轨道上排除单部门谱函数振幅的平方和。"""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb


def load(name):
    return run_path(str(Path(__file__).with_name(f"ch39-{name}.py")))


def gram(vectors,adj):
    return acb_mat([[(adj(v)*w).trace() for w in vectors] for v in vectors])


def orthonormal(words,adj,norm):
    basis = []
    for w in words:
        for previous in basis:
            w -= previous*(adj(previous)*w).trace()
        length = norm(w)
        assert length > 0
        basis.append(w/length)
    return basis


def certify(model,orbit,up,down):
    adj,norm = model["adj"],orbit["frobenius"]
    u,d = adj(up)*up,adj(down)*down
    identity = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    basis = orthonormal((identity,u,u*u,d,d*d),adj,norm)
    original_gram = gram(basis,adj)
    numerical_basis = np.array([w.tolist() for w in basis],dtype=complex)
    parts,_ = orbit["orbit"](model,up,down)
    # 依次对应频率 (1,-1)、(1,0)、(0,1)。
    pairs = ((0,1),(0,2),(1,2))
    cross = [parts[i]*adj(parts[j]) for i,j in pairs]

    def separates(x,y,half):
        phases = orbit["phases_at"](x,y)
        center = orbit["phase_sum"](parts,phases)
        # 浮点计算仅挑选候选直线和阈值，此后两侧严格不等式都用球算术复核。
        left = np.linalg.svd(np.array(center.tolist(),dtype=complex))[0][:,0]
        approximate_projection = np.eye(3)-np.outer(left,left.conj())
        projected = (approximate_projection@numerical_basis).reshape(5,9)
        threshold = arb(float(.998*np.linalg.svd(projected,compute_uv=False)[-1]**2))
        vector = acb_mat([[acb(complex(entry))] for entry in left])
        projection = identity-vector*adj(vector)/(adj(vector)*vector)[0,0]

        functional = identity*threshold-projection
        value = threshold*norm(center)**2 - norm(projection*center)**2
        harmonics = [(term*functional).trace()*phases[i]*phases[j].conjugate()
                     for term,(i,j) in zip(cross,pairs,strict=True)]
        mixed,x_only,y_only = harmonics
        slope = 2*(abs((mixed+x_only).imag)+abs((y_only-mixed).imag))
        curvature = 4*abs(mixed)+abs(x_only)+abs(y_only)
        if not value > half*slope+half**2*curvature:
            return False
        pencil = gram([projection*w for w in basis],adj)-original_gram*threshold
        for size in range(1,6):
            minor = acb_mat([[pencil[i,j] for j in range(size)] for i in range(size)])
            if not minor.det().real > 0:
                return False
        return True

    cells = [(0,0,0)]
    leaves = deepest = 0
    while cells:
        depth,i,j = cells.pop()
        half = arb.pi()/2**depth
        if separates((2*i+1)*half,(2*j+1)*half,half):
            leaves += 1
            deepest = max(deepest,depth)
        else:
            assert depth < 14, (depth,i,j)
            cells.extend((depth+1,2*i+di,2*j+dj) for di in (0,1) for dj in (0,1))
    print("完整固定读数轨道的分离终格数与最大深度",leaves,deepest)


if __name__ == "__main__":
    orbit = load("two-zero-orbit-condition")
    model,_,_ = orbit["calibration"]()
    gram_source = load("quadratic-gram-source")
    up,down = [gram_source["at_zero"](model,q,n) for q,n in ((-4,2345666),(2,964890))]
    certify(model,orbit,up,down)
