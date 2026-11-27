"""保持 N+2 完整 CKM 见证的上型四站目标，在原三平面上恢复 N+3 的轻上质量。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat, ctx


# N+3 按能量升序列出各目标线所取的 N+2 源编号。
SOURCE_AT_DESTINATION = (16, 6, 13, 7, 5, 12, 4, 8, 14, 10, 17, 3, 19, 11, 2, 1, 9, 18, 15)


def adj(matrix):
    return matrix.conjugate().transpose()


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    embedding = load("common-embedding-central-right-ckm")
    independent = load("independent-unitary-quark-right-holdout")
    full = load("one-generator-gravity-full-ckm")
    node = independent["CENTER"]+2
    _,labels,_ = embedding["verify_nodes"](
        root,((node-1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    original = labels[node-1,-4]
    target = full["targets"](root,independent,embedding)[0]
    rows = load("free-four-site-native-positive-spectrum")["native"](root,node)[0]
    coords = independent["coordinates"]
    cols = arb_mat([[*coords(row),arb(1)] for row in rows]).transpose()
    rhs = arb_mat([[v] for v in target.entries()+[sum(original)]])
    seed = independent["SEEDS"][0]
    unknown = [i for i,num in enumerate(seed) if num is None]
    known = arb_mat([[arb(num or 0)/10**12] for num in seed])
    system = arb_mat([[cols[i,j] for j in unknown] for i in range(7)])
    assert not system.det().contains(0)
    solution = system.solve(rhs-cols*known,algorithm="precond")
    for i,val in zip(unknown,solution.entries(),strict=True):
        known[i,0] = val
    y = known.entries()
    diff = cols*known-rhs
    assert all(val.contains(0) and abs(val)<arb("1e-38") for val in diff.entries())
    assert min(y)>0
    before = sorted((v.lower() for v in original),reverse=True)
    after = sorted((v.upper() for v in y),reverse=True)
    gap = min((sum(before[:k])-sum(after[:k])).lower() for k in range(1,19))
    assert gap>0
    print("第一步七个未知编号、十八个严格主序余量的最小下界",
          [i+1 for i in unknown],gap.str(12),flush=True)

    shell = load("three-ray-chiral-boundary")
    sites = (shell["soft"]["pole"],*shell["boundary_sites"]())
    assert sites==(14,8,7,6)
    end = load("adjacent-pvm-transport")["spectrum"](shell,-4,node+1)[1]
    s,t,u = [1/arb(i).sqrt() for i in (3,2,6)]
    four = acb_mat([[1,0,0,0],[0,s,t,u],[0,s,-t,u],[0,s,0,-2*u]])
    _,_,intrinsic,rotation = embedding["common_embedding"](*independent["GEOMETRY"])
    frame = four*rotation*intrinsic
    rays = adj(frame)*acb_mat([[end[i,j] for j in range(19)] for i in sites])
    weights = [y[j-1] for j in SOURCE_AT_DESTINATION]
    assert sorted(SOURCE_AT_DESTINATION)==list(range(1,20))
    weighted = load("common-unitary-quark-right-trace")["weighted"]
    center = arb(".00000615")

    def levels(prob):
        mixed = weights.copy()
        delta = prob*(weights[9]-weights[3])
        mixed[3] += delta
        mixed[9] -= delta
        block = weighted(rays,mixed)
        eig = block.eig(algorithm="rump")
        assert all(val.imag.contains(0) for val in eig)
        result = sorted((val.real for val in eig),key=lambda val:val.mid())
        assert all(a<b for a,b in zip(result,result[1:]))
        return result

    straight = levels(arb(0))
    assert arb(".00000621799") < straight[0] < arb(".00000621802")
    print("纯置换在原三平面上的三个上型振幅",*[v.str(14) for v in straight],flush=True)
    lo,hi = arb(".31278"),arb(".31279")
    left,right = levels(lo),levels(hi)
    assert left[0]>center>right[0]
    print("两线旋转的概率、轻上振幅与中心之差",
          lo.str(7),(left[0]-center).str(11),hi.str(7),(right[0]-center).str(11),flush=True)
    # 同时控制其余两个质量。两端矩阵的差由一对确认归一化的目的站投影给出，
    # 压缩不增加算子范数，故在整个概率小段上用 Weyl 界即可。
    perturbation = (hi-lo)*abs(weights[9]-weights[3])
    rest = [val+arb(0,perturbation.upper()) for val in left[1:]]
    assert arb(".0002071")<rest[0]<arb(".0002072") and arb(".056149")<rest[1]<arb(".056152")
    print("精确轻上见证的粲、顶包围",*[v.str(11) for v in rest],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
