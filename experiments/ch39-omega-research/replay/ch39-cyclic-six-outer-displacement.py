"""固定六个外站，内十三站允许任意三循环等变酉输运。"""
from itertools import product
from pathlib import Path
from runpy import run_path
from flint import acb_mat, arb, ctx

def adj(x):
    return x.conjugate().transpose()

def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/('ch39-'+name+'.py')))
    ctx.prec=256
    graph=lib('three-ray-chiral-boundary')
    sites=graph['geo']['sites']
    outer=[i for i,s in enumerate(sites) if sorted(s)==[1,5,29]]
    inner=[i for i in range(19) if i not in outer]
    assert outer==list(range(6))
    # 以三循环的平凡字符正交投影替换源的第一列。
    rot=[sites.index((y,z,x)) for x,y,z in sites]
    P=acb_mat([[arb(sum(j==k for k in (i,rot[i],rot[rot[i]])))/3
                for j in range(19)] for i in range(19)])
    source=fifth*columns
    col=acb_mat([[source[i,0]] for i in range(19)])
    retained=(P*col)
    src_in=acb_mat([[retained[i,0]] for i in inner])
    src_out=acb_mat([[retained[i,0]] for i in outer])
    ins=lib('three-shell-rank-one-obstruction')
    # 该检验同时从三个壳的完整弱核导出候选与投影上界。
    ins["experiment"](root,fifth,columns,levels)
    bases=[adj(lib('adjacent-pvm-transport')['spectrum'](graph,6,n)[1])
           for n in range(2_306_888_100_092, 2_306_888_100_095)]
    norm=lambda m:sum(abs(y).upper()**2 for y in m.entries()).sqrt()
    survivors=lib('two-block-low-line-gram')['SURVIVORS']
    bounds=[]
    displacements={}
    for ids in product(*ins['CANDIDATES']):
        if ids in survivors:
            continue
        S=acb_mat([[s[k,i] for i in range(19)] for s,k in zip(bases,ids,strict=True)])
        # 谱行与其平凡字符投影的差直接纳入上界，无需先假设它严格为零。
        e=norm(S-S*P)
        SP=S*P
        A=acb_mat([[SP[i,j] for j in inner] for i in range(3)])
        L=acb_mat([[SP[i,j] for j in outer] for i in range(3)])
        B=L*src_out
        dual=adj(A)*(A*adj(A)).inv()
        # ||dual||_op<=||dual||_F，且 dual*A 是正交投影。
        H=dual*B
        q=(sum(abs(v)**2 for v in H.entries()).sqrt()-norm(src_in))/norm(dual)-e*norm(col)
        bounds.append(q.lower())
        excess=sum(abs(v)**2 for v in H.entries()).sqrt()-norm(src_in)
        # δ=||(U-I)P_out||=||P_out(U-I)||；酉性给出内、外两向同一偏差上界。
        gap=(excess-norm(dual)*(arb("2.853e-8")+e*norm(col)))/(norm(dual*L)*norm(retained)+norm(src_out))
        displacements[ids]=gap.lower()
    assert len(bounds)==16
    print("十六种低线的第一重列投影范数统一下界",min(bounds).str(14),flush=True)
    print("十六种低线的外站变动下界、三循环九式见证选线的下界",
          min(displacements.values()).str(14), displacements[(0,0,15)].str(14), flush=True)
    assert min(bounds)>arb('.0064')
    assert min(displacements.values())>arb('.0097') and displacements[(0,0,15)]>arb('.5094')


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    src = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root, *src["verify"](root, return_operator="spectral"))
