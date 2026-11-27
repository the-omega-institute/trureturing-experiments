"""两枚实际零点的五荷原生谱线恰可辨识全部三个字符对角块。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def norm_infinity(matrix):
    return max(sum(abs(z).upper() for z in row) for row in matrix.tolist())


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    graph = lib("three-ray-chiral-boundary")
    frames, rotation = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    lib("eleven-weak-eight-sites")["check_cyclic_geometry"](graph,rotation)
    ctx.prec = 384
    heights = [acb.zeta_zero(2_306_888_100_085+j).imag for j in (7,8)]
    deltas = [arb(23)/25*arb(2)**-20/(h**2+graph["soft"]["point"]["study"]["offset"]) for h in heights]
    difference = deltas[1]-deltas[0]
    assert difference < 0
    sensitivity = abs(difference)/((1-deltas[0])**2+5).sqrt()
    assert sensitivity < arb("6.51e-43")
    print("两壳星边亏额差",difference.str(13),flush=True)
    design = [[],[],[]]
    for hindex,height in enumerate(heights):
        for charge in (1,-4,2,-3,6):
            lap = graph["laplacian"](charge,height)
            levels=[]
            for k,frame in enumerate(frames):
                block = frame.conjugate().transpose()*lap*frame
                if k==0:
                    for j in range(6):
                        w=1-deltas[hindex] if j==0 else arb(1)
                        assert (block[6,j]+arb(3).sqrt()*w).contains(0)
                values,vectors = block.eig(right=True,algorithm="rump")
                assert all(z.imag.contains(0) for z in values)
                levels.extend(z.real for z in values)
                size = block.ncols()
                for j in range(size):
                    norm = sum(abs(vectors[i,j])**2 for i in range(size)).sqrt()
                    assert norm>0
                    v=[vectors[i,j]/norm for i in range(size)]
                    row=[abs(z)**2 for z in v]
                    for a,b in combinations(range(size),2):
                        term=v[a]*v[b].conjugate()
                        row.append(term.real)
                        if k==0 and (a,b)==(5,6):
                            # 用极点本征方程作列变换，精确除去极小的两壳星边亏额差。
                            row.append(arb(0) if hindex==0 else (v[0]*v[6].conjugate()).imag)
                        else:
                            row.append(term.imag)
                    design[k].append(row)
            levels.sort(key=lambda x:x.mid())
            assert all(x<y for x,y in zip(levels,levels[1:]))
    for k,rows in enumerate(design):
        mat=arb_mat(rows)
        gram=mat.transpose()*mat
        size=mat.ncols()
        assert size==(49,36,36)[k]
        pre=gram.mid().inv()
        identity=arb_mat([[int(i==j) for j in range(size)] for i in range(size)])
        defect=norm_infinity(identity-pre*gram)
        print("字符、读数行数、列秩、预条件缺陷",k,mat.nrows(),size,defect.str(12),flush=True)
        assert defect.is_finite() and defect<arb(".001")


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
