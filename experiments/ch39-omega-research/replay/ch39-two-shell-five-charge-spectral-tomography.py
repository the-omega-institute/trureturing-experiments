"""两枚实际零点的五荷原生谱线恰可辨识全部三个字符对角块。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def norm_infinity(matrix):
    return max(sum(abs(z).upper() for z in row) for row in matrix.tolist())


def verify(root, *, first_zero=False):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    graph = lib("three-ray-chiral-boundary")
    frames, rotation = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    lib("eleven-weak-eight-sites")["check_cyclic_geometry"](graph,rotation)
    ctx.prec = 384
    numbers = (1,2_306_888_100_092) if first_zero else (2_306_888_100_092,2_306_888_100_093)
    heights = [acb.zeta_zero(j).imag for j in numbers]
    deltas = [arb(23)/25*arb(2)**-20/(h**2+graph["soft"]["point"]["study"]["offset"]) for h in heights]
    difference = deltas[1]-deltas[0]
    assert difference < 0
    sensitivity = abs(difference)/((1-deltas[0])**2+5).sqrt()
    if not first_zero:
        assert sensitivity < arb("6.51e-43")
    if first_zero:
        assert abs(difference)>arb("3.8118e-9")
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
    def rank_bound(rows,size):
        mat=arb_mat(rows)
        gram=mat.transpose()*mat
        assert mat.ncols()==size
        pre=gram.mid().inv()
        identity=arb_mat([[int(i==j) for j in range(size)] for i in range(size)])
        defect=norm_infinity(identity-pre*gram)
        assert defect.is_finite() and defect<1
        return defect,(norm_infinity(pre*mat.transpose())/(1-defect)).upper()

    for k,rows in enumerate(design):
        size=(49,36,36)[k]
        defect,error=rank_bound(rows,size)
        # 新原生谱线在同一列变换下有界，直接外推绕开逆列变换的放大。
        if first_zero:
            dimension=(7,6,6)[k]
            drift=deltas[0]/abs(difference) if k==0 else arb(0)
            target_norm=1+arb(dimension-1)/arb(2).sqrt()+drift/2
            prediction=(target_norm*error).upper()
            print("任意后续简单谱行的误差放大上界",k,prediction.str(12),flush=True)
            assert prediction < arb("2e5")
        transform = 1+1/abs(difference) if k==0 else arb(1)
        error=(error*transform).upper()
        print("索引、字符、原坐标逐项误差的放大上界",numbers,k,error.str(12),flush=True)
        if first_zero:
            assert error < arb(("8.961e12","2843","10570")[k])
        print("字符、读数行数、列秩、预条件缺陷",k,len(rows),size,defect.str(12),flush=True)
        assert defect.is_finite() and defect<arb(".001")

    # 每个设置在七维块只有六个独立于公共迹的额外读数。检查任取八个设置。
    worst=arb(0)
    for omit in combinations(range(10),2):
        for rows,dim in zip(design,(7,6,6),strict=True):
            chosen=[row for i,row in enumerate(rows) if i//dim not in omit]
            defect,_=rank_bound(chosen,dim**2)
            assert defect.is_finite() and defect<arb(".001"),(omit,dim)
            worst=max(worst,defect.upper())
    print("任取八种荷与高度组合，45 次删选中的最大预条件缺陷",worst.str(12),flush=True)


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    verify(root)
    verify(root,first_zero=True)
