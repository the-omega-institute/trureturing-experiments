"""利用三重轮换和三阶舒尔补，为全双钟环面的谱计数给出证书。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


ROOT=Path(__file__).parent
cut=run_path(str(ROOT/"ch39-monotone-heavy-mass-ratio.py"))
core=[i for i in range(19) if i not in (6,7,8)]
inner=[6,7,8]


def sub(matrix,rows,cols):
    return acb_mat([[matrix[i,j] for j in cols] for i in rows])


def scalar_range(constant,a,b,c):
    # 固定 u 后，可先对 v 做精确的圆幅极值，再包住 u 的圆弧。
    for step in (1,8,64,512):
        pos=neg=True
        for k in range(step):
            theta=(arb(2*k+1)+arb(0,1))*arb.pi()/step
            u=acb(0,theta).exp()
            center=constant+2*(a*u).real
            radius=2*abs(b+c*u.conjugate())
            pos=pos and center>radius
            neg=neg and center < -radius
            if not pos and not neg:
                break
        if pos or neg:
            return int(pos)
    raise ValueError("此分隔点的全环面符号未获证书")


def schur_coefficients(charge,level):
    mat=cut["cut_matrix"](charge)
    for i in range(19):
        mat[i,i]-=level
    fixed=sub(mat,core,core)
    inverse=fixed.inv()
    base=sub(mat,inner,inner)
    W=[sub(mat,core,inner),acb_mat(16,3),acb_mat(16,3)]
    for i,j,d,u in cut["SHELL"]["edges"]:
        if d in (724,804):
            W[1 if d==724 else 2][core.index(i),inner.index(j)]=-u**charge
    count=cut["count_above"](fixed,arb(0))
    omega=acb(-arb(1)/2,arb(3).sqrt()/2)
    terms=[]
    for k in range(3):
        v=acb_mat([[1],[omega**k],[omega**(2*k)]])
        vectors=[w*v for w in W]
        def gram(i,j):
            return (vectors[i].conjugate().transpose()*inverse*vectors[j])[0,0]/3
        const=(v.conjugate().transpose()*base*v)[0,0]/3-sum(gram(i,i) for i in range(3))
        assert const.imag.contains(0)
        terms.append((const.real,-gram(0,1),-gram(0,2),-gram(1,2)))
    return count,terms


def count_below(charge,level):
    count,terms=schur_coefficients(charge,level)
    return 19-count-sum(scalar_range(*coeffs) for coeffs in terms)


def verify():
    ctx.prec=192
    for q in (1,2,3,4,6):
        levels=(".005","14.9969","14.9991","15.0009","15.0031")
        counts=[count_below(q,arb(s)) for s in levels]
        assert counts[0]==0 and counts[2]==counts[3]
        assert counts[4]-counts[1]==(2 if q==1 else 0)
        print("荷与五个全环面阈计数",q,counts)
    # 允许三辐边去刚度的范数 <1e-7 以及均值的同阶偏移。
    h_floor=arb("14.996")*arb(".00089")
    planck=arb("2.435e18")
    assert planck*h_floor>arb("3e16")
    # 高于首个新阈之后，两个荷一 Dirac 塔模至多改变 Y-W 的带符号差 56/5 倍对数。
    neutral=cut["cut_matrix"](0)
    for i,j,d,_ in cut["SHELL"]["edges"]:
        if d in (724,804):
            neutral[i,j]=neutral[j,i]=-1
    assert cut["count_above"](neutral,arb(".006"))==18
    assert cut["count_above"](neutral,arb("14.99"))==cut["count_above"](neutral,arb("15.01"))
    em,weak=arb("127.930"),arb(".23129")
    initial_y,initial_w=3*(1-weak)*em/5,weak*em
    diff=arb(41)/10+arb(19)/6
    L0=2*arb.pi()*(initial_y-initial_w)/diff
    mu0=arb("91.1876")*L0.exp()
    assert mu0<arb("3e16")
    assert diff*(arb("3e16")/mu0).log()>arb(56)/5*(arb("1e17")/arb("3e16")).log()
    alpha=1/(initial_y-(arb(41)/10+7)*L0/(2*arb.pi()))
    assert arb(".07106579")<alpha<arb(".07106580")
    print("普朗克质量单位下全谱唯一电弱交点与强耦合",mu0.str(13),alpha.str(14))


if __name__=="__main__":
    verify()
