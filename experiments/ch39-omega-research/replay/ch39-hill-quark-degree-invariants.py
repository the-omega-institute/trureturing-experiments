"""在任意共同味基下，三荷同一谱多项式的最低次数约束。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, arb_mat, ctx


def moment_system(root):
    ray = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    common = run_path(str(root/"ch39-signed-three-source-quark-pair.py"))
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    interpolation = run_path(str(root/"ch39-hill-quark-common-polynomial.py"))
    product = interpolation["product"]
    coords = interpolation["four_coordinates"]
    ctx.prec = 448
    height = acb.zeta_zero(2_306_888_100_085).imag
    sites = [ray["soft"]["pole"],*ray["boundary_sites"]()]
    assert sites == [14,8,7,6]
    six = ray["laplacian"](6,height).eig(algorithm="rump")
    assert all(z.imag.contains(0) for z in six)
    nodes = sorted(((z.real-18)/8 for z in six),key=lambda x:x.mid())
    assert all(a < b for a,b in zip(nodes,nodes[1:]))
    slopes = [product(a-b for j,b in enumerate(nodes) if i!=j) for i,a in enumerate(nodes)]
    frame = common["adapted_frame"](native["gram_cholesky"](root))
    change = interpolation["compression"](frame,common)
    rank,card = [],[]
    for charge in (-4,2):
        levels, vectors = ray["laplacian"](charge,height).eig(right=True,algorithm="rump")
        assert all(x.imag.contains(0) for x in levels)
        kernel = arb_mat(6,12)
        labels = arb_mat(6,19)
        for i,z in enumerate(levels):
            t = (z.real-18)/8
            chi = product(t-v for v in nodes)
            assert not chi.contains(0)
            vec = acb_mat([[vectors[j,i]] for j in sites])
            norm2 = sum(abs(vectors[j,i])**2 for j in range(19))
            p = arb_mat([[x] for x in coords(vec*common["adj"](vec)*(1/norm2))])
            kernel += p*arb_mat([[chi*t**k for k in range(12)]])
            labels += p*arb_mat([[chi/((t-v)*den) for v,den in zip(nodes,slopes,strict=True)]])
        rank.extend(kernel.tolist())
        card.extend(labels.tolist())
    h = arb("5.66600e-35").union(arb("5.66604e-35"))
    k = arb(".635595244").union(arb(".635595246"))
    s = arb("3.291030650e-40").union(arb("3.291030654e-40"))
    powers = [((8*t+18)/(8*nodes[-1]+18))**6000 for t in nodes[:-1]]+[arb(1)]
    x = [h+k*v/(v+s) for v in powers]
    return arb_mat(rank), arb_mat(card), change, x, slopes



def midpoint(a):
    return arb_mat([[v.mid() for v in row] for row in a.tolist()])


def columns(a, ids):
    return arb_mat([[row[i] for i in ids] for row in a.tolist()])


def frobenius_map(change, heavy):
    rows = []
    for i,y in enumerate(heavy):
        for j,row in enumerate(change.tolist()):
            side = [z*(1 if j<3 else arb(2).sqrt())/y for z in row]
            rows.append(side+[0]*6 if i==0 else [0]*6+side)
    return arb_mat(rows)


def certify(C,L,B,values,count=8):
    # 中点最小二乘残差只用于选择分离方向；在四站坐标中精确消去核的前 count 列。
    midB = midpoint(B)
    c = columns(C,range(count))
    design = midB*midpoint(c)
    y = midB*midpoint(L)*arb_mat([[z.mid()] for z in values])
    normal = design.transpose()*design
    resid = y-design*normal.solve(design.transpose()*y)
    trial = midpoint(resid.transpose()*midB)
    # 在四站行空间中，投影到左零空间。C 的误差来自原磁联络，不含 Gram 转换的误差。
    from itertools import combinations
    pivot = max(combinations(range(12),count),
                key=lambda cols:abs(columns(midpoint(c).transpose(),cols).det()).mid())
    rest = [i for i in range(12) if i not in pivot]
    free = arb_mat([[trial[0,i]] for i in rest])
    z = -columns(c.transpose(),pivot).solve(columns(c.transpose(),rest)*free,algorithm="precond")
    parts = {**dict(zip(pivot,z.entries())),**dict(zip(rest,free.entries()))}
    dual = arb_mat([[parts[i] for i in range(12)]])
    # 一次合并 Lagrange 基系数，避免在同一截点误差的重复副本之间做差。
    v = (dual*L*arb_mat([[z] for z in values]))[0,0]
    norm2 = sum(a*a for a in B.transpose().solve(dual.transpose(),algorithm="precond").entries())
    return v*v/norm2


def phase_system(A,B,t,vnorm,dual,intercepts):
    def adj(m):
        return m.conjugate().transpose()
    pencil = A*(1/vnorm[0])+B*(t/vnorm[1])
    ev, evec = pencil.eig(right=True,algorithm="rump")
    assert all(v.imag.contains(0) for v in ev)
    order = sorted(range(3),key=lambda i:ev[i].real.mid())
    assert all(ev[i].real < ev[j].real for i,j in zip(order,order[1:]))
    y = acb_mat([[evec[i,order[0]]] for i in range(3)])
    def unit(u):
        r=sum(abs(z)**2 for z in u.entries()).sqrt()
        assert r>0
        return u*(1/r)
    y=unit(y)
    ay=A*y
    z=unit((ay-y*(adj(y)*ay)[0,0])*acb(0,1))
    x=[(y[(k+1)%3,0]*z[(k+2)%3,0]-y[(k+2)%3,0]*z[(k+1)%3,0]).conjugate()
       for k in range(3)]
    vals=[]
    for phase in (acb(1),acb(-1),acb(0,1)):
        frame=acb_mat([[phase*x[i],y[i,0],z[i,0]] for i in range(3)])
        assert (adj(frame)*frame).contains(acb_mat([[int(i==j) for j in range(3)] for i in range(3)]))
        data=[]
        for src in (A,B):
            M=adj(frame)*src*frame
            assert M[0,1].contains(0) and M[1,2].real.contains(0)
            data += [M[i,i].real for i in range(3)]+[M[0,2].real,M[0,2].imag,M[1,2].imag]
        vals.append(dual*arb_mat([[v] for v in data])-intercepts)
    offset=(vals[0]+vals[1])*(arb(1)/2)
    D=arb_mat(list(zip((vals[0]-offset).entries(),(vals[2]-offset).entries(),strict=True)))
    return [ev[i].real for i in order],D,offset


def physical_brackets(root,C,L,T,x,p):
    ckm=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    ctx.prec=448
    V=acb_mat(ckm["standard_matrix"](*ckm["centers"]))
    assert all(abs(z)>0 for z in V.entries())
    yup = [arb(z)**p for z in (".00000615",".00311",".8616")]
    yd = [arb(z)**p for z in (".0000135",".000268",".01401")]
    diag=lambda a:acb_mat([[a[i] if i==j else 0 for j in range(3)] for i in range(3)])
    A = V.conjugate().transpose()*diag(yup)*V
    B = diag(yd)
    dual=C.transpose().solve(arb_mat([[int(i==j) for j in (10,11)] for i in range(12)]),
                             algorithm="precond").transpose()
    alpha=(dual[0,0]/dual[1,0]).mid()
    dual=arb_mat([[(dual[0,k]-alpha*dual[1,k]) for k in range(12)],dual.tolist()[1]])
    shift=dual*L*arb_mat([[z.sqrt() if p==1 else z] for z in x])
    traces=[]
    for i in range(2):
        trace=[]
        for j in (0,1):
            v=arb_mat([dual.tolist()[i][6*j:6*j+6]])*T.inv()
            trace.extend((v*(1/shift[i,0].mid())).entries())
        traces.append(trace)
    readings=arb_mat(traces)
    shift=arb_mat([[s*(1/s.mid())] for s in shift.entries()])
    ends = (".19",".20") if p==1 else ("-1","-.99")
    levels=[]
    for index,t in enumerate(map(arb,ends)):
        ev,D,offset=phase_system(A,B,t,(yup[-1],yd[-1]),readings,shift)
        levels.append(ev)
        if p==1:
            arg=-D.solve(offset,algorithm="precond")
            value=sum(z*z for z in arg.entries())-1
            assert arb("4.52") < value < arb("4.54") if index==0 else arb("-.201") < value < arb("-.199")
            print("振幅，铅笔参数与椭圆条件",t.str(5),value.str(9))
        else:
            radius=(D[1,0]**2+D[1,1]**2).sqrt()
            bound=offset[1,0]+arb(0,radius.upper())
            assert bound < arb("-.0157") if index==0 else bound > arb(".1730")
            print("平方质量，铅笔参数与最高次系数的归一范围",t.str(5),bound.str(9))
    # 铅笔随 t 在厄米序下递增，右端最低值低于左端次低值，故全段最低线简单。
    assert levels[1][0] < levels[0][1]


def verify(root):
    from itertools import product
    C,L,T,x,slopes=moment_system(root)
    roots=[v.sqrt() for v in x]
    for p,cut in ((1,roots),(2,x)):
        leading=sum(c/w for c,w in zip(cut,slopes,strict=True))
        assert leading < -arb("4.34e8") if p==1 else leading < -arb("2.11e8")
        print("荷六单独插值的十八次系数",p,leading.str(11))
        heavy=[arb(".8616")**p,arb(".01401")**p]
        B=frobenius_map(T,heavy)
        if p==1:
            # 整体变号保留平方范数，故第十九层固定正。前十六层以球包含任意符号。
            for a,b in product((1,-1),repeat=2):
                values=[arb(0,v.upper()) for v in cut[:16]]+[a*cut[16],b*cut[17],cut[18]]
                lower=certify(C,L,B,values)
                assert lower>113
                print("第十七、十八层符号及 26 次加权迹平方和下界",a,b,lower.lower().str(12))
        else:
            lower=certify(C,L,B,cut)
            assert lower>42
            print("平方质量的 26 次加权迹平方和下界",lower.lower().str(12))
        windows=((".00000629",".00316",".8659"),(".0000137",".000271",".01410"))
        maximum=sum((arb(y)**p/h)**2 for rows,h in zip(windows,heavy,strict=True) for y in rows)
        assert maximum < (arb("2.024") if p==1 else arb("2.047"))
        physical_brackets(root,C,L,T,x,p)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
