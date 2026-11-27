"""算术射线强度与谱能序配对的三次幂，验证右外壳第十六层和弱差。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def strength_order(frame, state):
    powers = [abs(z)**2 for z in (adj(frame)*state).entries()]
    order = sorted(range(19), key=lambda i: powers[i].mid())
    assert all(powers[i] < powers[j] for i,j in zip(order,order[1:]))
    return order


def pairing(old, new):
    perm = [None]*19
    for i,j in zip(old,new,strict=True):
        perm[i] = j
    return perm


def supports(pair, colors, dest):
    src,dst = colors
    q = pair.index(dest)
    return [[i for i in range(19) if src[i]==dst[pair[h]] and dst[pair[i]]==src[q]]
            for h in (17,18)]


def polygon(v, perm, colors, dest=15):
    q = perm.index(dest)
    rows = supports(perm,colors,dest)
    coeff = [[v[q,perm[k]]*v[k,perm[h]] for k in ks] for h,ks in zip((17,18),rows,strict=True)]
    bounds = []
    for row in coeff:
        lengths = list(map(abs,row))
        total = sum(lengths,arb(0))
        excess = 2*max(lengths,default=arb(0))-total
        assert not lengths or not excess.contains(0)
        bounds.append((max(arb(0),excess.lower()),total))
    return rows, coeff, bounds


def pair_to_real(z1, z2, radius):
    a,b = abs(z1),abs(z2)
    x = (radius**2+a**2-b**2)/(2*radius)
    y = (a**2-x**2).sqrt()
    assert y>0
    return acb(x,y)/z1, acb(radius-x,-y)/z2


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    shell = lib("three-ray-chiral-boundary")
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    native = lib("native-rh-ray-pvm")
    algebra = lib("charge-six-projector-algebra-transport")
    anchor = lib("arithmetic-source-character-selector")["anchored_frame"]
    n = 2_306_888_100_086
    levels = [z.sqrt() for z in lib("adjacent-ordered-charge-six-cutoffs")["verify"](root)[2][1]]
    moments,heat,fields = lib("three-index-free-atoms")["broken_fields"](root,[n+1])
    _,_,bands,_ = next(row for row in native["shell_weights"](root,n+1) if row[0]==6)
    weights = [z.real for z in bands[0]]
    weak_high = lib("full-neutral-free-cutoff")["weak_column"]
    ctx.prec = 384
    weak_low = native["weak_series"]()
    mp,mz = arb("1.22089e19"),arb("91.1876")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    left,right = [spectrum(shell,6,k)[1] for k in (n,n+1)]
    groups = [algebra["groups"](v,rotation) for v in (left,right)]
    colors = [[next(k for k,group in enumerate(groups_at_end) if i in group)
               for i in range(19)] for groups_at_end in groups]
    raw = adj(left)*right
    ranks = list(range(19))
    # 同字符所有交叠模严格为正，才可由字符路径判定精确的支撑。
    assert min(abs(raw[i,j]).lower() for i in ranks for j in ranks
               if colors[0][i]==colors[1][j]) > arb(".0003")
    ranked = [sorted(ranks,key=lambda j: (abs(raw[i,j])**2).mid()) for i in ranks]
    optimal = [row[-1] for row in ranked]
    assert len(set(optimal))==19
    assert all(abs(raw[i,row[-1]])**2-abs(raw[i,j])**2 > arb(".2642")
               for i,row in enumerate(ranked) for j in row[:-1])
    assert all(colors[0][i]==colors[1][j] for i,j in enumerate(optimal))
    _,_,optimal_bounds = polygon(raw,optimal,colors)
    assert optimal_bounds[0][0]==0 and optimal_bounds[0][1]>arb(".447")
    assert optimal_bounds[1][1].is_zero()
    print("最优重叠配对、两源路径模长界", [j+1 for j in optimal], optimal_bounds,flush=True)

    gram = native["gram_cholesky"](root)
    sites = (14,8,7,6)
    for index,coeff in enumerate(gram[1:]+[[0,1,0,0]],1):
        state = acb_mat([[coeff[sites.index(i)] if i in sites else 0] for i in ranks])
        state *= 1/(adj(state)*state)[0,0].real.sqrt()
        s,t = [anchor(v,state)[0] for v in (left,right)]
        v = adj(s)*t
        ins,outs = [strength_order(frame,state) for frame in (s,t)]
        same = pairing(ins,outs)
        _,_,bounds = polygon(v,same,colors)
        lower = sum(lo**2 for lo,hi in bounds)
        assert lower>arb(".70")
        backward = pairing(ins,outs[::-1])
        _,_,bounds = polygon(v,backward,colors)
        reverse_lower = sum(lo**2 for lo,hi in bounds)
        assert reverse_lower>arb("1.3e-9")
        perm = pairing(ins,ranks)
        rows,coeffs,bounds = polygon(v,perm,colors)
        _,_,bounds18 = polygon(v,perm,colors,17)
        ceiling18 = levels[0]+sum((h-levels[0])*hi**2 for h,(_,hi) in zip(levels[17:],bounds18,strict=True))
        assert ceiling18<arb(".0121")
        available = sum(hi**2 for lo,hi in bounds)
        base = sum(abs(sum(row))**2 for row in coeffs)
        print("源编号、同强度与反序概率下界、源强度到能序的固定相位概率",
              index,lower.str(10),reverse_lower.str(10),arb(base).str(12),
              "第18层振幅上界",ceiling18.upper().str(10),flush=True)
        if index in (2,3):
            assert [j+1 for j in ranks if not any(supports(perm,colors,j))]==[2,3,5,8,9,16]
        if index==4:
            assert sum(lo**2 for lo,hi in bounds)>arb(".0475")
        if index!=1:
            continue
        assert rows==[[4,8,11,14],[]] and perm[17:]==[18,13] and perm.index(15)==12
        assert bounds[0][0]==0 and available>arb("4.75e-6") and available<arb("4.77e-6")
        middle,highest = levels[17:]
        low = levels[0]
        four = coeffs[0]
        a,b,c,d = map(abs,four)
        radius = (max(abs(a-b),abs(c-d))+min(a+b,c+d))/2
        first = pair_to_real(four[0],four[1],radius)
        second = pair_to_real(four[2],four[3],radius)
        print("四条路径长度、半径、最大概率",*[z.str(12) for z in (a,b,c,d,radius,available)],flush=True)

        def read(q, report=True):
            probability = (q/mp-low)/(middle-low)
            assert probability>0
            # 第一对合成 r，第二对合成 r u，使 |r(1+u)|^2=b。
            cosine = probability/(2*radius**2)-1
            sine = (probability/radius**2*(1-probability/(4*radius**2))).sqrt()
            turn = acb(cosine,sine)
            phases = [acb(1) for _ in ranks]
            for i,z in zip(rows[0],(*first,*(turn*z for z in second)),strict=True):
                phases[i] = z
            # 在算术源固定的坐标中，三次幂两次经过 V。首尾的单位相位不影响模方。
            values = []
            for j in ranks:
                last = perm.index(j)
                paths = supports(perm,colors,j)
                atoms = [sum((v[last,perm[k]]*phases[k]*v[k,perm[h]] for k in ks),acb(0))
                         for h,ks in zip((17,18),paths,strict=True)]
                if j==15:
                    if report:
                        assert (abs(atoms[0])**2).overlaps(probability) and atoms[1].is_zero()
                    values.append(q/mp)
                else:
                    values.append(low+(middle-low)*abs(atoms[0])**2+(highest-low)*abs(atoms[1])**2)
            assert min(x for j,x in enumerate(values) if j!=15)>arb("1e-12")
            cols = [moments(fields[0],heat,x**2) for x in values]
            weak = [weak_low((q/mz)**2) if j==15 else weak_high(x**2,p0,p1)
                    for j,x in enumerate(values)]
            k,w,delta = [sum(a*b for a,b in zip(weights,col,strict=True)) for col in (*zip(*cols),weak)]
            if report:
                print("第16层 GeV、概率、K W D",*[z.str(12) for z in (q,probability,k,w,delta)],flush=True)
            return k,w,delta

        left_q,right_q = arb(105),arb(155)
        assert read(left_q)[2]>0 and read(right_q)[2]<0
        assert read(arb("107.850464"))[2]>0 and read(arb("107.850466"))[2]<0
        # 用源低振幅本身给出一个确定的目标概率，无需输入额外的能标。
        q0 = mp*low*(1+middle-low)
        assert arb("104.651")<q0<arb("156.173") and arb("1.7454e-19")<low<arb("1.0131e-17")
        read(q0)
        k,w,_ = read(left_q.union(right_q),False)
        print("整个弱差变号段 K、W 区间",k.str(9),w.str(9),flush=True)
        assert not k.contains(1) and not w.contains(arb("1e-122"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
