"""由两端荷六完整谱投影检验共同稳定子与字符内按能阶匹配的极分解。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def groups(vectors, rotation):
    omega = acb(-arb(1)/2,arb(3).sqrt()/2)
    result = [[],[],[]]
    for j in range(19):
        phase = sum(vectors[i,j].conjugate()*vectors[rotation[i],j] for i in range(19))
        # 拉普拉斯对三循环不变，简单本征线具有且只具有一个字符。
        k, = [k for k,z in enumerate((1,omega,omega**2)) if abs(phase-z)<arb(".01")]
        result[k].append(j)
    assert [len(row) for row in result] == [7,6,6]
    return result


def verify_pair(spectrum, shell, rotation, source, target, gap_floor):
    ctx.prec = 384
    left = spectrum(shell,6,source)[1]
    right = spectrum(shell,6,target)[1]
    incoming, outgoing = [groups(v,rotation) for v in (left,right)]
    overlap = adj(right)*left
    eigen = []
    for old,new in zip(incoming,outgoing,strict=True):
        least = min(abs(overlap[i,j]).lower() for i in new for j in old)
        assert least > arb(".0003")
        matched = [overlap[i,j] for i,j in zip(new,old,strict=True)]
        assert all(abs(z)>arb(".78") for z in matched)
        # G 的旧壳字符基矩阵。每一列为对应新谱矢，乘以消去选相的重叠相位。
        g = acb_mat([[overlap[i,j].conjugate()*z/abs(z) for i,z in zip(new,matched,strict=True)]
                     for j in old])
        identity = acb_mat([[int(k==j) for j in old] for k in old])
        assert (adj(g)*g).contains(identity)
        eigen.append(g.eig(algorithm="rump"))
        print("字符维数、所有交叠模与配对交叠模的下界",len(old),least.str(12),
              min(abs(z).lower() for z in matched).str(12),flush=True)
    gap = min(abs(z-w).lower() for z in eigen[0] for w in eigen[1]+eigen[2])
    assert gap > arb(gap_floor)
    print("零点编号、平凡块与两字符块的最小本征弦距",source,target,gap.str(15),flush=True)
    return incoming,outgoing


def gravity(root, old, new):
    # 右内壳的三个位置沿字符内的升序配对送达右外壳。
    ordered = run_path(str(root/"ch39-adjacent-ordered-charge-six-cutoffs.py"))["verify"]
    _, positions = ordered(root)[2]
    assigned = [None]*19
    for source,target in zip(old,new,strict=True):
        for j,k in zip(source,target,strict=True):
            assigned[k] = positions[j]
    assert all(x is positions[0] for x in assigned[:17])
    assert assigned[17] is positions[17] and assigned[18] is positions[18]
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    high_weak = run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"]
    n = 2_306_888_100_087
    _,_,raw,_ = next(row for row in native["shell_weights"](root,n) if row[0]==6)
    c = [z.real for z in raw[0]]
    co = (sum(c[:17]),c[17],c[18])
    moments,heat,fields = run_path(str(root/"ch39-three-index-free-atoms.py"))["broken_fields"](root,[n])
    p0 = (arb("91.1876")/arb("1.22089e19"))**2
    p1 = (arb("1.176397e14")/arb("1.22089e19"))**2
    low,middle,high = [assigned[i] for i in (0,17,18)]
    thermal = [moments(fields[0],heat,x) for x in (low,middle,high)]
    weak = (native["weak_series"]()(low/p0),high_weak(middle,p0,p1),high_weak(high,p0,p1))
    k,w,d = [sum(a*b for a,b in zip(co,col,strict=True)) for col in (*zip(*thermal),weak)]
    assert arb("1.98859") < k < arb("1.98860")
    assert arb("-.85977") < w < arb("-.85975")
    assert arb(".03960") < d < arb(".03961")
    print("右外壳牛顿、真空、弱差读数",k.str(12),w.str(12),d.str(12),flush=True)


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    spectrum = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    sites = shell["geo"]["sites"]
    cycle = [sites.index((y,z,x)) for x,y,z in sites]
    n = 2_306_888_100_085
    old,new = verify_pair(spectrum,shell,cycle,n+1,n+2,".00151554")
    verify_pair(spectrum,shell,cycle,n-1,n-2,".00311815")
    gravity(root,old,new)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
