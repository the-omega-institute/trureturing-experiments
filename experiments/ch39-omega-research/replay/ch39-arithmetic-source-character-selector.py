"""四个原始 Gram 源与三个顺序枢轴的字符支撑，以及以源固定两壳谱矢相位的生成元。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def anchored_frame(eigenvectors, state):
    overlaps = adj(eigenvectors)*state
    minimum = min(abs(v).lower() for v in overlaps.entries())
    assert minimum > arb(".0004")
    phases = [v/abs(v) for v in overlaps.entries()]
    frame = acb_mat([[eigenvectors[i,j]*phases[j] for j in range(19)] for i in range(19)])
    return frame,minimum


def verify(root):
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    spectrum = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    algebra = run_path(str(root/"ch39-charge-six-projector-algebra-transport.py"))
    common = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))
    ordered = run_path(str(root/"ch39-adjacent-ordered-charge-six-cutoffs.py"))["verify"]
    weighted = run_path(str(root/"ch39-common-unitary-quark-right-trace.py"))["weighted"]
    up = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    n = 2_306_888_100_086
    _,labels,_ = common["verify_nodes"](root,((n,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    _,positions = ordered(root)[2]
    source_levels = [z.sqrt() for z in positions]
    moments,heat,fields = run_path(str(root/"ch39-three-index-free-atoms.py"))["broken_fields"](root,[n+1])
    _,_,band,_ = next(row for row in native["shell_weights"](root,n+1) if row[0]==6)
    read = [z.real for z in band[0]]
    weak_high = run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"]
    ctx.prec = 384
    weak_low = native["weak_series"]()
    p0=(arb("91.1876")/arb("1.22089e19"))**2
    p1=(arb("1.176397e14")/arb("1.22089e19"))**2
    native_sites = shell["geo"]["sites"]
    rotate = [native_sites.index((y,z,x)) for x,y,z in native_sites]
    algebra["verify_pair"](spectrum,shell,rotate,n,n+1,".00151554")
    gram = native["gram_cholesky"](root)
    sites = (shell["soft"]["pole"],*shell["boundary_sites"]())
    assert sites == (14,8,7,6)
    left,right = [spectrum(shell,6,k)[1] for k in (n,n+1)]
    src_up,dst_up = [spectrum(shell,-4,k)[1] for k in (n,n+1)]
    four = acb_mat([[dst_up[i,j] for j in range(19)] for i in sites])
    eye = up["eye"](19)
    groups = algebra["groups"]
    old,new = [groups(v,rotate) for v in (left,right)]
    assert set(old[0])-set(new[0]) == {5}
    assert set(new[0])-set(old[0]) == {3}
    assert 3 in old[2] and 5 in new[1]
    trivial = acb_mat([[left[i,j] for j in old[0]] for i in range(19)])
    projection = trivial*adj(trivial)
    residuals = [[int(k==i) for k in range(4)] for i in (1,2,3)]
    for index,coeff in enumerate(gram+residuals):
        vector = acb_mat([[coeff[sites.index(i)] if i in sites else 0] for i in range(19)])
        vector *= 1/(adj(vector)*vector)[0,0].real.sqrt()
        fixed = abs(vector[14,0])**2+abs(sum(vector[i,0] for i in sites[1:]))**2/3
        nontrivial = (1-fixed)/2
        if index == 0:
            assert all(x.is_zero() for x in gram[0][1:]) and fixed.contains(1)
            print("g0 只含平凡字符",flush=True)
            continue
        assert min(fixed,nontrivial)>arb(".00053")
        if index >= 4:
            assert (3*fixed).contains(1) and (3*nontrivial).contains(1)
        label = ("g"+str(index)) if index<4 else ("f"+str(index-3))
        old_frame,old_min = anchored_frame(left,vector)
        new_frame,new_min = anchored_frame(right,vector)
        x = new_frame*adj(old_frame)
        assert (adj(x)*x).contains(eye)
        u = adj(x*x)
        leakage = (eye-projection)*u*trivial
        # X 与其逆各只送一条平凡谱线到补空间，故逆平方的跨块秩至多二。
        # 由两个迹复原至多两个非零本征值，可避开另外五个零根。
        a = adj(leakage)*leakage
        tr = a.trace().real
        discriminant = 2*(a*a).trace().real-tr**2
        assert discriminant > 0
        p = (tr+discriminant.sqrt())/2
        assert p.contains(1) and p > arb(".999999999")
        up_amplitudes = up["read"](src_up,dst_up,labels[n,-4],u)
        up_levels = sorted((v.real for v in weighted(four,up_amplitudes).eig(algorithm="rump")),key=lambda z:z.mid())
        assert up_levels[0] > arb(".00000629")
        g=x*x*x
        grav = up["read"](left,right,source_levels,g)
        assert min(grav)**2 > p0
        thermal = [moments(fields[0],heat,a**2) for a in grav]
        weak = [weak_low(a**2/p0) if a**2 < 4096*p0 else weak_high(a**2,p0,p1) for a in grav]
        k,w,d = [sum(a*b for a,b in zip(read,col,strict=True)) for col in (*zip(*thermal),weak)]
        assert grav[15]>arb(".05") and d<arb("-.01656")
        print(label,"字符质量",fixed.str(13),nontrivial.str(13),
              "两端最小谱投影范数",old_min.str(8),new_min.str(8),
              "上型转移、四站最小值",p.str(12),up_levels[0].str(12),
              "荷六第16层、K、W、D",*[v.str(11) for v in (grav[15],k,w,d)],flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
