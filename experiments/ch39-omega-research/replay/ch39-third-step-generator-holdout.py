"""从两壳联合见证的整参数盒检验 N+4 两条预先指定的全相干续接。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def heavy_centers(blocks,levels,eye):
    root3 = arb(3).sqrt()
    omega = acb(-1,root3)/2
    trivial = acb_mat([[1],[acb(0,-arb("3.1")/root3)],[acb(0,-arb("3.1")/root3)],
                       [acb(0,-arb("3.1")/root3)]])*(1/arb("10.61").sqrt())
    chars = [acb_mat([[0],[1],[o],[o**2]])*(1/root3) for o in (omega,omega.conjugate())]
    directions = [trivial,*chars]
    rows = [[arb(1)]*3]
    for raw,lev,mass in zip(blocks,levels,map(arb,(".8616",".01401")),strict=True):
        assert lev[2]<mass<lev[3]
        below = mass*eye(4)-raw
        cofactor = below.inv()*below.det()
        assert all((adj(v)*cofactor*w)[0,0].contains(0) for j,v in enumerate(directions) for w in directions[j+1:])
        rows.append([(adj(v)*cofactor*v)[0,0].real for v in directions])
    system = arb_mat(rows)
    assert not system.det().contains(0)
    p = system.solve(arb_mat([[1],[0],[0]]),algorithm="precond").entries()
    assert min(p)>arb(".0789") and sum(p).contains(1)
    print("顶底两中心公共平面的三个法向概率",[v.str(12) for v in p],flush=True)


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    joint = lib("next-common-ckm-gravity")
    cx,rx,cy,ry,_ = joint["verify"](root)
    rank = lib("rank-two-generator-up-quark")
    emb = lib("common-embedding-central-right-ckm")
    n = rank["N"]
    _,trained,_ = emb["verify_nodes"](root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    eye = rank["eye"](19)
    local = lib("next-common-ckm")["local"]
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"]
    upstream = lib("one-generator-gravity-full-ckm")
    def whole(name,center,radius,directions):
        wide = [z+arb(0,radius.upper()) for z in center]
        return local(seed(root,rank["eye"],name),wide,directions,eye,False)[0]
    x = whole("ch39-full-ckm-gravity-cayley.tsv",cx,rx,upstream["DIRECTIONS"])
    y = whole("ch39-next-common-ckm-gravity-cayley.tsv",cy,ry,joint["DIRECTIONS"])
    reference = lib("next-coherent-full-ckm")
    quat = reference["QUATERNION"]
    length = arb(sum(i*i for i in quat)).sqrt()
    a,b,c,d = [arb(i)/length for i in quat]
    rot = acb_mat([[acb(d,-c),acb(-b,-a),0,0],[acb(b,-a),acb(d,c),0,0],[0,0,1,0],[0,0,0,1]])
    s,t,u = [1/arb(j).sqrt() for j in (3,2,6)]
    site_basis = acb_mat([[1,0,0,0],[0,s,t,u],[0,s,-t,u],[0,s,0,-2*u]])
    independent = lib("independent-unitary-quark-right-holdout")
    _,_,intrinsic,orig = emb["common_embedding"](*independent["GEOMETRY"])
    frame = site_basis*rot*orig*intrinsic
    assert (adj(frame)*frame).contains(rank["eye"](3))
    windows = lib("three-shell-outer-holdout")["WINDOWS"]
    sorted_eigen = lib("full-ckm-generator-next-shell")["sorted_eigen"]
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    sites = (shell["soft"]["pole"],*shell["boundary_sites"]())
    assert sites==(14,8,7,6)
    spectra = []
    for q in (-4,2):
        src,dst = [spec(shell,q,v)[1] for v in (n+1,n+4)]
        four = acb_mat([[dst[i,j] for j in range(19)] for i in sites])
        spectra.append((src,dst,four,trained[n+1,q]))
    angles = [arb(m)+arb(0,arb(e).upper()) for m,e in
              ((".22704",".00082"),(".04275",".00042"),(".00377",".00008"),("1.139",".023"))]
    ref = acb_mat(lib("shell-clock-ckm-blocks")["standard_matrix"](*[v.sin() for v in angles[:3]],angles[3]))
    jband = arb("3.23e-5")+arb(0,arb(".08e-5").upper())
    for rule,z in (("延用左增量",y*adj(x)*y),("固定 Y",y)):
        vecs,blocks,four_spectra = [],[],[]
        for iq,(op,(src,dst,four,w)) in enumerate(zip((adj(z*z)*adj(y*y)*adj(x*x),z*y*x),spectra,strict=True)):
            response = rank["read"](src,dst,w,op)
            raw = acb_mat([[four[i,j]*response[j] for j in range(19)] for i in range(4)])*adj(four)
            raw_levels,_ = sorted_eigen(raw)
            levels,vec = sorted_eigen(adj(frame)*raw*frame)
            excluded = [raw_levels[j]>arb(hi) or raw_levels[j+1]<arb(lo) for j,(lo,hi) in enumerate(windows[iq])]
            outside = [v>arb(hi) or v<arb(lo) for v,(lo,hi) in zip(levels,windows[iq],strict=True)]
            print(rule,"荷",(-4,2)[iq],"四站谱",[v.str(11) for v in raw_levels],"任意三平面排除",excluded,flush=True)
            print(rule,"三代质量",[v.str(11) for v in levels],"原平面出窗",outside,flush=True)
            assert all(outside) and excluded==[True,True,rule=="延用左增量"]
            blocks.append(raw)
            four_spectra.append(raw_levels)
            vecs.append(vec)
        pu,pd = vecs
        mixing = acb_mat([[(adj(a)*b)[0,0] for b in pd] for a in pu])
        prob = [[abs(mixing[i,j])**2 for j in range(3)] for i in range(3)]
        signed = (mixing[0,0]*mixing[1,1]*mixing[0,1].conjugate()*mixing[1,0].conjugate()).imag
        outside = [[not prob[i][j].overlaps(abs(ref[i,j])**2) for j in range(3)] for i in range(3)]
        print(rule,"九个模平方",[[v.str(10) for v in row] for row in prob],
              "逐项出窗",outside,"有符号 J",signed.str(12),flush=True)
        assert all(all(row) for row in outside) and not signed.overlaps(jband)
        if rule=="固定 Y":
            heavy_centers(blocks,four_spectra,rank["eye"])


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
