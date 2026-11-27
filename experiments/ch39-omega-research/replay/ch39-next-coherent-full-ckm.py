"""N+3 全相干独立酉输运：原目标的分离证书与共同旋转后的六质量、完整 CKM。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


QUATERNION = (-39764766276, -82568934973, -530684641478, 842600110159)
ROTATED_SEEDS = (
    (None,947553448968,1920027505,None,33100075377,947553448968,9185572984,
     1920027505,947553448968,996777656,None,8189807417,1920027505,947553448968,
     146907274974,None,None,None,None),
    (None,313681557,10512242,153774721094,55878013162,153774721094,None,
     153774721094,148334523689,4858469521,153774721094,None,223962931475,153774721094,
     42509077992,None,None,None,None),
)
ORIGINAL_DOWN = (
    None,219426817081,11453242159,153774721094,113277179137,369701149,None,
    153774721094,153774721094,10512242,153774721094,None,153771401142,
    153774721094,32681931295,None,None,None,None,
)


def adj(matrix):
    return matrix.conjugate().transpose()


def block_from_coordinates(coords):
    tr,xx,yy,zz,plus,minus = coords
    off = acb(xx,-yy)/2
    char = acb(0,minus-plus)/2
    return acb_mat([[(tr+zz)/2,off,0,0],[off.conjugate(),(tr-zz)/2,0,0],
                    [0,0,(plus+minus)/2,char],[0,0,char.conjugate(),(plus+minus)/2]])


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))

    embedding = load("common-embedding-central-right-ckm")
    independent = load("independent-unitary-quark-right-holdout")
    n = independent["CENTER"]
    _,source,_ = embedding["verify_nodes"](
        root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    goals = load("one-generator-gravity-full-ckm")["targets"](root,independent,embedding)
    projectors = load("free-four-site-native-positive-spectrum")["native"](root,n+3)
    ctx.prec = 384
    blocks = [block_from_coordinates(goal.entries()) for goal in goals]
    _,_,intrinsic,old_rotation = embedding["common_embedding"](*independent["GEOMETRY"])
    old_frame = old_rotation*intrinsic
    coords = independent["coordinates"]
    separator = (-8,-3,-8,-6)

    def linear(row):
        return sum(a*b for a,b in zip(separator,coords(row)[:4],strict=True))

    trivial_lines = []
    # 平凡签名为零的非平凡字符，由循环对称性使整个平凡 2x2 块恰为零。
    for i,row in enumerate(projectors[0]):
        if (row[0,0]+row[1,1]).real > 0:
            trivial_lines.append(i+1)
            assert linear(row)<arb("-.001")
        else:
            assert all(row[a,b].contains(0) for a in (0,1) for b in (0,1))
    assert len(trivial_lines)==7 and linear(blocks[0])>arb(".00051125")
    reconstruct = independent["reconstruct"]
    _,margin = reconstruct(projectors[1],blocks[1],source[n+1,2],ORIGINAL_DOWN)
    print("原上型分离式及原下型主序余量",linear(blocks[0]).str(13),margin.str(13),flush=True)

    norm = arb(sum(j*j for j in QUATERNION)).sqrt()
    x,y,z,w = [arb(j)/norm for j in QUATERNION]
    rotation = acb_mat([[acb(w,-z),acb(-y,-x),0,0],[acb(y,-x),acb(w,z),0,0],
                        [0,0,1,0],[0,0,0,1]])
    eye = load("rank-two-generator-up-quark")["eye"]
    assert (adj(rotation)*rotation).contains(eye(4))
    frame = rotation*old_frame
    assert (adj(frame)*frame).contains(eye(3))
    flavor = load("three-gram-source-flavor")["projections"]
    ps = []
    masses = [list(map(arb,s)) for s in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    for charge,rows,old,seed,reference in zip((-4,2),projectors,blocks,ROTATED_SEEDS,masses,strict=True):
        desired = rotation*old*adj(rotation)
        raw,margin = reconstruct(rows,desired,source[n+1,charge],seed)
        physical = adj(frame)*raw*frame
        delta = physical-adj(old_frame)*old*old_frame
        assert delta.contains(acb_mat(3,3)) and max(abs(v).upper() for v in delta.entries())<arb("1e-38")
        values,atoms = flavor(physical*physical)
        assert all(found.overlaps(expected**2) for found,expected in zip(values,reference,strict=True))
        ps.append(atoms)
        print("荷、主序余量、共同三平面质量",charge,margin.str(13),
              *[v.sqrt().str(13) for v in values],flush=True)
    up,down = ps
    mixing = acb_mat(load("shell-clock-ckm-blocks")["standard_matrix"](
        *[arb(a).sin() for a in (".22704",".04275",".00377")],arb("1.139")))
    assert all((u*d).trace().real.overlaps(abs(mixing[i,j])**2)
               for i,u in enumerate(up) for j,d in enumerate(down))
    signed = (up[0]*down[0]*up[1]*down[1]).trace().imag
    standard = (mixing[0,0]*mixing[1,1]*mixing[0,1].conjugate()*mixing[1,0].conjugate()).imag
    assert signed>0 and signed.overlaps(standard)
    print("两荷九个模平方及带号 J",signed.str(15),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
