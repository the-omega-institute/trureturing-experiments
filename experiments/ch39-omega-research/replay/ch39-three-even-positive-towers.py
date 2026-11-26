"""由三组精确的正电导逆谱解，用一个十四周期余弦公式实现三座偶荷塔。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


NUMBER = 2_306_888_100_085
DATA = {
    2: ((
        6999735515233161, 2170470427213028, 584060727770040, 843772918362397,
        6017482732827891, 196352910099546, 6969147911375079, 6563582593668786,
        142958658079434, 277485469586484, 2745046029600884, 5690230814199892,
        3846681550036405, 306814381910860, 5960958163902379, 6612760562448771,
        1659554907550815, 5625086470618308, 196629535150219, 3240318354052606,
        3750434146778018, 2799516062130604, 3753723788020724, 1962806854285594,
        3576021277920711, 3174613372201782, 4520013953486041, 1694133426982114,
        3761399345946166, 1069053904394441, 1943521593409928, 5332176981962071,
        906091476246077, 4734044252975594, 5543732951562311, 384999320035463,
        3965512255038389, 585867543270518, 607081260420220, 6139853322928003,
        596195774239599, 839991894104594, 673032574180682, 1013561191668629,
        1813199483125014,
    ), (19, 8, 15, 36, 41, 0, 33, 18, 27, 40, 14, 11, 38, 12, 43, 1, 42, 34, 22), (1, 2, 0)),
    4: ((
        7985153020905765, 461018954590175, 2102301385600881, 2303423386060196,
        7930675023166716, 145767388946421, 7915056601112822, 6447788992528374,
        1317939576880151, 125000002243318, 3898961132847858, 6010540723155218,
        1230047269806760, 750369242629327, 5432816005248421, 7990276101078083,
        938602058619528, 2955907014526387, 210337259075150, 3258231090577646,
        4068949220410112, 4085302096826132, 2082899423310404, 3183348262672802,
        4483322490743161, 2049159620720859, 3072827989739741, 2517010777181162,
        2459223361262521, 2111685519485216, 1799516300289676, 4786501545534823,
        378996451471702, 5491530562055153, 3826521465264121, 125004773924647,
        4720879504778644, 400581680699686, 993569615905840, 5816693519679903,
        152927438870491, 1323594394023685, 1183143231741809, 165171961913659,
        1501908519034566,
    ), (19, 16, 43, 36, 17, 21, 40, 15, 30, 2, 34, 38, 37, 13, 18, 12, 9, 4, 11), (2, 0, 1)),
    6: ((
        6499996246379350, 5972762544732165, 618213888495223, 551039568883480,
        6017674465951536, 154517868764868, 6498948672637201, 6292474916637366,
        607602486647025, 220163338831895, 2376558767583506, 6499996521456488,
        153864754301776, 153898511964528, 6499336927445409, 6460849030624627,
        2083773501779342, 6499999999998914, 283296943389151, 347499366782084,
        4186105981095821, 1103715262470633, 4411501426682724, 3005391415574410,
        5422919343636346, 511155343914897, 3292587752727661, 1218681444903943,
        5999790224795687, 646139007970598, 505889201797343, 3396010068664034,
        686529570832008, 6275543458970670, 5570967736685033, 1580660804886402,
        4322627519544732, 979348915747654, 992507492954779, 6044035300568289,
        641871621043011, 1019801535664083, 2641462194199057, 1854820482026086,
        1534385939162316,
    ), (8, 0, 35, 43, 28, 20, 15, 22, 40, 38, 2, 21, 6, 39, 12, 3, 11, 1, 34), (2, 0, 1)),
}


def orbits(graph):
    points = graph["sites"]
    rotate = [points.index((y,z,x)) for x,y,z in points]
    edges = {(i,j) for i,j,_,_ in graph["edges"]}
    groups = []
    while edges:
        edge = min(edges)
        initial = edge
        group = []
        for _ in range(3):
            group.append(edge)
            edges.remove(edge)
            i,j = edge
            edge = tuple(sorted((rotate[i],rotate[j])))
        assert edge==initial
        groups.append(group)
    assert len(groups)==45
    return groups


def certify(native, groups, reference, c, charge):
    numerators,variable,perm = DATA[charge]
    weights = [arb(z)/10**15 for z in numerators]
    radius = arb("2e-11")
    assert len(weights)==45 and len(set(variable))==19
    phases = {}
    matrix = acb_mat(19,19)
    for k,group in enumerate(groups):
        for i,j in group:
            norm = abs(native[i,j])
            assert norm.overlaps(reference[k])
            phase = -native[i,j]/norm
            phases[i,j] = phase
            matrix[i,i] += weights[k]
            matrix[j,j] += weights[k]
            matrix[i,j] = -weights[k]*phase
            matrix[j,i] = matrix[i,j].conjugate()
    original = native.eig(multiple=True)
    assert len(original)==19 and all(z.imag.contains(0) and z.real>0 for z in original)
    folded = sorted((z.real*abs(z.real-c) for z in original),key=float)
    assert all(a<b for a,b in zip(folded,folded[1:])) and folded[2]<9<10<folded[3]
    i,j,k = perm
    target = [(c-(c*c-4*folded[i]).sqrt())/2,(c+(c*c-4*folded[j]).sqrt())/2,
              (c+(c*c+4*folded[k]).sqrt())/2]
    target.extend((c+(c*c+128*z).sqrt())/2 for z in folded[3:])
    assert all(a<b for a,b in zip(target,target[1:]))
    values,vec = matrix.eig(right=True,algorithm="rump")
    order = sorted(range(19),key=lambda k:float(values[k].real))
    eig = [values[k].real for k in order]
    assert all(z.imag.contains(0) and z.real>0 for z in values)
    gap = min((b-a).lower() for a,b in zip(eig,eig[1:]))
    derivative = arb_mat(19,19)
    for row,col in enumerate(order):
        norm = sum(abs(vec[i,col])**2 for i in range(19))
        assert norm>0
        for j,g in enumerate(variable):
            derivative[row,j] = sum(abs(vec[i,col]-phases[i,k]*vec[k,col])**2 for i,k in groups[g])/norm
    pre = derivative.inv()
    norm = max(sum(abs(pre[i,j]) for j in range(19)).upper() for i in range(19))
    correction = max(abs(z).upper() for z in (pre*arb_mat([[a-b] for a,b in zip(eig,target)])).entries())
    change = 6*19*radius
    assert gap>2*change
    contraction = norm*19*12*change/(gap-2*change)
    assert contraction<1 and correction+radius*contraction<radius
    boxes = [w+arb(0,radius) if k in variable else w for k,w in enumerate(weights)]
    assert all(w>0 for w in boxes)
    ratio = max((w/v).max(v/w).upper() for w,v in zip(boxes,reference))
    assert ratio<arb({2:"7",4:"8",6:"6.5"}[charge])
    print("荷、中点最小隙、逆雅可比范数、预乘残差、整盒压缩率、电导反差上界",charge,
          *[z.str(15) for z in (gap,norm,correction,contraction,ratio)])
    print("原三轻质量、新第四与第十九层",*[z.str(13) for z in (*folded[:3],32*folded[3],32*folded[-1])])
    return [(w/v).log() for w,v in zip(boxes,reference)]


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 320
    h = acb.zeta_zero(NUMBER).imag
    groups = orbits(shell["geo"])
    neutral = shell["laplacian"](0,h)
    c = sum(neutral[i,i].real for i in range(19))/18
    reference = [abs(neutral[group[0][0],group[0][1]]) for group in groups]
    assert all(all(abs(neutral[i,j]).overlaps(reference[k]) for i,j in g) for k,g in enumerate(groups))
    logs = [certify(shell["laplacian"](q,h),groups,reference,c,q) for q in (2,4,6)]
    # 偶数荷 q=2k 时使用七阶循环群上的实偶 Fourier 展开。
    freq = [[(2*k*m*arb.pi()/7).cos() for m in (1,2,3)] for k in (1,2,3)]
    inverse = arb_mat([[4*freq[k][m]/7 for k in range(3)] for m in range(3)])
    design = arb_mat([[z-1 for z in row] for row in freq])
    eye = arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    assert all(z.contains(0) for z in (design*inverse-eye).entries())
    logtable = arb_mat(logs)
    coefficients = inverse*logtable
    for k in range(14):
        parity = arb(1+(-1)**k)/2
        co = [parity*((m*k*arb.pi()/7).cos()-1) for m in (1,2,3)]
        out = arb_mat([co])*coefficients
        index = min(k,14-k)//2-1
        for orbit in range(45):
            expected = arb(0) if k==0 or k%2 else logs[index][orbit]
            assert out[0,orbit].overlaps(expected) and out[0,orbit].exp()>0
    minor = arb_mat([[logs[i][j] for j in range(3)] for i in range(3)]).det()
    assert arb("-2.784")<minor<arb("-2.783")
    print("三套对数电导在前三轨道上的三阶子式",minor.str(12))


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
