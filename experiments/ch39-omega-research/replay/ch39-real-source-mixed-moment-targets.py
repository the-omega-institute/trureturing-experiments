"""质量窗和三个混合角窗上的两个混合矩。相位可在整圆上任取。"""
from itertools import product
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


def normal_masses(masses, inverse):
    # 每项都直接约去相同部门的标尺。
    exponent = arb(-2 if inverse else 2)/3
    return [(masses[i]**2/(masses[(i+1)%3]*masses[(i+2)%3]))**exponent
            for i in range(3)]


def kernel(x, y, z, cosine):
    # x,y,z 是三个角的正弦平方，以稳定的正数式给出两个接近零的对角元。
    inter = 2*(x*(1-x)*y*(1-y)*z).sqrt()*cosine
    return [[x+(1-x)*z, 1-x*(1-z), 1-z],
            [1-x*(1-y)-(1-x)*y*z-inter, x+y-x*y-x*y*z+inter, 1-y*(1-z)],
            [1-x*y-(1-x)*(1-y)*z+inter, 1-(1-x)*y-x*(1-y)*z-inter, y+(1-y)*z]]


def moment(up, down, mix):
    return sum(up[i]*down[j]*mix[i][j] for i,j in product(range(3),repeat=2))


def angle_derivatives(x, y, z, cosine):
    # 三个偏导仅对 sin²(theta) 求导。相位干涉的符号与 kernel 保持同序。
    a = (x*(1-x)*y*(1-y)*z).sqrt()*cosine
    ix = a*(1-2*x)/(x*(1-x))
    iy = a*(1-2*y)/(y*(1-y))
    iz = a/z
    return [
        [[1-z,-(1-z),0],[-(1-y)+y*z-ix,1-y-y*z+ix,0],
         [-y+(1-y)*z+ix,y-(1-y)*z-ix,0]],
        [[0,0,0],[x-(1-x)*z-iy,1-x-x*z+iy,-(1-z)],
         [-x+(1-x)*z+iy,-(1-x)+x*z-iy,1-z]],
        [[1-x,x,-1],[-(1-x)*y-iz,-x*y+iz,y],
         [-(1-x)*(1-y)+iz,-x*(1-y)-iz,1-y]],
    ]


def target_window(model, inverse=False, free_phase=True, normalized=False):
    """返回混合矩最小值、最大值的球区间。传入原质量模型中的 ref、tol。"""
    centers, widths = model["ref"], model["tol"]
    pairs = [(c-w,c+w) for c,w in zip(centers,widths,strict=True)]
    full_mass = [lo.union(hi) for lo,hi in pairs[:6]]
    angles = [(lo.sin()**2,hi.sin()**2) for lo,hi in pairs[6:9]]
    assert all(0<lo<hi<arb.pi()/2 for lo,hi in pairs[6:9])
    phase = (arb(-1),arb(1)) if free_phase else (pairs[9][1].cos(), pairs[9][0].cos())
    assert free_phase or 0 < pairs[9][0] < pairs[9][1] < arb.pi()
    x,y,z = [lo.union(hi) for lo,hi in angles]
    cosine = phase[0].union(phase[1])
    w = normal_masses(full_mass[:3], inverse)
    t = normal_masses(full_mass[3:], inverse)
    mix = kernel(x,y,z,cosine)
    total = moment(w,t,mix)
    signs = []
    rows = [sum(t[j]*mix[i][j] for j in range(3)) for i in range(3)]
    cols = [sum(w[i]*mix[i][j] for i in range(3)) for j in range(3)]
    for spectrum, readings in ((w,rows),(t,cols)):
        for i in range(3):
            if normalized:
                # 去掉迹积时，其偏导为该行读数与加权均值之差，消去自项。
                derivative = sum(spectrum[k]*(readings[i]-readings[k]) for k in range(3) if k!=i)
            else:
                derivative = 3*spectrum[i]*readings[i]-total
            signs.append((1-2*int(inverse))*derivative)
    signs += [moment(w,t,d) for d in angle_derivatives(x,y,z,cosine)]
    # d/d cos(delta) 的非零公共因子严格为负。
    signs.append(-(w[1]-w[2])*(t[0]-t[1]))
    assert all(s>0 or s<0 for s in signs), (inverse,free_phase,signs)
    pairs = pairs[:6]+angles+[phase]
    endpoints = []
    for maximum in (False,True):
        pt = [p[int((s>0)==maximum)] for p,s in zip(pairs,signs,strict=True)]
        up, down = normal_masses(pt[:3],inverse), normal_masses(pt[3:6],inverse)
        value = moment(up,down,kernel(*pt[6:]))
        endpoints.append(value/(sum(up)*sum(down)) if normalized else value)
    return tuple(endpoints)


def mixed_windows(model, free_phase=True, normalized=False):
    """返回 (原矩的下上界, 逆矩的下上界)，均为向外舍入的端点。"""
    limits = [target_window(model, inverse, free_phase, normalized) for inverse in (False,True)]
    return [(m.lower(),M.upper()) for m,M in limits]


def cross_rows(up, down):
    """混合矩的 9×3 系数矩阵，每行是上下部门各一行的双线性叉积。"""
    return acb_mat([[up[i,(k+1)%3]*down[j,(k+2)%3]-up[i,(k+2)%3]*down[j,(k+1)%3]
                     for k in range(3)] for i,j in product(range(3),repeat=2)])


def verify(replay):
    model = run_path(str(replay/"ch39-two-zero-quark-window.py"))
    ctx.prec = 224
    for free in (True,False):
        windows = mixed_windows(model,free)
        print("相位", "自由" if free else "参考窗", [tuple(x.str(19) for x in p) for p in windows])
        normalized = mixed_windows(model,free,normalized=True)
        print("除以两部门迹积", [tuple(x.str(15) for x in p) for p in normalized])
    # 独立以标准复 CKM 计算全部模方，在参考中心投影到同一组读数。
    standard = run_path(str(replay/"ch39-shell-clock-ckm-blocks.py"))["standard_matrix"]
    centers = model["ref"]
    v = standard(*(x.sin() for x in centers[6:9]),centers[9])
    stable = kernel(*(x.sin()**2 for x in centers[6:9]),centers[9].cos())
    for i,j in product(range(3),repeat=2):
        assert stable[i][j].overlaps(1-abs(v[i][j])**2)
    for inverse,(lo,hi) in zip((False,True),mixed_windows(model),strict=True):
        s,d = [normal_masses(centers[i:i+3],inverse) for i in (0,3)]
        d1 = acb_mat([[s[i]*int(i==j) for j in range(3)] for i in range(3)])
        d2 = acb_mat(v)*acb_mat([[d[i]*int(i==j) for j in range(3)] for i in range(3)])*acb_mat(v).conjugate().transpose()
        raw = (d1.trace()*d2.trace()-(d1*d2).trace()).real
        assert raw.overlaps(moment(s,d,stable)) and lo < raw < hi
        print("逆质量" if inverse else "正质量","参考中心",raw.str(18))

    # 用含复相位的两矩阵和非坐标轴实谱架复核观测端与源端的权重方向。
    up = acb_mat([[1,acb(1,2),-1],[0,2,acb(0,1)],[0,0,arb(1)/2]])
    down = acb_mat([[1,0,0],[acb(2,1),1,0],[1,acb(0,-1),1]])
    frame = acb_mat([[3,4,0],[-4,3,0],[0,0,5]])/5
    spectrum = [arb(3),arb(2),arb(1)/6]
    source = frame*acb_mat([[spectrum[i].sqrt()*int(i==j) for j in range(3)]
                            for i in range(3)])*frame.transpose()
    adj = lambda v: v.conjugate().transpose()
    hu, hd = [adj(m*source)*(m*source) for m in (up,down)]
    for inverse in (False,True):
        a,b = (adj(up.inv()),adj(down.inv())) if inverse else (up,down)
        gu,gd = (hu.inv(),hd.inv()) if inverse else (hu,hd)
        measured = (gu.trace()*gd.trace()-(gu*gd).trace()).real
        cols = cross_rows(a,b)*frame
        calculated = sum(abs(cols[j,i])**2*(spectrum[i] if inverse else 1/spectrum[i])
                         for i,j in product(range(3),range(9)))
        assert calculated.overlaps(measured)
    print("两个混合矩与复行叉积的正逆谱权恒等式通过")


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
