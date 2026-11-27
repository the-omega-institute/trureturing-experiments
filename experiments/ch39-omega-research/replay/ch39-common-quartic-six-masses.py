"""在同一编号内，用同一单位区间四次有理函数与同一复超平面精确重现六个质量。"""
from pathlib import Path
from runpy import run_path
from math import prod

from flint import acb, acb_mat, arb, arb_mat, arb_series, ctx


FLOOR = "1e-12"


CASES = (
    (".99999",".0018186078866516415","19.651365591150167","12.910548347316272",
     "19.0882494660982",".2469313030911994"),
    (".8941438363071211",".01035140019469654","19.026242037122625","18.754637443352735",
     "19.262487702375","19.080462544049897"),
)


def read(coeff,peaks,nodes,rows,targets,*,matrices=False):
    amplitude,scale,z1,z2 = coeff[:4]
    spread = coeff[4] if len(coeff)>4 else arb(0)
    residues,axes,operators = [],[],[]
    for ts,projections,mu in zip(nodes,rows,targets,strict=True):
        weights = []
        for t in ts:
            u = scale**2*((t-z1)**2+spread**2)*(t-z2)**2
            v = (t-peaks[0])*(t-peaks[1])
            weights.append(arb(FLOOR)+(amplitude-arb(FLOOR))*u/(u+v*v))
        p,h,re,im,plus,minus = [sum(w*c for w,c in zip(weights,row,strict=True)) for row in projections]
        axis = (2*re,-2*im,p-h)
        mag = sum(t*t for t in axis).sqrt()
        axes.append([t/mag for t in axis])
        eig = [(p+h-mag)/2,(p+h+mag)/2,plus,minus]
        residues.append([prod(l-x for x in mu)/prod(l-v for j,v in enumerate(eig) if j!=i)
                         for i,l in enumerate(eig)])
        if matrices:
            operators.append(acb_mat([[p,acb(re,im),0,0],[acb(re,-im),h,0,0],
                                      [0,0,plus,0],[0,0,0,minus]]))
    if matrices:
        return residues,axes,operators
    return residues,sum(a*b for a,b in zip(*axes,strict=True))


def constraints(coeff,points,nodes,rows,target):
    f,der = [],[]
    for j in range(2):
        var = [arb_series([x,int(i==j)]) for i,x in enumerate(points)]
        r,_ = read(coeff,var,nodes,rows,target)
        v = [a-b for a,b in zip(r[0][2:],r[1][2:],strict=True)]
        f = [x[0] for x in v]
        der.append([x[1] for x in v])
    return arb_mat([[x] for x in f]),arb_mat(der).transpose()


def ckm_branches(coeff,box,nodes,rows,masses,*,verbose=True):
    residues,axes,operators = read(coeff,box,nodes,rows,masses,matrices=True)
    a,b = axes
    c = sum(u*v for u,v in zip(a,b,strict=True))
    s = (1-c*c).sqrt()
    u = [(v-c*z)/s for z,v in zip(a,b,strict=True)]
    v = [a[(i+1)%3]*u[(i+2)%3]-a[(i+2)%3]*u[(i+1)%3] for i in range(3)]
    dx,dy = [(r[1]-r[0])/(r[1]+r[0]) for r in residues]
    coord = (dy-c*dx)/s
    delta = (1-dx*dx-coord*coord).sqrt()
    ru = residues[0]
    res = []
    ident = acb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    for sign in (-1,1):
        h = [dx*i+coord*j+sign*delta*k for i,j,k in zip(a,u,v,strict=True)]
        f = ((1+h[2])/2).sqrt()
        amp = (ru[0]+ru[1]).sqrt()
        normal = acb_mat([[f*amp],[acb(h[0],h[1])*amp/(2*f)],[ru[2].sqrt()],[ru[3].sqrt()]])
        state = []
        for op,mu in zip(operators,masses,strict=True):
            vs = []
            for mass in mu:
                vector = (op-ident*mass).solve(normal)
                adj = vector.conjugate().transpose()
                vs.append(vector*(1/(adj*vector)[0,0].real.sqrt()))
            state.append(acb_mat([[col[i,0] for col in vs] for i in range(4)]))
        mix = state[0].conjugate().transpose()*state[1]
        j = (mix[0,0]*mix[1,1]*mix[0,1].conjugate()*mix[1,0].conjugate()).imag
        res.append((abs(mix[0,2]),j,[[abs(mix[i,k])**2 for k in range(3)] for i in range(3)]))
        if verbose:
            print("镜像标记、|Vub|、有符号 J",sign,abs(mix[0,2]).str(12),j.str(12),flush=True)
    assert all(res[0][2][i][k].overlaps(res[1][2][i][k]) for i in range(3) for k in range(3))
    assert res[0][1].overlaps(-res[1][1]) and not res[0][1].contains(0)
    return res


def reference_windows():
    angles=[(arb(c)-arb(e)).union(arb(c)+arb(e)) for c,e in
            ((".22704",".00082"),(".04275",".00042"),(".00377",".00008"),("1.139",".023"))]
    a,b,g,phase=angles
    req_vub=g.sin()
    req_j=a.sin()*a.cos()*b.sin()*b.cos()*g.sin()*g.cos()**2*phase.sin()
    return req_vub,req_j



def certified_cases(root, *, radius_text="1e-15"):
    native=run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ordered=run_path(str(root/"ch39-bounded-common-rational-amplitude.py"))["ordered"]
    for index,case in enumerate(CASES,2_306_888_100_085):
        ps=native["native"](root,index)
        ctx.prec=288
        ctx.cap=2
        target=[list(map(arb,row)) for row in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
        h=acb.zeta_zero(index).imag
        nodes=[ordered(shell,q,h) for q in (-4,2)]
        rows=[]
        for sect in ps:
            coord=[(p[0,0].real,p[1,1].real,p[0,1].real,p[0,1].imag,
                    p[2,2].real+p[2,3].imag,p[2,2].real-p[2,3].imag) for p in sect]
            rows.append(list(zip(*coord,strict=True)))
        pars=list(map(arb,case));coeff,center=pars[:4],pars[4:]
        assert 0<coeff[0]<1 and coeff[1]>0
        for _ in range(4):
            f,j=constraints(coeff,center,nodes,rows,target)
            corrected=j.solve(f)
            center=[(x-v).mid() for x,v in zip(center,corrected.entries(),strict=True)]
        radius=arb(radius_text)
        box=[x+arb(0,radius.upper()) for x in center]
        assert all(not (h-z).contains(0) for h in box for z in coeff[2:])
        f,j=constraints(coeff,center,nodes,rows,target)
        inv=j.inv()
        _,all_j=constraints(coeff,box,nodes,rows,target)
        defect=arb_mat([[1,0],[0,1]])-inv*all_j
        gain=max(sum(abs(z).upper() for z in row) for row in defect.tolist())
        assert max(abs(z).upper() for z in (inv*f).entries())+gain*radius<radius
        r,c=read(coeff,box,nodes,rows,target)
        assert all(0<v<1 for row in r for v in row)
        sx,sy=[(x[1]-x[0])/(x[1]+x[0]) for x in r]
        margin=(1-c*c)*(1-sx*sx)-(sy-c*sx)**2
        assert margin>0
        yield index,coeff,box,nodes,rows,target,gain,margin


def verify(root):
    for index,coeff,box,nodes,rows,target,gain,margin in certified_cases(root):
        branches = ckm_branches(coeff,box,nodes,rows,target)
        req_vub,req_j = reference_windows()
        assert all(not x.overlaps(req_vub) and not j.overlaps(req_j) for x,j,_ in branches)
        print("编号、两个峰位、压缩映射范数上界、共同法向相容余量",index,
              *[x.str(15) for x in box],gain.str(8),margin.str(10),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
