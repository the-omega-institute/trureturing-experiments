"""六质量三个角的完整闭窗排除任意壳实共同源协方差，允许两部门任意正标尺。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat

N=10
zero=arb_mat(1,N)
class Affine:
    """实数 c + g·u + 误差，其中每个 |u_i| ≤ 1，余项半径为 e。"""
    def __init__(self,c,g=zero,e=0):
        self.c=arb(c); self.g=g; self.e=arb(e)
        self.rad=(sum(abs(v) for v in g.entries())+self.e).upper()
    @staticmethod
    def lift(x):
        return x if isinstance(x,Affine) else Affine(x)
    def bound(self):
        return self.c+arb(0,self.rad)
    def __neg__(self):
        return Affine(-self.c,-self.g,self.e)
    def __add__(self,v):
        v=Affine.lift(v)
        return Affine(self.c+v.c,self.g+v.g,self.e+v.e)
    __radd__=__add__
    def __sub__(self,v):return self+-Affine.lift(v)
    def __rsub__(self,v):return -self+v
    def __mul__(self,v):
        v=Affine.lift(v)
        return Affine(self.c*v.c,self.g*v.c+v.g*self.c,
                 self.e*abs(v.c)+v.e*abs(self.c)+self.rad*v.rad)
    __rmul__=__mul__
    def recip(self):
        low=abs(self.bound()).lower()
        assert low>0
        ci=1/self.c
        return Affine(ci,-self.g*ci**2,ci**2*(self.e+self.rad**2/low))
    def __truediv__(self,v):
        return self*Affine.lift(v).recip()
    def sqrt(self):
        s=self.c.sqrt()
        t=self.bound().sqrt()
        return Affine(s,self.g/(2*s),(self.e+self.rad**2/(t+s)**2)/(2*s))
    def sincos(self):
        s,c=self.c.sin_cos()
        return (Affine(s,self.g*c,self.e*abs(c)+self.rad**2/2),
                Affine(c,-self.g*s,self.e*abs(s)+self.rad**2/2))
class ComplexAffine:
    def __init__(self,r=0,i=0):self.r,self.i=Affine.lift(r),Affine.lift(i)
    def conj(self):return ComplexAffine(self.r,-self.i)
    def __add__(self,v):return ComplexAffine(self.r+v.r,self.i+v.i)
    def __sub__(self,v):return ComplexAffine(self.r-v.r,self.i-v.i)
    def __mul__(self,v):
        if isinstance(v,ComplexAffine): return ComplexAffine(self.r*v.r-self.i*v.i,self.r*v.i+self.i*v.r)
        return ComplexAffine(self.r*v,self.i*v)
    def norm2(self):return self.r*self.r+self.i*self.i
def inputs(values):
    ret=[]
    for i,v in enumerate(values):
        g=arb_mat(1,N);g[0,i]=v.rad().upper()
        ret.append(Affine(v.mid(),g))
    return ret


def corrections(coords):
    """由秩一主项和三次几何级数包围去除坐标相位后的两个极端本征向量。"""
    u=coords[:3];d=coords[3:6]
    mixing=acb_mat(standard(*(v.sin() for v in coords[6:9]),coords[9]))
    out=[]
    for inverse in (False,True):
        row, pivot = (0,1) if inverse else (2,0)
        anchor=mixing[row,pivot].conjugate()
        f=[]
        for i in range(3):
            ratio=(u[row]/u[i] if inverse else u[i]/u[row]) if i!=row else arb(1)
            f.append([acb(1) if (i,j)==(row,pivot) else mixing[i,j].conjugate()/anchor*ratio*
                       ((d[j]/d[pivot] if inverse else d[pivot]/d[j]) if j!=pivot else arb(1)) for j in range(3)])
        # 呈主导的行把对应总 Gram 写成 zz* + E，E 为另两行的 Gram。
        z=acb_mat([[x] for x in f[row]])
        E=lib["gram"]([f[i] for i in range(3) if i!=row])
        es=[1/v if inverse else v for v in eigen]
        lam=(sum(abs(entry)**2 for entry in z.entries())+E.trace().real)*es[0 if inverse else 2]/sum(es)
        err=E.trace().real/lam
        assert 0 <= err < arb(".05")
        errpad=err**4/(1-err)*sum(abs(x)**2 for x in z.entries()).sqrt()
        remainder_matrix=E/lam;term=remainder_matrix*z
        correction=term
        for _ in range(2):
            term=remainder_matrix*term
            correction+=term
        pad=arb(0,errpad.upper())
        correction+=acb_mat([[acb(pad,pad)]]*3)
        t=correction[pivot,0]
        out.append([acb(0) if j==pivot else (correction[j,0]-z[j,0]*t)/(1+t) for j in range(3)])
    return out

def aff_moments(coords):
    """保留共同输入造成的一阶相关性，建立归一化正、逆质量二次矩。"""
    v=inputs(coords); corr=corrections(coords)
    s,a=v[6].sincos();t,b=v[7].sincos();z,c=v[8].sincos();sd,cd=v[9].sincos()
    d=v[3:6]
    p=a*b*z/(s*t)
    den=(1-p)*(1-p)+p*(1-cd)*2
    q=ComplexAffine(1-p*cd,-p*sd)*(den*s*t).recip()
    vtop=[ComplexAffine(1), (ComplexAffine(s/a)-q*(t/a))*(d[0]/d[1]), q*(b*c*d[0]/d[2])]
    vlow=[ComplexAffine(a*d[0]/(s*d[1])),ComplexAffine(1),ComplexAffine(cd,sd)*(z*d[2]/(s*c*d[1]))]
    basis=[]
    for col, delta in zip((vtop,vlow),corr):
        col=[entry+ComplexAffine(ds.real,ds.imag) for entry,ds in zip(col,delta)]
        length=sum(v.norm2() for v in col).sqrt().recip()
        basis.append([v*length for v in col])
    top,low=basis
    mid=[(top[(i+1)%3]*low[(i+2)%3]-top[(i+2)%3]*low[(i+1)%3]).conj() for i in range(3)]
    cols=[low,mid,top]
    diags=([d[0]/d[2],d[1]/d[2],Affine(1)],[Affine(1),d[0]/d[1],d[0]/d[2]])
    out=[]
    for w in diags:
        moment=[[ComplexAffine() for _ in range(3)] for _ in range(3)]
        for k in range(3):
            for j in range(k,3):
                moment[k][j]=sum((cols[k][i].conj()*cols[j][i]*(w[i]*w[i]) for i in range(3)),ComplexAffine())
        out.append(moment)
    return out
def test(coords):
    """复用同一行预条件算子，给任一可能圆上解的六个实坐标逐项夹逼。"""
    rows=lib["coefficients"]([v.mid() for v in coords],eigen,standard,normals,adj)
    mat=arb_mat([row[1:] for row in rows]);rhs=arb_mat([[-row[0]] for row in rows])
    P=mat.mid().inv()
    centre=(P*rhs).mid()
    moments=aff_moments(coords)
    pre_normals=[[acb_mat([[sum(P[i,3*g+k]*normals[g][k][a][b] for k in range(3)) for b in range(3)] for a in range(3)]) for g in range(2)] for i in range(6)]
    rows=[]
    for ns in pre_normals:
        v=[sum(m[a][a].r*n[a,a].real for n,m in zip(ns,moments) for a in range(3))]
        for a,b in ((0,1),(0,2),(1,2)):
            z=sum((m[a][b]*ComplexAffine(n[a,b].real,n[a,b].imag) for n,m in zip(ns,moments)),ComplexAffine())
            v.extend([z.r*2,z.i*(-2)])
        rows.append(v)
    D=[abs((r[0]+sum(x*y for x,y in zip(r[1:],centre.entries()))).bound()).upper() for r in rows]
    error=arb_mat([[int(i==j)-r[j+1].bound() for j in range(6)]for i,r in enumerate(rows)])
    cap=[(1+abs(x)).upper() for x in centre.entries()]
    for _ in range(7):
        shrink=[sum(abs(error[i,j])*cap[j] for j in range(6)).upper() for i in range(6)]
        cap=[min(r, (d+q).upper()) for r,d,q in zip(cap,D,shrink)]
    box=[v+arb(0,r) for v,r in zip(centre.entries(),cap)]
    return any(not (p*p+q*q).contains(1) for p,q in zip(box[::2],box[1::2]))
def all_boxes(sign):
    center=[0.]*9+[sign*1.140625]
    radii=[1.]*9+[.109375]
    priority=[.01,.01,.01,4,2,2,1,2,4,72]
    stack=[(center,radii,0)]
    n=deep=0
    while stack:
        c,r,de=stack.pop()
        b=[arb(v)+arb(0,s) for v,s in zip(c,r)]
        coords=[v+w*x for v,w,x in zip(model["ref"][:9],model["tol"][:9],b[:9],strict=True)]+[b[-1]]
        if test(coords):
            n+=1;deep=max(de,deep)
            continue
        assert de<30
        j=max(range(10),key=lambda i:priority[i]*r[i])
        rr=r.copy();rr[j]/=2
        for sg in (1,-1):
            cc=c.copy();cc[j]+=sg*rr[j];stack.append((cc,rr,de+1))
    print("相位符号、排除终格、最大深度",sign,n,deep)


def independent_check():
    """在两条可达符号分支上，用 T*T 直接求谱核对仿射矩的所有上三角系数。"""
    at_zero = load("quadratic-gram-source")["at_zero"]
    au, ad = at_zero(model, -4, 2345666), at_zero(model, 2, 964890)
    direct, _ = load("two-zero-sharp-quark-fit")["normalized_traces"](au, ad)
    radius = arb(".0574339450674740").union(arb(".0574339450674742"))
    signs = (-1, 1, 1, -1, 1, 1, 1, 1, 1)
    values = [v + s*w*radius for v, s, w in zip(model["ref"][:9], signs, model["tol"][:9], strict=True)]
    delta = load("nine-window-cp-cosine")["phase"](values, direct).acos()
    for sign in (-1, 1):
        point = values + [sign*delta]
        mixing = acb_mat(standard(*(v.sin() for v in point[6:9]), point[9]))
        target = acb_mat([[point[i]*mixing[i,j]/point[3+j] for j in range(3)] for i in range(3)])
        _, frame = model["eig"](adj(target)*target)
        def column(j,p):
            anchor = frame[p,j]
            return [frame[i,j]*anchor.conjugate()/abs(anchor) for i in range(3)]
        low, top = column(0,1), column(2,0)
        middle = [(top[(i+1)%3]*low[(i+2)%3]-top[(i+2)%3]*low[(i+1)%3]).conjugate() for i in range(3)]
        frame = acb_mat([low, middle, top]).transpose()
        for inverse, moment in enumerate(aff_moments(point)):
            down_masses = point[3:6]
            weights = [(down_masses[0]/x)**2 if inverse else (x/down_masses[2])**2 for x in down_masses]
            expected = adj(frame)*acb_mat([[w if i==j else 0 for j in range(3)] for i,w in enumerate(weights)])*frame
            for i in range(3):
                for j in range(i,3):
                    assert expected[i,j].real.overlaps(moment[i][j].r.bound())
                    assert expected[i,j].imag.overlaps(moment[i][j].i.bound())


if __name__ == "__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    def load(name):
        return run_path(str(root / ("ch39-"+name+".py")))
    # 两迹隐函数证书已把全制备盒上任何可达相位包入所用封闭区间。
    load("two-zero-nine-phase-range")["verify"]()
    assert arb(33)/32 < arb("1.0413") < arb("1.2364") < arb(5)/4
    model,up,down=load("two-zero-orbit-condition")["calibration"]()
    adj=model["adj"]
    chart,_=load("shell-pencil-six-mass-window")["shell_frame"](model)
    relative=up*down.inv()
    eigen,basis=model["eig"](adj(relative)*relative)
    eigen[0]=abs(relative.det())**2/(eigen[1]*eigen[2])
    assert 0 < eigen[0] < eigen[1] < eigen[2]
    lib=load("shell-real-ten-observation-window")
    normals=[lib["axes"](adj(chart)*m*basis) for m in (down.inv(),adj(down))]
    standard=model["load"]("shell-clock-ckm-blocks")["standard_matrix"]
    independent_check()
    for sign in (1,-1):
        all_boxes(sign)
