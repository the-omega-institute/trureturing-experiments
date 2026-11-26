"""在整个双角环面及独立软亏额区间上排除两枚实际零点的质子因子统一检出。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx
import numpy as np

SCALE = 10**7
UNIT = 2**20
SIGNS = ((1, 1), (1, -1), (-1, 1), (-1, -1))

def data(root):
    def load(name):
        return run_path(str(root/('ch39-'+name+'.py')))
    geo=load('mincut-real-triangle-tomography')
    sh=load('three-ray-chiral-boundary')
    fixed=load('height-independent-flavor-vacuum')
    model=load('mincut-global-winding')
    angles=load('mincut-reference-ckm-action')
    flavor=load('chiral-neutral-flavor-phase')
    ctx.prec=192
    x=arb('-.01046')
    y=x-angles['reference_angles'](flavor)[0]
    weights=fixed['minimum_order_weights'](model,model['terms'](sh),x,y)
    groups=[(1,0),(0,1),(1,-1)]
    types=load('fixed-word-proton-faces')['signals'](root)[0]
    labels=[[geo['unit'](k) for g,k,s in types if g==group] for group in groups]
    soft = [[bool(k[2]) for k in types if k[0]==group] for group in groups]
    assert list(map(len, labels)) == [14, 14, 4] and max(weights) == 105
    assert min(weights.values()) > 0 and (sum(weights.values())-1).contains(0)
    cap = arb(23)/25*arb(2)**-20/14**2
    assert acb.zeta_zero(1).imag > 14 and sh['soft']['point']['study']['offset'] > 0
    assert cap*sum(weights.values()) < arb('5e-9')
    assert [row.count(True) for row in soft] == [1, 1, 0]
    return weights, labels, [row.index(True) for row in soft[:2]]

def bounds_table(weights,labels,n,offset=.5,radius=1):
    # 对 F 及其导数在每个格心用复球离散傅里叶变换求值，只向外量化。
    shell2=sum(q*q*w for q,w in weights.items())
    shift=acb(0,2*arb.pi()*arb(str(offset))/n).exp()
    st=arb.pi()*radius/n
    val=np.empty((n,len(labels),4),np.int64)
    for j,u in enumerate(labels):
        coef=[acb(0)]*n
        deriv=[acb(0)]*n
        for q,w in weights.items():
            c=w*(u*shift)**q
            coef[(-q)%n]+=c
            deriv[(-q)%n]+=acb(0,q)*c
        v=acb.dft(coef)
        d=acb.dft(deriv)
        for k in range(n):
            g=d[k].real
            # 一项软面的乘子误差 <5e-9，统一扩大函数值。斜率区间仅用于硬面割线。
            a=v[k].real+arb(0,st*abs(g)+shell2*st**2/2+arb('5e-9'))
            b=g+arb(0,st*shell2)
            val[k,j]=[int((a*SCALE).lower().floor().unique_fmpz()),int((a*SCALE).upper().ceil().unique_fmpz()),
                         int((b*SCALE).lower().floor().unique_fmpz()),int((b*SCALE).upper().ceil().unique_fmpz())]
    assert np.max(np.abs(val)) < 10**10 and 4*10**10*UNIT < 2**63
    return val


def intersect(A,B,lo,hi):
    # A,B 为带符号分量的逐面上端，且 t 范围按 UNIT 向外取整。
    for j in range(A.shape[1]):
        a,b=A[:,j],B[:,j]
        impossible=(a<=0)&(b<=0)
        lo=np.where(impossible,UNIT+1,lo)
        lower=(a<=0)&(b>0)
        upper=(a>0)&(b<=0)
        den=np.maximum(1,b-a)
        lo=np.maximum(lo,np.where(lower,(-a*UNIT)//den,0))
        den=np.maximum(1,a-b)
        hi=np.minimum(hi,np.where(upper,(a*UNIT+den-1)//den,UNIT))
    return lo,hi

def signed(a,s):
    return a[:,:,1] if s>0 else -a[:,:,0]

def query(tab,i,j,ss,lo=0,hi=UNIT):
    return intersect(signed(tab[i],ss[0]),signed(tab[j],ss[1]),lo,hi)

def near_cut(tab,i,j,skip):
    n=len(tab)
    distance=(j-i+n//2)%n-n//2
    sel=abs(distance)<=max(3,min(24,n//500))
    bad=np.zeros(len(i),dtype=bool)
    if not np.any(sel):return bad
    x,y,d=i[sel],j[sel],distance[sel]
    # 短弧覆盖两角及中间的全部割线点；只用硬面避免逐高度软质量的差分。
    chosen=np.arange(tab.shape[1])!=skip
    rows=tab[x][:,chosen,:2]
    low=tab[x][:,chosen,2].copy()
    high=tab[x][:,chosen,3].copy()
    for k in range(1,int(max(abs(d)))+1):
        jj=(x+np.sign(d)*np.minimum(k,abs(d)))%n
        low=np.minimum(tab[jj][:,chosen,2],low)
        high=np.maximum(tab[jj][:,chosen,3],high)
    sec=np.stack([low,high],axis=-1)
    can=np.zeros(len(x),dtype=bool)
    for ss in SIGNS:
        lo,hi=intersect(signed(rows,ss[0]),signed(sec,ss[1]),0,UNIT)
        can|=lo<=hi
    bad[sel]=~can
    return bad

def axis_step(prior,tab,ss,skip):
    n=len(tab)
    if prior is None:
        i,j=np.unravel_index(np.arange(n*n),(n,n))
    else:
        a,b=prior[:,:2].T
        i=np.repeat(2*a,4)+np.tile([0,0,1,1],len(a))
        j=np.repeat(2*b,4)+np.tile([0,1,0,1],len(b))
    out=[]
    for k in range(0,len(i),32768):
        a,b=i[k:k+32768],j[k:k+32768]
        lo,hi=query(tab,a,b,ss)
        ok=lo<=hi
        a,b,lo,hi=a[ok],b[ok],lo[ok],hi[ok]
        good=~near_cut(tab,a,b,skip)
        out.append(np.column_stack([a[good],b[good],lo[good],hi[good]]))
    return np.concatenate(out) if out else np.empty((0,4),dtype=np.int64)

def combine(xx,yy,D,ss):
    n=len(D)
    combined=[]
    # 父盒如通过四个钟差面的区间约束，才细分其四维笛卡尔积。
    for x1,x2,lo,hi in xx:
        jp=(yy[:,2]<=hi)&(yy[:,3]>=lo)
        if not np.any(jp):continue
        ys=yy[jp]
        y1,y2,yl,yh=ys.T
        i,j=(x1-y1)%n,(x2-y2)%n
        small,large=query(D,i,j,ss,np.maximum(lo,yl),np.minimum(hi,yh))
        for a,b in zip(y1[small<=large],y2[small<=large]):
            combined.append((x1,x2,a,b))
    return np.asarray(combined,dtype=np.int64).reshape(-1,4)

def refine_coupled(parents, ss, tabs, D, soft_indices):
    n = len(D)
    bits = np.array([[i>>3, (i>>2)&1, (i>>1)&1, i&1] for i in range(16)])
    rows=[]
    for k in range(0,len(parents),2048):
        children=(2*parents[k:k+2048,None,:]+bits[None,:,:]).reshape(-1,4)
        a,b,e,f=children.T
        lo,hi=query(tabs[0],a,b,ss)
        lo,hi=query(tabs[1],e,f,ss,lo,hi)
        lo,hi=query(D,(a-e)%n,(b-f)%n,ss,lo,hi)
        poss=lo<=hi
        keep=children[poss]
        ok=np.ones(len(keep),dtype=bool)
        for g,tab in enumerate(tabs):
            if len(keep):
                pairs,ids=np.unique(keep[:,2*g:2*g+2],axis=0,return_inverse=True)
                ok &=~near_cut(tab,*pairs.T,soft_indices[g])[ids]
        rows.extend(keep[ok])
    return np.asarray(rows,dtype=np.int64).reshape(-1,4)


def verify(root):
    weights, labels, soft_indices = data(root)
    # 这两个角在十四个第一钟面组上可分离，单类凸包一般不足。
    phases = [acb(0,2*arb.pi()*k/16384).exp() for k in (242,19)]
    scores = [[sum(w*(z*u)**q for q,w in weights.items()).real for u in labels[0]] for z in phases]
    margin = min(((a+arb('.939')*b)/arb('1.939')).lower() for a,b in zip(*scores))-arb('5e-9')
    assert margin > arb('.002')
    print("仅十四个第一钟组的严格分离反例余量",margin.str(13),flush=True)
    states = {(g, ss):None for g in (0,1) for ss in SIGNS}
    for n in (256,512,1024,2048,4096,8192):
        tabs = [bounds_table(weights, group, n) for group in labels[:2]]
        for g,tab in enumerate(tabs):
            for ss in SIGNS:
                states[g,ss] = axis_step(states[g,ss],tab,ss,soft_indices[g])
        print("每圆分格数、两个单钟类的四种符号余盒数",n,
              [[len(states[g,ss]) for ss in SIGNS] for g in (0,1)],flush=True)
    # 两个格心之差落在整数格点，对应角差的半宽为单格半宽的两倍。
    D = bounds_table(weights,labels[2],n,offset=0,radius=2)
    coupled = {ss:combine(states[0,ss],states[1,ss],D,ss) for ss in SIGNS}
    print("每圆分格数及四维余盒数",n,[len(parents) for parents in coupled.values()],flush=True)
    for n in (16384,32768,65536):
        tabs = [bounds_table(weights,group,n) for group in labels[:2]]
        D = bounds_table(weights,labels[2],n,offset=0,radius=2)
        coupled = {ss:refine_coupled(parents,ss,tabs,D,soft_indices) for ss,parents in coupled.items()}
        print("每圆分格数及四维余盒数",n,[len(parents) for parents in coupled.values()],flush=True)
    assert not any(len(rows) for rows in coupled.values())
    print("全部两点与独立软亏额的严格分离方向均已排除。")


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
