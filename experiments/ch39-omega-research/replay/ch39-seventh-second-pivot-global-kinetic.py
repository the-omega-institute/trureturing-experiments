"""第七壳原生第二枢轴的五个荷六端点条件与三处色动能的全域冲突。"""
from pathlib import Path
from runpy import run_path
from flint import arb,ctx
import numpy as np


def build_table(root):
    lib,source,ws,(mom,heat,field),data,rows=run_path(str(root/'ch39-fifth-arbitrary-pivot-obstructions.py'))['inputs'](root,number=2306888100092)
    core=lib('fourth-third-pivot-global-kinetic-obstruction')
    lo,hi,_,_,start,_=core['build_grid'](root,(0,),number=2306888100092,energies=(10**11,10**9,10**7))
    ticks=core['cutoff_ticks']();p0=data['p0'];pre=core['SCALE'];sc=10**10
    def bracket(v):
     return [int((v.lower()*sc).floor().unique_fmpz()),int((v.upper()*sc).ceil().unique_fmpz())]
    c=ws[1];stage=[int(j) for j in ticks]
    pos=[p0*(arb(j)/512).exp() for j in stage]
    ctx.prec=160
    assert pos[0]<source[0]<source[-1]<pos[-1]
    real=[j for j,t in enumerate(pos) if source[0]<t<source[-1]]
    coords=[source[0],*[pos[j] for j in real],source[-1]]
    kk=sum(abs(v) for _,v,_ in field)/(2*arb.pi())
    ww=2*sum(abs(g) for g,_,_ in field)/(4*arb.pi())
    cols=[]
    # 角核的每个节点已含左右弦余项。将首末两段截到源谱允许的端点。
    inds=[0,1,2,4,6,8];factor=[12,3,arb('.01'),10,10,10]
    for i,x in enumerate(coords):
     if i==0:j=real[0]-1
     elif i==len(coords)-1:j=real[-1]
     else:j=real[i-1]
     lam=(x/pos[j]).log()/(pos[j+1]/pos[j]).log()
     left=[arb(int(a)).union(arb(int(b)))/pre for a,b in zip(lo[j,inds],hi[j,inds])]
     right=[arb(int(a)).union(arb(int(b)))/pre for a,b in zip(lo[j+1,inds],hi[j+1,inds])]
     gauge=[((1-lam)*a+lam*b)*f for a,b,f in zip(left,right,factor)]
     if i+1<len(coords):upper=coords[i+1];rr=(upper/x).log()
     else:upper=x;rr=arb(0)
     if i>0:rr=max(rr.upper(),(x/coords[i-1]).log().upper())
     th=mom(field,heat,x)
     ex=[arb(0, (v*rr**2/8).upper()) for v in [kk*upper,ww*upper**2,upper.sqrt()/4]]
     cols.append([*gauge[:2],th[0]+ex[0],th[1]+ex[1],gauge[2],x.sqrt()+ex[2],*gauge[3:]])
    basis=[start[0]*12,start[1]*3,-arb(1),-arb('1e-122'),arb(0),-sum(x.sqrt() for x in source),start[4]*10,start[6]*10,start[8]*10]
    bands=np.array([[bracket(v) for v in col] for col in cols],dtype=np.int64)
    table=dict(bands=bands,sc=sc,c=np.array([bracket(w) for w in c]),b=np.array([bracket(v) for v in basis]))
    return table,coords

def certify(table,tree):
    bands=table['bands'].astype(object); C=table['c'].astype(object);B=table['b'].astype(object)
    sc=int(table['sc']);N=len(bands)-1;results=[];colour_losses=[]
    def visit(node,domain):
        if len(node)==4:
            i,m,l,r=node;i-=1
            a,b=domain[i];assert a<m<b
            left=domain.copy();right=domain.copy();left[i]=(a,m);right[i]=(m,b)
            visit(l,left);visit(r,right);return
        lam=np.array(node,dtype=object)
        assert len(lam)==9 and all(v>=0 for v in lam[6:])
        signed=np.where(lam>=0,1,0)
        total=sum(B[j,signed[j]]*lam[j] for j in range(9))*sc
        tr=lam[5];lam[5]=0
        lower=bands[:,np.arange(9),1-signed]@lam
        upper=bands[:,np.arange(9),signed]@lam
        trace=bands[:,5,1 if tr>=0 else 0]*tr*sc
        for i,(a,b) in enumerate(domain):
            l,h=C[i];assert l*h>0
            q=upper[a:b+1] if l>0 else lower[a:b+1]
            total+=max(np.maximum(l*q,h*q)+trace[a:b+1])
        assert total<0,total
        results.append(-total)
        coefficient=sum(node[6:])
        if coefficient:
            colour_losses.append(arb(-total)/(coefficient*sc**2))
    visit(tree,[(0,N)]*19)
    return len(results),arb(min(results))/(10**7*sc**2),min(v.lower() for v in colour_losses)


# 中间结点依次记录层号、公共截点序号、左右子树，叶记录九个整数乘子。
TREE=[16, 1504, [18, 1472, [70447, -400000, 277987, 107462, 4205317, 0, 0, 0, 4249], [18, 2099, [71495, -400000, 359275, 138885, 4747364, 0, 0, 9951, 0], [-100000, -400000, 10000000, 3643748, 100000000, 63611, 200000, 200000, 200000]]], [1, 1536, [18, 1472, [70564, -400000, 368357, 142397, 4205081, 0, 0, 0, 4624], [1, 768, [18, 2099, [73261, -400000, 257322, 99474, 4286781, 0, 14690, 0, 1], [8296, -45859, 32772, 11714, 34145837, 5517, 0, 200000, 200000]], [4, -23, 19, 7, 100000000, 0, 0, 0, 2]]], [2, -12, 11, 4, 100000000, 0, 0, 0, 1]]]

def verify(root):
    table,coords=build_table(root)
    leaves,deficit,colour_loss=certify(table,TREE)
    assert leaves==8 and deficit>arb(".00261") and colour_loss>arb("2.6238")
    # 已证端点正例的原解盒，在 10^11 GeV 单独核对色逆耦合。
    fixed=run_path(str(root/"ch39-seventh-second-pivot-endpoint.py"))["ANCHORS"]
    coeff=table["c"];sc=int(table["sc"]);total=int(table["b"][6,1])*sc
    for i,x in enumerate(map(arb,fixed)):
        width=arb("1e-9") if i==0 else (arb("1e-5") if i in (9,17,16,15,18) else arb(0))
        box=x*(1+arb(0,width.upper()))
        j=next(j for j,y in enumerate(coords) if not y<box.lower())
        k=next(j for j,y in enumerate(coords) if y>box.upper())
        data=table["bands"][max(0,j-1):k+1,6,:]
        low,high=int(data[:,0].min()),int(data[:,1].max())
        total+=max(a*b for a in (low,high) for b in map(int,coeff[i]))
    colour=arb(total)/sc**2
    print("端点正例在十的十一次方 GeV 的色逆耦合上界",colour.str(14))
    print("全部连续分支、最小严格分离余量",leaves,deficit.str(14))
    print("全域三个色读数中至少一处的统一亏量",colour_loss.str(14))
    assert colour<0
    return colour,deficit


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
