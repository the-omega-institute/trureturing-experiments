"""保持 357 条格点时，三枚零点之间的连续重谱与正交汇解支。"""
from concurrent.futures import ProcessPoolExecutor
from math import factorial
from os import cpu_count
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parent
ORDER = 48
S = arb(2)**-17
DIAGONAL = [arb(2)**-6, arb(2)**-3, arb(2)**-4, arb(1)]
RADIUS = arb(4)/125
NS = 1024


def vector(values):
    return arb_mat([[v] for v in values])


def columns(vectors):
    return arb_mat([v.entries() for v in vectors]).transpose()


def midpoint(matrix):
    return arb_mat([[matrix[i,j].mid() for j in range(matrix.ncols())] for i in range(matrix.nrows())])


def basis(x, derivative=0):
    """依次取 r²,r,1,r^-1,...,r^-ORDER, log r；导数对 r=x/S。"""
    r = x/S
    powers = [r**2, r, arb(1)]
    for j in range(ORDER+derivative):
        powers.append(powers[-1]/r)
    if derivative == 0:
        return vector([*powers[:ORDER+3],r.log()])
    if derivative == 1:
        return vector([(2-i)*powers[i+1] for i in range(ORDER+3)]+[1/r])
    return vector([(2-i)*(1-i)*powers[i+2] for i in range(ORDER+3)]+[-1/r**2])


def spectrum(shell, lib, neutral, h, half):
    ordered = lib['chiral']['ordered']
    kappa = arb(23)/25*arb(2)**-20
    offset = shell['soft']['point']['study']['offset']
    delta = kappa/(h**2+offset)
    band = h+arb(0,half.upper())
    deltas = kappa/(band**2+offset)
    pole = shell['soft']['pole']
    speed2 = sum(2*(shell['geo']['theta'](d)-shell['geo']['theta'](244))**2
                 for i,j,d,_ in shell['geo']['edges'] if pole not in (i,j) and d in (724,804))
    exact, intervals = {}, {}
    for q in (0,1,4,2,3,6):
        if q == 0:
            levels = neutral(shell,delta)
            enclosed = neutral(shell,deltas)
        else:
            eig = shell['laplacian'](q,h).eig(multiple=True)
            assert len(eig)==19 and all(z.imag.contains(0) for z in eig)
            levels = [z.real for z in eig]
            margin = 4*abs(deltas-delta)+q*half*speed2.sqrt()
            enclosed = [v+arb(0,margin.upper()) for v in levels]
        for result,values,loss in ((exact,levels,delta),(intervals,enclosed,deltas)):
            result[q] = ordered([v*abs(v-15+loss/3) for v in values])
    point, boxes = lib['fields'](exact),lib['fields'](intervals)
    assert len(point)==len(boxes)==99
    assert all(z.is_zero() if i==0 else z>0 for i,(_,_,z) in enumerate(boxes))
    return point,boxes


def series(rows, gauge):
    """三类热核 xE₂、x²E₃、E₁ 的截断展开系数矩阵。"""
    mat=arb_mat(5,ORDER+4)
    pi4=4*arb.pi()
    for (g,a,z),bs in zip(rows,gauge,strict=True):
        y=z/S
        log=y.log()+arb.const_euler() if not y.is_zero() else arb(0)
        v=[2*a/pi4,-g/pi4,*[-b/pi4 for b in bs]]
        # 除公共因子外 (场 × 读数 × 谱矩) 的非多项式部分。
        mat[0,1]+=v[0]*S
        mat[0,2]+=v[0]*S*y*(log-1)
        mat[0,ORDER+3]-=v[0]*S*y
        mat[1,0]+=v[1]*S**2/2
        mat[1,1]-=v[1]*S**2*y
        mat[1,2]+=v[1]*S**2*y**2*(arb(3)/2-log)/2
        mat[1,ORDER+3]+=v[1]*S**2*y**2/2
        for j in range(2,5):
            mat[j,2]-=v[j]*log
            mat[j,ORDER+3]+=v[j]
        # 常规幂项按 r^-n 收集，各行恰截到 n=ORDER。
        e=arb(1)
        for n in range(1,ORDER+1):
            e *= -y/n
            # E1 的 r^-n 项为 -e/n，xE2 项为 S*y*e/[n(n+1)]，
            # x²E3 项为 -S²*y²*e/[n(n+1)(n+2)]。
            coeff=[S*y*e/(n*(n+1)), -S**2*y**2*e/(n*(n+1)*(n+2)), -e/n]
            for j in range(5):
                mat[j,n+2]+=v[j]*coeff[min(j,2)]
    return mat


def physical_derivative(z, m):
    """在一个质量区间内对一个有符号谱矩向量求三项核的 z 导数。"""
    y=z/S
    log=y.log()+arb.const_euler()
    k=log*m[2]-m[ORDER+3]
    w=-S*m[1]-S*y*(log-1)*m[2]+S*y*m[ORDER+3]
    u=-m[2]/z
    factor=arb(1)
    for n in range(1,ORDER+1):
        factor *= -y/n
        term=factor*m[n+2]
        k+=term/n
        w-=S*y*term/(n*(n+1))
        u-=term/z
    return k,w,u


def verify_cell(index, setup):
    shell,lib,neutral,edges,labels,fields,initial,beta,logscale,seed=setup
    assert len(edges)==len(labels)==361
    lo,hi=arb('1310310.5534'),arb('1310311.7259')
    half=(hi-lo)/(2*NS)
    h=lo+(2*index+1)*half
    rows, boxes=spectrum(shell,lib,neutral,h,half)
    fields=[(0,0,0)]+fields
    coeff=series(rows,fields)
    base=vector([0,0,*[v-b*logscale for v,b in zip(initial,beta,strict=True)]])
    target=vector([1,arb('1e-122'),0,0])
    proj=arb_mat([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,-1,0],[0,0,1,0,-1]])
    scale=arb(7)/1000
    def line(i,t,d=0):
        a,b,w=edges[i]
        return basis(scale*(a+t*(b-a)),d)*(w*(scale*(b-a)/S)**d)
    const=vector([0]*(ORDER+4))
    for i,digit in enumerate(labels):
        if i not in (17,338):
            const+=line(i,1 if i in (0,16) else arb(int(digit,12))/11)
    d0,d16=[line(i,0)-line(i,1) for i in (0,16)]
    def measure(p):
        return const+d0*p[0]+d16*p[1]+line(17,p[2])+line(338,p[3])
    def deriv(p):
        return columns([d0,d16,line(17,p[2],1),line(338,p[3],1)])
    p=[arb(v) for v in seed]
    for _ in range(2):
        step=(proj*coeff*deriv(p)).solve(proj*(base+coeff*measure(p))-target)
        p=[(v-step[i,0]).mid() for i,v in enumerate(p)]
    m=measure(p)
    dp=deriv(p)
    j=coeff*dp
    pre=midpoint(midpoint(proj*j).inv())
    assert not pre.det().contains(0)
    normal=arb_mat([[pre[i,k]/DIAGONAL[i] for k in range(4)] for i in range(4)])
    a=normal*proj
    # 在根处 U₃−LF=U₃。择 L 使中心处的一阶参数导数接近零。
    L=midpoint(arb_mat([[j[4,k] for k in range(4)]])*pre)
    free=arb_mat([[0,0,0,0,1]])-L*proj
    transform=arb_mat([[a[i,k] for k in range(5)] for i in range(4)]+[free.entries()])
    center=base+coeff*m
    rv=normal*(proj*center-target)
    cell=list(p)
    for i in range(4):
        cell[i]=p[i]+arb(0,(RADIUS*DIAGONAL[i]).upper())
    assert all(0<t<1 for t in cell)
    span=deriv(cell)
    # 在位置中心展开雅可比列，一次均值界避免不同幂的区间失去相关。
    dm=[[span[k,i] for k in range(span.nrows())] for i in range(4)]
    disps=[cell[i]-p[i] for i in range(4)]
    J=transform*j
    for i,ix in ((2,17),(3,338)):
        q=transform*coeff*line(ix,cell[i],2)*disps[i]
        for k in range(5):
            J[k,i]+=q[k,0]
    residual=[rv[i,0] for i in range(4)]+[(free*center+L*target)[0,0]+beta[2]*logscale]
    # 逐场将三个核的质量导数先作同一线性组合，再乘该场的质量偏移。
    # 对 E1 的指数常数先乘总权，避免误差乘上单线绝对权之和。
    for (g,alpha,z),(_,_,zb),bs in zip(rows[1:],boxes[1:],fields[1:],strict=True):
        dz=zb-z
        multipliers=transform*arb_mat([[2*alpha,0,0],[0,-g,0],*[ [0,0,-b] for b in bs]])/(4*arb.pi())
        for jcol,weights in enumerate([m.entries(),*dm]):
            kernels=vector(physical_derivative(zb,weights))
            v=multipliers*kernels*dz
            for k in range(5):
                if jcol==0: residual[k]+=v[k,0]
                else: J[k,jcol-1]+=v[k,0]
    # 给真正热核加上截断尾。|z/x|≤5，保留的所有对数项均为精确项。
    # 第 n 项的三个系数绝对值分别不超过 z、z²、1 乘 (z/x)^n/n!。
    # 位置导数再乘 n |x'/x|；对 n>=ORDER+1，两个级数的几何比都小于 1/2。
    floor=min((scale*a).lower() for a,_,_ in edges)
    zmax=max(z.upper() for _,_,z in boxes)
    assert zmax<1 and zmax/floor < 5 and all(0<a<=b and b<=1 for a,b,_ in edges)
    tail=2*arb(5)**(ORDER+1)/factorial(ORDER+1)
    tv=sum(abs(w) for _,_,w in edges)
    row_bound=sum(max(abs(2*alpha),abs(g),*(abs(b) for b in bs))
                  for (g,alpha,_),bs in zip(boxes,fields,strict=True))/(4*arb.pi())
    for k in range(5):
        scaling=sum(abs(transform[k,j]) for j in range(5))
        eps=(tail*tv*scaling*row_bound).upper()
        residual[k]+=arb(0,eps)
        for j in range(4):
            d=2 if j<2 else (ORDER+1)*scale/floor
            J[k,j]+=arb(0,(eps*d).upper())
    contraction=max(sum(abs(int(i==j)-J[i,j]*DIAGONAL[j]) for j in range(4)).upper()
                    for i in range(4))
    displacement=max(abs(x).upper() for x in residual[:4])
    assert beta[2]<beta[1]<0<beta[0] and logscale>0
    assert contraction<1 and displacement+RADIUS*contraction<RADIUS,(index,contraction,displacement)
    root_radius=(displacement/(1-contraction)).upper()
    floor_u=(residual[4]-root_radius*sum(abs(J[4,j]*DIAGONAL[j]) for j in range(4))).lower()
    assert floor_u>arb('1.49'), (index,floor_u)
    return p,displacement,contraction,root_radius,floor_u


def setup():
    def load(name):
        return run_path(str(ROOT/f'ch39-{name}.py'))
    edges,_,_,_=load('two-copy-gram-cutoff-two-exceptions')['model_inputs']()
    lib=load('chiral-minimal-neutral-gravity')
    model,_,_,initial,_=lib['chiral']['source']['mean']['inputs']()
    shell=load('three-ray-chiral-boundary')
    labels=''.join(load('two-copy-four-split-lines')['LABELS'])
    fields=load('signed-gravity-unification')['gauge_weights'](model)
    ctx.prec=256
    return (shell,lib,load('two-zero-cutoff-gravity')['neutral_levels'],edges,labels,
            fields,initial,[b/(2*arb.pi()) for b in model['beta']],(arb(10**16)/model['mz']).log(),
            load('two-copy-two-split-lines')['SEED'])


def certify_block(interval):
    data=setup()
    start,stop=interval
    results=[]
    for i in range(start,stop):
        p,*bounds=verify_cell(i,data)
        # 跨进程只传递向外取整的数值串，主进程重新解析为包含原值的球。
        results.append([v.str(35) for v in [*p,*bounds]])
    return results


def verify():
    ctx.prec=256
    lo,hi=arb('1310310.5534'),arb('1310311.7259')
    zeros=[acb.zeta_zero(n).imag for n in (2_345_665,2_345_666,2_345_667)]
    assert lo < zeros[0] < zeros[1] < zeros[2] < hi
    workers=min(cpu_count() or 1,8)
    blocks=[(NS*i//workers,NS*(i+1)//workers) for i in range(workers)]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results=[list(map(arb, row)) for block in executor.map(certify_block,blocks) for row in block]
    assert len(results)==NS
    gap=max(max((abs(x-y)/d).upper() for x,y,d in zip(a[:4],b[:4],DIAGONAL,strict=True))
            for a,b in zip(results,results[1:]))
    displacement=max(row[4].upper() for row in results)
    contraction=max(row[5].upper() for row in results)
    radius=max(row[6].upper() for row in results)
    floor_u=min(row[7].lower() for row in results)
    assert gap+radius < RADIUS and floor_u > arb('1.50866')
    print('实高度闭区间、覆盖格数、实际零点高度',lo.str(13),hi.str(13),NS,[z.str(16) for z in zeros])
    print('缩放残差、收缩率、根半宽、相邻中心距离的上界、逆耦合下界',
          *[n.str(12) for n in (displacement,contraction,radius,gap,floor_u)])


if __name__=='__main__':
    verify()
