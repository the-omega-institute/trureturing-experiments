"""三个实际壳上共同复味几何与相同中心物理矩阵的必要铅笔区间。"""
from pathlib import Path
from runpy import run_path
from flint import arb, acb_mat, ctx


def setup(root):
    lib=run_path(str(root/"ch39-neighbour-all-frame-support.py"))
    geom=run_path(str(root/"ch39-complex-four-site-unit-quark-amplitudes.py"))
    ck=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    ctx.prec=256
    V=acb_mat(ck["standard_matrix"](*[arb(i).sin() for i in (".22704",".04275",".00377")],arb("1.139")))
    mass=[list(map(arb,row)) for row in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    up,down=[geom["diag"]([m/row[2] for m in row]) for row in mass]
    A,B=geom["adj"](V)*up*V,down
    assert not (A*B-B*A).det().contains(0)
    floors=run_path(str(root/"ch39-triple-native-bloch-radius.py"))["TRIPLE_RATIO_FLOORS"]
    return lib,geom,V,mass,A,B,tuple(map(arb,floors))


def physical_blocks(y,basis,pencil,V,masses,geom,lib,ratios):
    adj=geom["adj"]
    remote=1-basis
    yrep=[V*y,y]
    stats=[lib["moments"](h,spec) for h,spec in zip(yrep,masses,strict=True)]
    _,prob,var=stats[remote]
    if not var>0:
        return None
    norm=var.sqrt()
    h=yrep[remote].entries()
    spec=masses[remote]
    shifts=[sum(prob[j]*(spec[i]-spec[j]) for j in range(3) if j!=i) for i in range(3)]
    z=acb_mat([[shift*hj/norm] for shift,hj in zip(shifts,h,strict=True)])
    x=acb_mat([[((spec[(i+2)%3]-spec[(i+1)%3])*h[(i+1)%3]*h[(i+2)%3]).conjugate()/norm] for i in range(3)])
    xs=[x,adj(V)*x] if remote==0 else [V*x,x]
    zs=[z,adj(V)*z] if remote==0 else [V*z,z]
    coeff=pencil*masses[basis][2]/masses[remote][2]
    w=[None,None]
    w[remote]=norm
    w[basis]=-coeff*norm
    result=[]
    for j,(xq,zq,mass) in enumerate(zip(xs,zs,masses,strict=True)):
        p=sum(abs(xq[i,0])**2*mass[i] for i in range(3)).intersection(mass[0].union(mass[2]))
        zz=sum(abs(zq[i,0])**2*mass[i] for i in range(3)).intersection(mass[0].union(mass[2]))
        k=abs(sum(xq[i,0].conjugate()*zq[i,0]*mass[i] for i in range(3)))
        s=stats[j][0].intersection(mass[0].union(mass[2]))
        bounds=mass[0].union(mass[2])
        p=p.intersection((sum(mass)-s-zz).intersection(bounds))
        zz=zz.intersection((sum(mass)-s-p).intersection(bounds))
        if not p>0 or not s>0: return None
        result.append((p,s,zz,k,w[j],mass[0]*mass[1]*mass[2]))
    return result,coeff,ratios


def possible(box,lo,hi,basis):
    """Q 为几何参数 a²。只用包围区间的端点舍去不可能的参数。"""
    blocks,coeff,ratios=box
    radius=(1-lo)
    Q=lo.union(hi)
    if not radius>0: return False
    x_window=[]
    u_low=v_low=arb(0)
    for (p,s,zz,k,w,det),eps in zip(blocks,ratios,strict=True):
        c=eps/(1+eps)**2
        k2=max(arb(0),k.lower())**2/radius
        disc=(1-4*c*(1+k2/p**2)).upper()
        if disc<0: return False
        upper=(p*(1-2*c+disc.sqrt())/(2*c)).upper()
        lower=((c*p**2+k2)/(c*upper)).lower()
        if hi<1:
            upper=min(upper,((det/(s*c*(1-hi))).sqrt()-p).upper())
        if lower>upper:return False
        # F=(u-v)w。给定几何区间时用两条仿射式分别包围两端。
        low=(zz-upper+Q*(upper-s)).lower()
        high=(zz-lower+Q*(lower-s)).upper()
        x_window.append(low.union(high))
        u_low=max(u_low,(w/s).lower())
        v_low=max(v_low,(-w/s).lower())
    if u_low*v_low>hi: return False
    other=1-basis
    wr=blocks[other][4]
    lower=x_window[other].lower()
    upper=x_window[other].upper()
    if u_low>0: lower=max(lower,((u_low-hi/u_low)*wr).lower())
    if v_low>0: upper=min(upper,((hi/v_low-v_low)*wr).upper())
    if lower>upper:
        return False
    return (-coeff*lower.union(upper)).overlaps(x_window[basis])


def covered(box,basis,depth_limit):
    """同一铅笔区间内，覆盖全部 Q=a²；覆盖失败就继续缩窄铅笔参数。"""
    sub=[(arb(0),arb(1),0)]
    tested=0
    while sub:
        lo,hi,depth=sub.pop()
        tested+=1
        if not possible(box,lo,hi,basis):
            continue
        if depth==depth_limit:
            return 0
        mid=(lo+hi)/2
        sub.extend(((lo,mid,depth+1),(mid,hi,depth+1)))
    return tested


def verify(root):
    lib,geom,V,mass,A,B,ratios=setup(root)
    for basis,(a,b) in enumerate(((A,B),(B,A))):
        for branch in range(3):
            stack=[(arb(-1),arb(1),0)]
            count=depth=geometry_count=0
            while stack:
                lo,hi,level=stack.pop()
                count+=1
                depth=max(depth,level)
                lines=lib["lines"](a,b,lo,hi,(branch,),geom)
                box=physical_blocks(lines[branch],basis,lo.union(hi),V,mass,geom,lib,ratios) if branch in lines else None
                tested=covered(box,basis,min(12,level+2)) if box else 0
                if tested:
                    geometry_count+=tested
                    continue
                assert level<23,(basis,branch,lo,hi)
                mid=(lo+hi)/2
                stack.extend(((lo,mid,level+1),(mid,hi,level+1)))
            print("坐标图、本征支、铅笔区间、成功覆盖的几何区间、最大铅笔深度",
                  basis,branch,count,geometry_count,depth,flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
