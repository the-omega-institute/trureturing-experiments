"""左邻实际壳上与平凡平面任意酉转角无关的字符尺度禁区。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def interval(lo, hi):
    return arb(lo).union(arb(hi))


def dot(x, y):
    return sum(a*b for a, b in zip(x, y, strict=True))


def cross(x, y):
    return [x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3] for i in range(3)]


def certify_invariants(root, scales=None, collect=None):
    pencil = run_path(str(root/"ch39-triple-geometry-pencil.py"))
    lib, geom, V, mass, A, B, _ = pencil["setup"](root)
    u, v = interval(".003839",".003841"), interval(".03839",".03841")
    boxes = [interval(".003109",".003123"),interval(".0002808",".0002830"),
             interval(".9959",".9962"),interval(".8122",".8132"),interval(".3959",".3969")]
    if scales is not None:
        u, v = scales
    a2 = u*v
    b = (1-a2).sqrt()
    count = hits = deepest = 0
    for chart, (base, other) in enumerate(((A,B),(B,A))):
        stack = [(arb(-1),arb(1),(0,1,2),0)]
        while stack:
            lo, hi, ids, depth = stack.pop()
            count += 1
            deepest = max(deepest,depth)
            t = lo.union(hi)
            lines = lib["lines"](base,other,lo,hi,ids,geom)
            keep = []
            for idx in ids:
                box = pencil["physical_blocks"](lines[idx],chart,t,V,mass,geom,lib,(),complex_cross=True) if idx in lines else None
                if box is None:
                    keep.append(idx)
                    continue
                for sign in (-1,1):
                    chars = [(s-sign*w/u,s+sign*w/v) for d,s,z,k,w,det in box[0]]
                    if any(c<0 or c>1 for pair in chars for c in pair):
                        continue
                    if collect is not None and any(det-a2*d*h[0]*h[1]<0 for (d,s,z,k,w,det),h in zip(box[0],chars,strict=True)):
                        continue
                    bloch = []
                    for d,s,z,k,w,det in box[0]:
                        low = (z-a2*s-(u-v)*sign*w)/(1-a2)
                        trace = d+low
                        if not trace>0:
                            break
                        bloch.append((trace,[2*k.real/(b*trace),-2*k.imag/(b*trace),(d-low)/trace]))
                    if len(bloch)==2:
                        (tu,x),(td,y) = bloch
                        radius = dot(x,x).sqrt()
                        if radius>0:
                            para = dot(x,y)/radius
                            perp2 = dot(y,y)-para**2
                            if perp2>0:
                                values = (tu,td,radius,para,perp2.sqrt())
                                ok = collect(values,x[2]/radius) if collect is not None else all(bound.contains(val) for bound,val in zip(boxes,values,strict=True))
                                if ok:
                                    hits += 1
                                    continue
                    keep.append(idx)
                    break
            if keep:
                assert depth<34,(chart,lo,hi,keep)
                mid=(lo+hi)/2
                stack.extend(((lo,mid,keep,depth+1),(mid,hi,keep,depth+1)))
    assert hits>0
    print("铅笔盒、落入不变量界的本征支取向盒、最大二分深度",count,hits,deepest)
    return boxes


def native_coordinates(row):
    return [row[0,0].real+row[1,1].real,2*row[0,1].real,-2*row[0,1].imag,
            row[0,0].real-row[1,1].real]


def forms(rows):
    columns=[native_coordinates(row) for row in rows]
    trivial=[v for v in columns if v[0]>0]
    assert len(trivial)==7
    directions=set()
    for tri in combinations(trivial,3):
        minors=[arb_mat([[x for j,x in enumerate(row) if j!=i] for row in tri]).det() for i in range(4)]
        n=[(-1)**j*float(x.mid()) for j,x in enumerate(minors)]
        norm=max(map(abs,n))
        if norm:
            co=tuple(round(10**8*x/norm) for x in n)
            directions.update((co,tuple(-x for x in co)))
    result=[]
    for co in sorted(directions):
        length=arb(sum(x*x for x in co)).sqrt()
        c=[x/length for x in co]
        support=sum(max(arb(0),dot(c,x).upper()) for x in columns)
        result.append((c[0],c[1:],support))
    return result


def native_rays(up):
    cs=[native_coordinates(row) for row in up]
    rays=[(j+1,c[0],[x/c[0] for x in c[1:]]) for j,c in enumerate(cs) if c[0]>0]
    assert len(rays)==7
    return rays


def exclude_rays(boxes,rays,support,allowed=()):
    tu, td, radius, para, perp = boxes
    defect=1-radius**2
    results=[]
    down_norm2=para**2+perp**2
    for j,trace,ray in rays:
        if j in allowed:
            continue
        distance=min((1-dot(ray,r)).lower() for k,_,r in rays if k!=j)
        L=2/distance
        if distance*(1-trace/tu)>defect:
            results.append((j,0,None))
            continue
        assert L*defect<1
        root_disc=(1-L*defect).sqrt()
        cos_loss=((L-1)*defect**2/(radius*(radius+root_disc)*(1+root_disc)*(1+radius))).upper()
        sin_bound=(2*cos_loss).sqrt()
        first=cross(ray,[0,0,1])
        norm=dot(first,first).sqrt()
        first=[x/norm for x in first]
        second=cross(ray,first)
        normals=[]
        for a,n,sup in support:
            if a+dot(n,n).sqrt()*down_norm2.sqrt()>sup/td:
                normals.append((a-sup/td,para*dot(n,ray),perp*dot(n,first),perp*dot(n,second),dot(n,n)*down_norm2))
        stack=[(0,0)]
        visited=0
        worst=arb(100)
        while stack:
            depth,l=stack.pop()
            visited+=1
            sn,co=(arb(2*l+1,1)*arb.pi()/(2**depth)).sin_cos()
            margins=[]
            for base,k0,ka,kb,nr in normals:
                k=k0+ka*co+kb*sn
                tang=max(arb(0),(nr-k*k).upper()).sqrt()
                margins.append((base+k-abs(k)*cos_loss-tang*sin_bound).lower())
            margin=max(margins)
            if margin>0:
                worst=min(worst,margin)
            else:
                assert depth<19,(j,depth,l,margin)
                stack.extend(((depth+1,2*l),(depth+1,2*l+1)))
        results.append((j,visited,worst))
    return results


def verify(root):
    boxes=certify_invariants(root)
    fixed=run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    up,down=fixed["native"](root,2_306_888_100_084)
    ctx.prec=256
    for j,visited,worst in exclude_rays(boxes,native_rays(up),forms(down)):
        if visited:
            print("上型最近谱层、下型方向角覆盖盒、最小严格违约量",j,visited,worst.str(12))
        else:
            print("上型最近谱层由容量排除",j)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
