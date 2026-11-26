"""两个角和相位同时触界时，外部权重增至三分之一附近的支撑与第五项味检验。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat

edges = run_path(str(Path(__file__).with_name("ch39-zero-port-three-lower-edges.py")))
prev, box, study, sharp = (edges[key] for key in ("prev", "box", "study", "sharp"))
S, T, K, h = (edges[key] for key in ("S", "T", "K", "h"))
dot = edges["dot"]
# 每行给出退出与进入的编号、两个端点的包围区间（万分数）及最优味系数角点。
branches = (
    (11406, 11405, (2808,2810), (2848,2851), (1,-1,-1,-1)),
    (11405, 11389, (2848,2851), (3119,3122), (1,-1,-1,-1)),
    (11389, 11388, (3119,3122), (3127,3130), (1,-1,-1,-1)),
    (11388, 11372, (3127,3130), (3231,3234), (1,-1,-1,-1)),
    (11372, 11304, (3231,3234), (3433,3436), (1,-1,1,-1)),
)


def moment_radius(u,d):
    return (K*u-d+S*T)**2/(4*arb("1.141").cos()**2*S*T*K*u)


def certify_branch(spec, spectral, scale, strength):
    leaving,entering,left,right,signs = spec
    ids = (11431,11432,11415,11430,leaving,entering)
    xs,costs,inputs,rows,derivatives = zip(*(spectral[n] for n in ids))
    theta = arb(0).union(1)
    coords = [x+theta*scale*c for x,c in zip(xs,costs)]
    dif = edges["rows_relative"](rows)
    corner = [c+s*arb(2)**-28 for c,s in zip(box["centers"],signs)]
    best = [box["flavor_data"](t,corner)[0] for t in inputs]
    # 跨角点处只有第三个味系数在盒内移动。
    sliding = corner.copy()
    sliding[2] = box["centers"][2]+arb(0,arb(2)**-28)
    slide_rows = None
    if leaving==11388:
        sliding[2] = box["centers"][2]
        slide_rows = [box["tight_data"](t,sliding,[0,0,arb(2)**-28,0])[0] for t in inputs]
    combined = None
    sign_det = None
    radius = None
    for i in range(2*left[0],2*right[1]):
        r = arb(i).union(i+1)/20000
        weights,ref = edges["candidate"](coords,rows,r)
        u,d = ref
        phase = [0,0,-h*(K*(1-2*u)+d-S*T),h*(1-u),1]
        def features(row):
            a,b,u,d,j = row
            return [a+S*u,b+T*u,dot(phase,row)]
        matrix = arb_mat([[coords[i]-coords[1],*features(dif[i])] for i in range(2,6)])
        speed = edges["solve_centered"](matrix.transpose(),
                                        arb_mat([[v] for v in [-coords[1],*features(dif[0])]]))
        assert speed[2,0]<0<speed[3,0]
        midpoint = r.mid()
        target_rows = slide_rows if leaving==11388 and i>=2*(right[0]-2) else best
        _,(u_mid,d_mid) = edges["candidate"](coords,target_rows,midpoint)
        # 中点求值再沿 r 作区间积分，保留消去权重时相抵的量。
        velocity = [speed[k,0] for k in range(4)]
        du,dd = [rows[1][k]-rows[0][k]+dot(velocity,[row[k] for row in dif[2:]]) for k in (2,3)]
        u_enclosure = u_mid+(r-midpoint)*du
        numerator = K*u_mid-d_mid+S*T+(r-midpoint)*(K*du-dd)
        radius_squared = numerator**2/(4*arb("1.141").cos()**2*S*T*K*u_enclosure)
        radius = radius_squared if radius is None else radius.union(radius_squared)
        if leaving==11388 and i>=2*(right[0]-2):
            # 在相邻角点间，第四个独立方程仍取真实相位切线。
            mean = [sum(w*ds[2][k] for w,ds in zip(weights,derivatives)) for k in range(5)]
            shift = edges["solve_centered"](matrix.transpose(),arb_mat([[v] for v in [0,*[-z for z in features(mean)]]]))
            assert shift[2,0]>0
        for tangent,v in ((True,u),(False,sharp["lo"]),(False,sharp["hi"])):
            det,beta,coeff = edges["dual_plane"](coords,costs,rows,ref,v)
            assert 0 not in det
            if sign_det is None:
                sign_det = det>0
            assert sign_det==(det>0)
            assert 1-theta*scale*beta>0 and 1+(1-theta)*scale*beta>0
            if tangent:
                for k,sign in enumerate(signs):
                    change = sum(w*dot(coeff,ds[k]) for w,ds in zip(weights,derivatives))
                    assert sign*change>0,(ids,i,k,change.str(8))
            plane = [beta,*coeff]
            combined = plane if combined is None else [a.union(b) for a,b in zip(combined,plane)]
    print("Pair",leaving,entering,"radius squared",radius.str(10),flush=True)
    assert radius<arb(".9999"),radius.upper().str(9)
    low,high = (tuple(f"{n//10000}.{n%10000:04d}" for n in side)
                for side in (left,right))
    for t in (theta,arb(0),arb(1)):
        v = [x+t*scale*c for x,c in zip(xs,costs)]
        def solution(r):
            return edges["candidate"](v,best,r)
        start = prev["bracket"](lambda r: -solution(r)[0][5],*low)
        end = prev["bracket"](lambda r: solution(r)[0][4],*high)
        assert start<end
        samples = [("entry",start),("exit",end)]
        if start<arb(1)/3<end:
            samples.append(("one third",arb(1)/3))
        for label,r in samples:
            ws,(u,d) = solution(r)
            optimum = dot(ws,costs)
            print("Direction",t.str(2),label,r.str(18),"port",((1-t)*strength*optimum).str(17),
                  "gauge",(t*scale*optimum).str(17),"radius squared",moment_radius(u,d).str(14),flush=True)
    beta,*coeff = combined
    def reading(obs):
        return dot(coeff,[a-b for a,b in zip(obs,rows[1])])
    def radial_cost(x,c):
        return c-costs[1]-beta*(x+theta*scale*c-coords[1])
    assert radial_cost(0,0)<reading(rows[0])
    upper = sharp["compare_torus"](coeff[4],reading)
    def slack(x,c,obs):
        return radial_cost(x,c)-reading(obs)
    return ids,beta,radial_cost,upper,slack


def single_unitary_bound(spectral,strength):
    ids = (11431,11432,11415,11430,11372,11304)
    xs,costs,inputs,rows,ds = zip(*(spectral[n] for n in ids))
    lo,hi = sharp["lo"],sharp["hi"]
    factor = arb("1.141").sin()*(S*T*K).sqrt()
    def lower(u):
        return factor*(1-u)*u.sqrt()
    slope = (lower(hi)-lower(lo))/(hi-lo)
    def cut(row):
        return row[4]-lower(lo)-slope*(row[2]-lo)
    # 在原代价减去 m 倍第五项必要条件的代价上，重新解支撑平面。
    multiplier = arb(1)/100
    altered = [c-multiplier*cut(o) for c,o in zip(costs,rows)]
    r = arb(1)/3
    ws,ref = edges["candidate"](xs,rows,r)
    deficit = -(cut(rows[0])+r*(cut(rows[1])-cut(rows[0]))+
                sum(w*(cut(o)-cut(rows[1])) for w,o in zip(ws[2:],rows[2:])))
    assert deficit>0,deficit.str(15)
    planes = []
    signs = []
    for u in (lo,hi):
        det,beta,coeff = edges["dual_plane"](xs,altered,rows,ref,u)
        assert 0 not in det
        signs.append(det>0)
        coeff = [v+multiplier*k for v,k in zip(coeff,(0,0,-slope,0,1))]
        planes.append([beta,*coeff])
    assert signs[0]==signs[1]
    beta,*coeff = [a.union(b) for a,b in zip(*planes)]
    def reading(obs):
        return dot(coeff,[a-b for a,b in zip(obs,rows[1])])
    def radial_cost(x,c):
        return c-costs[1]-beta*(x-xs[1])
    def slack(x,c,obs):
        return radial_cost(x,c)-reading(obs)
    assert slack(0,0,rows[0])<0
    extra=strength*multiplier*deficit
    print("Uniform unitary port surcharge at one third exceeds",extra.lower().str(17),flush=True)
    upper=sharp["compare_torus"](coeff[4],reading)
    return ids,beta,radial_cost,upper,slack


def verify():
    cost,strength,scale,*_ = prev["setup"]()
    delta = arb(2)**-28
    cube = [c+arb(0,delta) for c in box["centers"]]
    ids = {11431,11432,11415,11430}
    for leaving,entering,*_ in branches:
        ids.update((leaving,entering))
    spectral = {}
    for n in sorted(ids):
        height = acb.zeta_zero(n).imag
        x = arb(0) if n==11431 else study["center"]/(height**2+study["offset"])-1
        c = arb(0) if n==11431 else cost(x)[0]
        t = study["native"]["plane_response"](height)
        obs,ds = box["tight_data"](t,box["centers"],[delta]*4)
        spectral[n] = (x,c,t,obs,ds)
    checks = [certify_branch(spec,spectral,scale,strength) for spec in branches]
    checks.append(single_unitary_bound(spectral,strength))
    theta = arb(0).union(1)
    for n in (8499,16001):
        height = acb.zeta_zero(n).imag
        x = study["center"]/(height**2+study["offset"])-1
        c,dx = cost(x)
        slope = dx/(1+theta*scale*dx)
        for support,beta,radial_cost,upper,_ in checks:
            assert (slope>beta if n==8499 else slope<beta)
            assert radial_cost(x,c)>upper,(n,support,upper.str(8))
    floor = arb(1)
    for start in range(8500,16001,256):
        count = min(256,16001-start)
        zeros = acb.zeta_zeros(start,count)
        assert len(zeros)==count
        for n,zero in enumerate(zeros,start):
            if n==11431:
                continue
            x = study["center"]/(zero.imag**2+study["offset"])-1
            c = cost(x)[0]
            t = study["native"]["plane_response"](zero.imag)
            obs = box["observations"](box["flavor"]["mixing"](t,cube))
            for support,_,_,_,slack in checks:
                if n in support:
                    continue
                gap=slack(x,c,obs)
                if not gap>0:
                    obs=box["tight_data"](t,box["centers"],[delta]*4)[0]
                    gap=slack(x,c,obs)
                assert gap>0,(support,n,gap.str(9))
                floor=min(floor,gap.lower())
    print("Full spectrum strict finite slack",floor.str(13))


if __name__=="__main__":
    verify()
