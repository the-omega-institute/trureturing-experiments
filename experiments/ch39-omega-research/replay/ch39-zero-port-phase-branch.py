"""CP 相位触及下边后，用四个非中心编号延拓取等分布。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

box = run_path(str(Path(__file__).with_name("ch39-zero-port-flavor-box.py")))
study, sharp = box["study"], box["sharp"]
ids = (11431, 11432, 11430, 11406, 11416)
T = arb(".041775")**2
half_tan = sharp["tau"]/2


def dot(left, right):
    return sum(a*b for a,b in zip(left,right))


def branch(coords, rows):
    # 以外部总权重 r 和 11416 的绝对权重 z 为参数。
    margins = [row[1]-T*(1-row[2]) for row in rows]
    _, right, left, old, new = margins
    v0, v1, v2, v3, v4 = coords
    matrix = arb_mat([[v2-v1, v3-v1], [left-right, old-right]])
    rhs = arb_mat([[-v0, v0-v1, v1-v4],
                   [-margins[0], margins[0]-right, right-new]])
    sol = matrix.solve(rhs)
    def weights(r,z):
        left_w = sol[0,0]+r*sol[0,1]+z*sol[0,2]
        old_w = sol[1,0]+r*sol[1,1]+z*sol[1,2]
        return [1-r,r-left_w-old_w-z,left_w,old_w,z]
    def coefficients(col):
        c0,c1,c2,c3,c4 = col
        return [c0+sol[0,0]*(c2-c1)+sol[1,0]*(c3-c1),
                c1-c0+sol[0,1]*(c2-c1)+sol[1,1]*(c3-c1),
                c4-c1+sol[0,2]*(c2-c1)+sol[1,2]*(c3-c1)]
    mean = [coefficients(col) for col in zip(*rows)]
    def evaluate(r,z):
        return [c+r*d+z*e for c,d,e in mean]
    def phase(r,z):
        a,b,u,v,j = evaluate(r,z)
        real2 = (1-T)*u*(1-u)-(1-u)*v+(T-(1-T)*u)*a
        return j-half_tan*real2
    def gradient(r,z):
        a,b,u,v,j = evaluate(r,z)
        cu = (1-T)*(1-2*u)+v-(1-T)*a
        ca = T-(1-T)*u
        normal = [-half_tan*ca,0,-half_tan*cu,half_tan*(1-u),1]
        return [dot(normal,col) for col in list(zip(*mean))[1:]]
    return weights,evaluate,phase,gradient,mean,coefficients


def setup():
    ctx.prec = 240
    cost, strength, tangent = box["joint"]["port_cost"](include_tangent=True)
    scale = -strength/tangent
    zeros = [acb.zeta_zero(i).imag for i in ids]
    xs = [study["center"]/(h*h+study["offset"])-1 for h in zeros]
    xs[0]=arb(0)
    costs = [arb(0)]+[cost(x)[0] for x in xs[1:]]
    inputs = [study["native"]["plane_response"](h) for h in zeros]
    delta=arb(2)**-28
    corner=[v+s*delta for v,s in zip(box["centers"],(1,-1,-1,-1))]
    rows = [box["flavor_data"](t,corner)[0] for t in inputs]
    return cost,strength,scale,xs,costs,inputs,rows


def root_z(model,r):
    _,_,phase,gradient,_,_ = model
    lower,upper = arb("-.002"),arb(".01")
    full = lower.union(upper)
    assert phase(r,lower)<0<phase(r,upper)
    assert gradient(r,full)[1]>0
    for _ in range(6):
        mid = full.mid()
        full = full.intersection(mid-phase(r,mid)/gradient(r,full)[1])
    return full


def bounds(model, r, costs):
    weights,evaluate,_,gradient,mean,coefficients = model
    z = root_z(model,r)
    dr,dz = gradient(r,z)
    speed = -dr/dz
    assert speed>0
    da,_,du,dv,_ = [d+speed*e for _,d,e in mean]
    curvature = -2*half_tan*du*((1-T)*(du+da)-dv)/dz
    assert curvature>0 and coefficients(costs)[2]>0
    lower = arb(".224955")**2
    first = [a+lower*u for a,u in zip(mean[0],mean[2])]
    assert first[1]+speed*first[2]<0
    a,b,u,v,j = evaluate(r,z)
    assert sharp["lo"]<u<sharp["hi"] and j>0
    assert -sharp["hi"]+(1+sharp["lo"])*arb(".224955")**2>0
    assert a<arb(".225065")**2*(1-u)
    ws = weights(r,z)
    assert all(w>0 for w in ws[:4])
    return ws,(a,b,u,v,j)


def bracket(func,lo,hi):
    left,right = arb(lo),arb(hi)
    assert func(left)>0>func(right)
    for _ in range(76):
        mid=(left+right)/2
        val=func(mid)
        if 0 in val:
            break
        if val>0:
            left=mid
        else:
            right=mid
    return left.union(right)


def endpoints(coords,rows,costs,strength,theta,scale):
    model = branch(coords,rows)
    weights,evaluate,_,_,_,_ = model
    start = bracket(lambda r: -root_z(model,r),".212",".213")
    def first_margin(r):
        a,_,u,_,_ = evaluate(r,root_z(model,r))
        return a-arb(".224955")**2*(1-u)
    end = bracket(first_margin,".276",".280")
    for name,r in (("entry",start),("exit",end)):
        z = root_z(model,r)
        ws = weights(r,z)
        print(name,"outside weight",r.str(19),"port tolerance", ((1-theta)*dot(ws,costs)*strength).str(19),
              "gauge tolerance", (theta*scale*dot(ws,costs)).str(19),
              "11416 weight",z.str(14),flush=True)
    return start,end


def phase_normal(u,reference):
    a,b,uc,v,j = reference
    cu = (1-T)*(1-u-uc)+v-(1-T)*a
    return [-half_tan*(T-(1-T)*u),0,-half_tan*cu,half_tan*(1-u),1]


def support_plane(u,reference,coords,costs,rows):
    normal = phase_normal(u,reference)
    def differences(row):
        delta = [a-b for a,b in zip(row,rows[1])]
        return [delta[1]+T*delta[2],dot(normal,delta)]
    matrix = arb_mat([[coords[i]-coords[1],*differences(rows[i])] for i in (2,3,4)])
    rhs = arb_mat([[costs[i]-costs[1]] for i in (2,3,4)])
    middle = arb_mat([[matrix[i,j].mid() for j in range(3)] for i in range(3)])
    dual = middle.solve(rhs)
    dual += matrix.solve(rhs-matrix*dual)
    beta,lam,mu = [dual[i,0] for i in range(3)]
    assert lam>0 and mu>0
    # 保存系数，以便一次检查整条分支及两个相位割平面。
    coeff = [mu*a+lam*b for a,b in zip(normal,(0,1,T,0,0))]
    return matrix.det(), beta, coeff


def merged_planes(coords,costs,inputs,theta,scale):
    radius=arb(2)**-28
    data=[box["tight_data"](t,box["centers"],[radius]*4) for t in inputs]
    rows,ds=zip(*data)
    model=branch(coords,rows)
    # 取等方程对物理范围的 z 至多有一根。
    extent = arb(".212").union(arb(".280"))
    assert model[3](extent,arb(0).union(1))[1]>0
    det_sign=None
    union=None
    for i in range(136):
        r=(arb(424+i).union(425+i))/2000
        ws,ref=bounds(model,r,costs)
        for tangent,u in ((False,sharp["lo"]),(False,sharp["hi"]),(True,ref[2])):
            det,beta,coeff=support_plane(u,ref,coords,costs,rows)
            assert 0 not in det
            if det_sign is None:
                det_sign=det>0
            assert (det>0)==det_sign
            assert 1-theta*scale*beta>0 and 1+(1-theta)*scale*beta>0
            if tangent:
                for j,sgn in enumerate((-1,1,1,1)):
                    deriv=-sum(w*dot(coeff,der[j]) for w,der in zip(ws,ds))
                    assert deriv*sgn>0,(r.str(7),j,deriv.str(8))
            vals=[beta,*coeff]
            union=vals if union is None else [a.union(b) for a,b in zip(union,vals)]
    beta,*coeff=union
    def slack(x,c,obs):
        delta=[a-b for a,b in zip(obs,rows[1])]
        return c-costs[1]-beta*(x+theta*scale*c-coords[1])-dot(coeff,delta)
    assert slack(0,0,rows[0])<0
    return beta,slack


def verify():
    cost,strength,scale,xs,costs,inputs,rows=setup()
    theta=arb(0).union(1)
    coords=[x+theta*scale*c for x,c in zip(xs,costs)]
    for t in (theta,arb(0),arb(1)/2,arb(1)):
        print("Budget direction",t.str(3),flush=True)
        endpoints([x+t*scale*c for x,c in zip(xs,costs)],rows,costs,strength,t,scale)
    beta,slack=merged_planes(coords,costs,inputs,theta,scale)
    radius=arb(2)**-28
    coeff_box=[c+arb(0,radius) for c in box["centers"]]
    whole=arb(0).union(1)
    outer_obs=[whole,arb(0).union(arb(104)/867),whole,whole,arb(-1).union(1)]
    for n in (8499,16001):
        h=acb.zeta_zero(n).imag
        x=study["center"]/(h*h+study["offset"])-1
        c,dx=cost(x)
        speed=dx/(1+theta*scale*dx)
        assert (speed>beta if n==8499 else speed<beta)
        assert slack(x,c,outer_obs)>0
    floor=arb(1)
    for start in range(8500,16001,256):
        count=min(256,16001-start)
        zeros=acb.zeta_zeros(start,count)
        assert len(zeros)==count
        for n,zero in enumerate(zeros,start):
            if n in ids:
                continue
            h=zero.imag
            x=study["center"]/(h*h+study["offset"])-1
            c=cost(x)[0]
            t=study["native"]["plane_response"](h)
            obs=box["observations"](box["flavor"]["mixing"](t,coeff_box))
            gap=slack(x,c,obs)
            if not gap>0:
                gap=slack(x,c,box["tight_data"](t,box["centers"],[radius]*4)[0])
            assert gap>0,(n,gap.str(10))
            floor=min(floor,gap.lower())
    print("Strict finite support slack",floor.str(12))


if __name__ == "__main__":
    verify()
