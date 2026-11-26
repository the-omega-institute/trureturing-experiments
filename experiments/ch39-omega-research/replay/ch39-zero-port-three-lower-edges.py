"""前两个角与 CP 相位同时触界时，从 11415 进入到 11416 退出的端口界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat

prev = run_path(str(Path(__file__).with_name("ch39-zero-port-phase-branch.py")))
box, study, sharp = (prev[k] for k in ("box", "study", "sharp"))
S, T, h = sharp["S"], prev["T"], prev["half_tan"]
K = (1-S)*(1-T)
zero_ids = (*prev["ids"], 11415)


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def rows_relative(values):
    return [[v-base for v,base in zip(row,values[1])] for row in values]


def solve_centered(matrix,rhs):
    middle=arb_mat([[matrix[i,j].mid() for j in range(matrix.ncols())]
                    for i in range(matrix.nrows())])
    guess=middle.solve(arb_mat([[rhs[i,j].mid() for j in range(rhs.ncols())]
                                for i in range(rhs.nrows())]))
    return guess+matrix.solve(rhs-matrix*guess)


def candidate(coords,values,r):
    dif=rows_relative(values)
    # 用 11432 承担尚未分给另四个外部编号的权重，消去总质量方程。
    mid=r.mid()
    start=[a+mid*(b-a) for a,b in zip(values[0],values[1])]
    a,b,uc,dc,j=start
    domain=sharp["lo"].union(sharp["hi"])
    def solve_at(u):
        cu=K*(1-u-uc)+dc-S*T
        phase=[0,0,-h*cu,h*(1-u),1]
        mat=arb_mat([[coords[k]-coords[1] for k in range(2,6)],
                     [a+S*v for a,b,v,d,j in dif[2:]],
                     [b+T*v for a,b,v,d,j in dif[2:]],
                     [dot(phase,row) for row in dif[2:]]])
        drift=dif[0]
        rhs=arb_mat([[-mid*coords[1],-coords[1]],
                     [S*(1-uc)-a,drift[0]+S*drift[2]],
                     [T*(1-uc)-b,drift[1]+T*drift[2]],
                     [-j+h*(1-uc)*(K*uc-dc+S*T),dot(phase,drift)]])
        ws=solve_centered(mat,rhs)
        by_mass=[ws[i,1] for i in range(4)]
        exterior=[ws[i,0]+(r-mid)*by_mass[i] for i in range(4)]
        def reading(k,center):
            column=[row[k] for row in dif[2:]]
            base=center+dot([ws[i,0] for i in range(4)],column)
            return base+(r-mid)*(dot(by_mass,column)-drift[k])
        U,D=reading(2,uc),reading(3,dc)
        return mat,exterior,U,D
    mat,exterior,U,D=solve_at(domain)
    assert domain.contains_interior(U), (str(domain),str(U))
    # 求固定点映射 u ↦ U 的导数。只有第四行依赖 u。
    deriv=mat.solve(arb_mat([[0],[0],[0],[h*((D-dc)-K*(U-uc))]]))
    assert abs(dot([deriv[i,0] for i in range(4)],[row[2] for row in dif[2:]]))<1
    for _ in range(5):
        domain=domain.intersection(U)
        _,exterior,U,D=solve_at(domain)
    full=[1-r,r-sum(exterior),*exterior]
    assert all(w>0 for w in full[:4])
    assert dot(full,[row[4] for row in values])>0
    return full,(U,D)


def dual_plane(coords,costs,values,reference,u):
    uc,dc=reference
    phase=[0,0,-h*(K*(1-u-uc)+dc-S*T),h*(1-u),1]
    dif=rows_relative(values)
    mat=arb_mat([[coords[i]-coords[1],a+S*v,b+T*v,dot(phase,dif[i])]
                 for i,(a,b,v,d,j) in enumerate(dif) if i>=2])
    dual=solve_centered(mat,arb_mat([[c-costs[1]] for c in costs[2:]]))
    beta,first,second,cp=[dual[i,0] for i in range(4)]
    # 完整相位比代入两个角下边后的相位多出两个线性余量与一个非负乘积。
    assert cp>0
    assert first+cp*h*(T-(1-T)*u)>0
    assert second+cp*h*(S-(1-S)*u)>0
    coeff=[cp*z+first*x+second*y for z,x,y in zip(phase,(1,0,S,0,0),(0,1,T,0,0))]
    return mat.det(),beta,coeff


def verify():
    cost,strength,scale,xs,costs,inputs,_=prev["setup"]()
    height=acb.zeta_zero(11415).imag
    x=study["center"]/(height*height+study["offset"])-1
    xs.append(x)
    costs.append(cost(x)[0])
    inputs.append(study["native"]["plane_response"](height))
    theta=arb(0).union(1)
    coords=[x+theta*scale*y for x,y in zip(xs,costs)]
    delta=arb(2)**-28
    cube=[c+arb(0,delta) for c in box["centers"]]
    data=[box["tight_data"](t,box["centers"],[delta]*4) for t in inputs]
    rows,ds=zip(*data)
    determinant_sign=None
    planes=None
    signs=(1,-1,-1,-1)
    for k in range(100):
        # 覆盖 27.7% 到 28.2%，允许与转折相邻的辅助解带符号。
        outside=arb(5540+k).union(5541+k)/20000
        weights,ref=candidate(coords,rows,outside)
        uc,dc=ref
        phase=[0,0,-h*(K*(1-2*uc)+dc-S*T),h*(1-uc),1]
        dif=rows_relative(rows)
        def features(row):
            a,b,u,d,j=row
            return [a+S*u,b+T*u,dot(phase,row)]
        mat=arb_mat([[coords[i]-coords[1],*features(dif[i])] for i in range(2,6)])
        speeds=solve_centered(mat.transpose(),arb_mat([[v] for v in [-coords[1],*features(dif[0])]]))
        assert speeds[2,0]<0<speeds[3,0]
        for tangent,u in ((True,uc),(False,sharp["lo"]),(False,sharp["hi"])):
            det,beta,coeff=dual_plane(coords,costs,rows,ref,u)
            assert 0 not in det
            if determinant_sign is None:
                determinant_sign=det>0
            assert (det>0)==determinant_sign
            assert 1-theta*scale*beta>0 and 1+(1-theta)*scale*beta>0
            if tangent:
                for j,s in enumerate(signs):
                    deriv=-sum(w*dot(coeff,der[j]) for w,der in zip(weights,ds))
                    assert -s*deriv>0,(outside.str(9),j,deriv.str(9))
            point=[beta,*coeff]
            planes=point if planes is None else [a.union(b) for a,b in zip(planes,point)]
    corner=[c+sign*delta for c,sign in zip(box["centers"],signs)]
    best=[box["flavor_data"](t,corner)[0] for t in inputs]
    for t in (theta,arb(0),arb(1)/2,arb(1)):
        v=[x+t*scale*y for x,y in zip(xs,costs)]
        def distribution(r):
            return candidate(v,best,r)[0]
        entry=prev["bracket"](lambda r: -distribution(r)[5],".277",".278")
        exit_point=prev["bracket"](lambda r: distribution(r)[4],".280",".282")
        for name,r in (("entry",entry),("28 percent",arb(7)/25),("exit",exit_point)):
            ws=distribution(r)
            optimum=dot(ws,costs)
            print("Direction",t.str(3),name,"outside",r.str(19),
                  "gauge tolerance",(t*scale*optimum).str(19),
                  "port tolerance",((1-t)*strength*optimum).str(19),
                  "weights",*[w.str(12) for w in ws],flush=True)
    beta,*coeff=planes
    def slack(x,c,obs):
        diff=[a-b for a,b in zip(obs,rows[1])]
        return c-costs[1]-beta*(x+theta*scale*c-coords[1])-dot(coeff,diff)
    assert slack(0,0,rows[0])<0
    for lower in (S,T):
        assert -sharp["hi"]+(1+sharp["lo"])*lower>0
    whole=arb(0).union(1)
    outer=[whole,arb(0).union(arb(104)/867),whole,whole,arb(-1).union(1)]
    for n in (8499,16001):
        height=acb.zeta_zero(n).imag
        x=study["center"]/(height**2+study["offset"])-1
        c,dx=cost(x)
        speed=dx/(1+theta*scale*dx)
        assert (speed>beta if n==8499 else speed<beta)
        assert slack(x,c,outer)>0
    floor=arb(1)
    for start in range(8500,16001,256):
        count=min(256,16001-start)
        zeros=acb.zeta_zeros(start,count)
        assert len(zeros)==count
        for n,zero in enumerate(zeros,start):
            if n in zero_ids:
                continue
            x=study["center"]/(zero.imag**2+study["offset"])-1
            c=cost(x)[0]
            t=study["native"]["plane_response"](zero.imag)
            obs=box["observations"](box["flavor"]["mixing"](t,cube))
            gap=slack(x,c,obs)
            if not gap>0:
                gap=slack(x,c,box["tight_data"](t,box["centers"],[delta]*4)[0])
            assert gap>0,(n,gap.str(12))
            floor=min(floor,gap.lower())
    print("Strict finite slack",floor.str(13))


if __name__ == "__main__":
    verify()
