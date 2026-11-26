"""四个味作用系数改变时，四点取等分布及其端口代价的变化。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, arb_mat, ctx

root = Path(__file__).parent
joint = run_path(str(root/"ch39-zero-port-ckm.py"))
sharp = joint["sharp"]
study = sharp["study"]
flavor = study["flavor"]
centers = [arb(n)/2048 for n in flavor["codes"]]


def observations(matrix):
    a, b, u, v = [abs(matrix[i,j])**2 for i,j in ((0,1),(1,2),(0,2),(2,0))]
    j = (matrix[0,1]*matrix[1,2]*matrix[0,2].conjugate()*matrix[1,1].conjugate()).imag
    return [a,b,u,v,j]


def flavor_data(t, coeff):
    identity, basis = flavor["identity"], flavor["basis"]
    terms = (t, flavor["arcs"][32], flavor["arcs"][152], identity)
    action = sum((s*m for s,m in zip(coeff,terms)),acb_mat(2,2))
    mean = action.trace().real/2
    d = action-mean*identity
    r = (d[0,0].real**2+abs(d[0,1])**2).sqrt()
    assert r>0
    sinc = r.sin()/r
    phase = acb(0,mean).exp()
    wave = phase*(r.cos()*identity+acb(0,sinc)*d)
    matrix = flavor["radial"]+basis*wave*basis.transpose()
    derivatives = []
    for term in terms:
        mu = term.trace().real/2
        e = term-mu*identity
        dr = (d[0,0].real*e[0,0].real+(d[0,1].conjugate()*e[0,1]).real)/r
        ds = (r.cos()-sinc)*dr/r
        motion = acb(0,mu)*wave+phase*(-r.sin()*dr*identity+acb(0,ds)*d+acb(0,sinc)*e)
        change = basis*motion*basis.transpose()
        dat = [2*(matrix[i,j].conjugate()*change[i,j]).real for i,j in ((0,1),(1,2),(0,2),(2,0))]
        original = [matrix[0,1],matrix[1,2],matrix[0,2].conjugate(),matrix[1,1].conjugate()]
        varied = [change[0,1],change[1,2],change[0,2].conjugate(),change[1,1].conjugate()]
        product = acb(0)
        for i in range(4):
            value = varied[i]
            for j in range(4):
                if j != i:
                    value *= original[j]
            product += value
        derivatives.append([*dat,product.imag])
    return observations(matrix), derivatives


def tight_data(t, center, radii):
    base, _ = flavor_data(t,center)
    box = [v+arb(0,r) for v,r in zip(center,radii)]
    rough, derivatives = flavor_data(t,box)
    values = [a.intersection(b+sum(abs(derivatives[i][j])*arb(0,r) for i,r in enumerate(radii)))
              for j,(a,b) in enumerate(zip(rough,base))]
    return values, derivatives


def phase_system(weights, xs, point_rows, S):
    tau = sharp["tau"]/2
    means = [sum(w*row[j] for w,row in zip(weights,point_rows)) for j in range(5)]
    a,b,u,v,j = means
    cu = (1-S)*(1-2*u)+v-(1-S)*b
    cb = S-(1-S)*u
    delta = [0,-tau*cb,-tau*cu,tau*(1-u),1]
    phase = j-tau*((1-S)*u*(1-u)-(1-u)*v+(S-(1-S)*u)*b)
    residual = arb_mat([[sum(weights)-1],[sum(w*x for w,x in zip(weights,xs))],
                        [a-S*(1-u)],[phase]])
    matrix = arb_mat([[1]*4, xs,
                      [row[0]+S*row[2] for row in point_rows],
                      [sum(c*y for c,y in zip(delta,row)) for row in point_rows]])
    return residual, matrix, delta


def contact_weights(rows, S):
    from flint import arb_poly
    xs,obs = zip(*rows)
    a,b,u,v,j = zip(*obs)
    mat = arb_mat([[1]*4,xs,[a0-S*(1-u0) for a0,u0 in zip(a,u)],u])
    base = mat.solve(arb_mat([[1],[0],[0],[sharp["lo"]]]))
    speed = mat.solve(arb_mat([[0],[0],[0],[sharp["hi"]-sharp["lo"]]]))
    ws=[arb_poly([base[i,0],speed[i,0]]) for i in range(4)]
    U=arb_poly([sharp["lo"],sharp["hi"]-sharp["lo"]])
    B,V,J=[sum(w*y for w,y in zip(ws,col)) for col in (b,v,j)]
    eq=J-(sharp["tau"]/2)*((1-S)*U*(1-U)-(1-U)*V+(S-(1-S)*U)*B)
    lo, hi = arb(0),arb(1)
    assert eq(lo)*eq(hi)<0
    for _ in range(100):
        mid=(lo+hi)/2
        val=eq(mid)
        assert 0 not in val
        if val*eq(lo)>0:lo=mid
        else: hi=mid
    return [w(lo.union(hi)) for w in ws]


def prepare_case(share, support, cost, strength, tangent, second_angle=False):
    S = arb(".041775")**2 if second_angle else sharp["S"]
    def ordered(row):
        return [row[1],row[0],*row[2:]] if second_angle else row
    scale = -strength/tangent
    theta = arb(0).union(1)
    assert scale>0
    remaining = 1-share
    heights = [acb.zeta_zero(n).imag for n in support]
    inputs = [study["native"]["plane_response"](h) for h in heights]
    unshifted = [study["center"]/(h*h+study["offset"])-1 for h in heights]
    costs = [remaining*cost(x)[0] for x in unshifted]
    unshifted = [remaining*x for x in unshifted]
    xs = [x+theta*scale*y for x,y in zip(unshifted,costs)]
    middle = [x+scale*y/2 for x,y in zip(unshifted,costs)]
    center_input = study["native"]["plane_response"](acb.zeta_zero(11431).imag)
    def mixed(row, center):
        return ordered([share*b+remaining*a for a,b in zip(row,center)])
    center_obs = flavor_data(center_input,centers)[0]
    central = [mixed(flavor_data(t,centers)[0],center_obs) for t in inputs]
    weights = contact_weights(list(zip(middle,central)),S)
    mid = [w.mid() for w in weights]
    r = arb(2)**-28
    data = [tight_data(t,centers,[r]*4) for t in inputs]
    mid_box, mid_derivatives = tight_data(center_input,centers,[r]*4)
    values = [mixed(row,mid_box) for row,ds in data]
    derivatives = [[mixed(d,e) for d,e in zip(ds,mid_derivatives)] for row,ds in data]
    box = [w+arb(0,arb(".002")) for w in mid]
    residual, _, _ = phase_system(mid,xs,values,S)
    _, initial_jac, _ = phase_system(box,xs,values,S)
    pre = arb_mat([[initial_jac[i,j].mid() for j in range(4)] for i in range(4)]).inv()
    ident = arb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    contraction = ident-pre*initial_jac
    ratios = [sum(abs(contraction[i,j]) for j in range(4)) for i in range(4)]
    step0 = pre*residual
    print("contraction",max(ratios).str(10),flush=True)
    assert max(ratios)<1
    assert all(abs(step0[i,0])+ratios[i]*arb(".002") < arb(".002") for i in range(4))
    for k in range(6):
        _, jac, delta = phase_system(box,xs,values,S)
        step = jac.solve(residual)
        image = [w-step[i,0] for i,w in enumerate(mid)]
        if k==0:
            assert all(a.contains_interior(b) for a,b in zip(box,image))
        box = [a.intersection(b) for a,b in zip(box,image)]
    assert all(w>0 for w in box)
    _, jac, delta = phase_system(box,xs,values,S)
    multipliers = jac.transpose().solve(arb_mat([[y] for y in costs]))
    signs = []
    for j in range(4):
        mean = [sum(w*point[j][k] for w,point in zip(box,derivatives)) for k in range(5)]
        deriv = -multipliers[2,0]*(mean[0]+S*mean[2])
        deriv -= multipliers[3,0]*sum(a*b for a,b in zip(delta,mean))
        print("Derivative",j,deriv.str(15))
        assert 0 not in deriv
        signs.append(1 if deriv>0 else -1)
    expected_signs = [-1, 1, 1, 1] if second_angle else ([-1, 1, -1, -1] if share>0 else [-1]*4)
    assert signs == expected_signs
    print("Minimum corner signs",[-s for s in signs],"weights",*[a.str(10) for a in box])
    mean = [values[0][i]+sum(box[k]*(values[k][i]-values[0][i]) for k in (1,2,3)) for i in range(5)]
    reference = [mean[i] for i in (2,1,3,4)]
    assert sharp["lo"] < mean[2] < sharp["hi"]
    other_low,other_high = (".224955",".225065") if second_angle else (".041775",".041885")
    assert arb(other_low)**2 < mean[1]/(1-mean[2]) < arb(other_high)**2
    assert mean[4]>0
    endpoint_data = []
    determinant_signs = []
    for u in (sharp["lo"],sharp["hi"]):
        uc,bc,vc,jc = reference
        cu = (1-S)*(1-u-uc)+vc-(1-S)*bc
        cb = S-(1-S)*u
        def diff_features(row, base, cu=cu, cb=cb, u=u):
            a,b,k,v,j = [s-t for s,t in zip(row,base)]
            return a+S*k, j-sharp["tau"]/2*(cu*k-(1-u)*v+cb*b)
        mat = arb_mat([[1,x-xs[0],*diff_features(o,values[0])] if i else [1,0,0,0]
                       for i,(x,o) in enumerate(zip(xs,values))])
        rhs=arb_mat([[v-costs[0]] for v in costs])
        center_mat=arb_mat([[mat[i,j].mid() for j in range(4)] for i in range(4)])
        center_dual=center_mat.solve(rhs)
        dual=center_dual+mat.solve(rhs-mat*center_dual)

        assert 0 not in mat.det() and dual[2,0]>0 and dual[3,0]>0
        determinant_signs.append(mat.det()>0)
        a0,p0 = diff_features(ordered(mid_box),values[0])
        center_slack = -costs[0]+dual[1,0]*xs[0]-dual[2,0]*a0-dual[3,0]*p0
        if share>0:
            assert center_slack<0
        assert 1-theta*scale*dual[1,0]>0 and 1+(1-theta)*scale*dual[1,0]>0
        def normalized(obs):
            a,b,k,v,j=obs
            p=j-jc-sharp["tau"]/2*(cu*(k-uc)-(1-u)*(v-vc)+cb*(b-bc))
            return dual[2,0]/dual[3,0]*(a-S*(1-k))+p
        if second_angle:
            assert cu>0 and cb>0
            max_phase=1-jc+sharp["tau"]/2*(cu*uc+(1-u)*(1-vc)+cb*bc)
            outside_bound = dual[2,0]/dual[3,0]*arb(104)/867+max_phase
        else:
            outside_bound = sharp["compare_torus"](arb(1),normalized)
        upper = dual[3,0]*(remaining*outside_bound+share*normalized(ordered(mid_box)))
        candidate_mean = sum(w*y for w,y in zip(box,costs))
        for n in (8499,16001):
            x, _ = sharp["row"](acb.zeta_zero(n))
            value, dx = cost(x)
            value = remaining*value
            x = remaining*x
            derivative = dx/(1+theta*scale*dx)
            assert (derivative>dual[1,0] if n==8499 else derivative<dual[1,0])
            assert value-dual[1,0]*(x+theta*scale*value) > candidate_mean+upper
        endpoint_data.append((dual, diff_features))
    assert determinant_signs[0] == determinant_signs[1]
    corner = [center-sign*r for sign,center in zip(signs,centers)]
    corner_center = flavor_data(center_input,corner)[0]
    corner_rows=[mixed(flavor_data(t,corner)[0],corner_center) for t in inputs]
    for t in (arb(0),arb(1)/2,arb(1)):
        coord = [x+t*scale*y for x,y in zip(unshifted,costs)]
        corner_weights=contact_weights(list(zip(coord,corner_rows)),S)
        optimum=sum(w*v for w,v in zip(corner_weights,costs))
        if t == 0:
            expected = arb("1.74424814946858e-26") if share>0 else arb("5.868177693464756e-25")
            if second_angle:
                expected = arb("4.93918790720680305e-30")
            assert abs(optimum*strength-expected)<arb("1e-39")
        print("Central weight cap", share.str(3), "budget direction",t.str(3),"gauge tolerance",
              (t*scale*optimum).str(23),"port tolerance",((1-t)*strength*optimum).str(23),flush=True)
    slopes=[row[0]+S*row[2]-S for row in values]
    mat=arb_mat([[1]*4,xs,slopes,[row[2] for row in values]])
    base=mat.solve(arb_mat([[1],[0],[0],[sharp["lo"]]]))
    speed=mat.solve(arb_mat([[0],[0],[0],[sharp["hi"]-sharp["lo"]]]))
    from flint import arb_poly
    curves=[arb_poly([base[i,0],speed[i,0]]) for i in range(4)]
    u=arb_poly([sharp["lo"],sharp["hi"]-sharp["lo"]])
    b,v,j=[sum(w*row[k] for w,row in zip(curves,values)) for k in (1,3,4)]
    phase=j-(sharp["tau"]/2)*((1-S)*u*(1-u)-(1-u)*v+(S-(1-S)*u)*b)
    assert 0 not in phase.derivative()(arb(0).union(1))
    def lower_slack(x, observed, cost_at_x):
        value = remaining*cost_at_x-costs[0]
        coord_difference = remaining*x-unshifted[0]+theta*scale*value
        data = mixed(observed, mid_box)
        bounds = []
        for dual, features in endpoint_data:
            a, phase = features(data,values[0])
            bounds.append((value-dual[1,0]*coord_difference-dual[2,0]*a-dual[3,0]*phase).lower())
        return min(bounds)
    return lower_slack


def verify():
    ctx.prec = 224
    cost, strength, tangent = joint["port_cost"](include_tangent=True)
    cases = [(arb(0), sharp["support"], False), (arb(1)/2, (10762,11159,11432,11831), False),
             (arb(3)/4, (11406,11416,11430,11432), True)]
    checks = [(support,prepare_case(share,support,cost,strength,tangent,second)) for share,support,second in cases]
    radius = arb(2)**-28
    coeff_box = [c+arb(0,radius) for c in centers]
    floor=arb(1)
    for start in range(8500,16001,256):
        count = min(256,16001-start)
        zeros = acb.zeta_zeros(start,count)
        assert len(zeros)==count
        for n, zero in enumerate(zeros,start):
            if n==11431:
                continue
            h = zero.imag
            x = study["center"]/(h*h+study["offset"])-1
            value = cost(x)[0]
            t = study["native"]["plane_response"](h)
            obs = observations(flavor["mixing"](t,coeff_box))
            for support, slack in checks:
                if n in support:
                    continue
                gap=slack(x,obs,value)
                if not gap>0:
                    gap=slack(x,tight_data(t,centers,[radius]*4)[0],value)
                assert gap>0, (n,str(gap))
                floor=min(floor,gap)
    print("Uniform strict slack for all three occupation thresholds:",floor.str(14))


if __name__ == "__main__":
    verify()
