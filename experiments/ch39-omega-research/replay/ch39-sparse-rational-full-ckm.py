"""用稀疏有理谱响应在一个真实壳上实现完整中心夸克观测。"""
from pathlib import Path
from runpy import run_path
from itertools import combinations
from flint import arb, acb, acb_mat, arb_mat, arb_series, ctx


def setup(root,n):
    native=run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    model=run_path(str(root/"ch39-complex-four-site-unit-quark-amplitudes.py"))
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    order=run_path(str(root/"ch39-bounded-common-rational-amplitude.py"))["ordered"]
    ckm=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    ps=native["native"](root,n)
    ctx.prec=256
    h=acb.zeta_zero(n).imag
    nodes=[order(shell,q,h) for q in (-4,2)]
    v=acb_mat(ckm["standard_matrix"](*[arb(i).sin() for i in (".22704",".04275",".00377")],arb("1.139")))
    mass=[list(map(arb,s)) for s in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    A,B=[model["diag"]([x/row[2] for x in row]) for row in mass]
    A=model["adj"](v)*A*v
    coords=[[(p[0,0].real,p[1,1].real,p[0,1].real,p[0,1].imag,p[2,2].real+p[2,3].imag,p[2,2].real-p[2,3].imag) for p in rows] for rows in ps]
    return native,model,A,B,mass,nodes,coords


# 谱序先按荷 −4 的十九层、再按荷 2 的十九层排列；程序索引从零起。
CASES=(
    (2306888100085,".123",(9,28,16,35),(2,3),
     ("3.890343e-6",".880742900",".28563560",".984586890",".128542358",
      ".004617737",".037612878","1.41604533e-8","5.03638969e-8")),
    (2306888100086,".0499143",(9,28,16,35),(2,3),
     ("3.866638e-6",".43402028",".128611824",".924019496",".296581452",
      ".002002098",".040199762","2.60150506e-9","1.89350008e-7")),
)


def frequency_data(nodes,coords,picks):
    all_nodes=nodes[0]+nodes[1]
    centers=[all_nodes[i] for i in picks]
    assert all(not (a-b).contains(0) for a,b in combinations(centers,2))
    data=[]
    for q,spec in enumerate(nodes):
        sector=[]
        for i,energy in enumerate(spec):
            index=19*q+i
            if index in picks:
                sector.append((coords[q][i],picks.index(index),None))
            else:
                gaps=[energy-z for z in centers]
                assert all(not z.contains(0) for z in gaps)
                sector.append((coords[q][i],None,[1/g**2 for g in gaps]))
        data.append(sector)
    return data


def native_moments(variables,data,wide):
    count=len(variables)-3-len(wide)
    background=variables[0]
    amplitudes=variables[1:1+count]
    u,v=variables[1+count:3+count]
    widths=[arb(2)**-80 for _ in amplitudes]
    for index,weight in zip(wide,variables[3+count:],strict=True):
        widths[index]=weight
    output=[]
    for sector in data:
        result=[arb(0)]*6
        for column,anchor,gaps in sector:
            if anchor is not None:
                value=amplitudes[anchor]
            else:
                shares=[a*b for a,b in zip(widths,gaps,strict=True)]
                value=(background+dot(shares,amplitudes))/(1+sum(shares))
            result=[z+value*c for z,c in zip(result,column,strict=True)]
        output.append(result)
    return output,u,v


def dot(a,b):
    return sum(x*y for x,y in zip(a,b,strict=True))


def equations(variables,data,target,wide):
    native,u,v=native_moments(variables,data,wide)
    return block_invariants(native,u,v,target)


def block_invariants(native,u,v,target):
    q=u*v
    b=(1-q).sqrt()
    invariant=[]
    raw_vectors=[]
    goal_vectors=[]
    for (x,y,re,im,minus,plus),(p,s,T,kr,ki,w) in zip(native,target,strict=True):
        d=(T-q*s-(u-v)*w)/(1-q)
        invariant.extend((plus-s+w/u,minus-s-w/v,x+y-p-d))
        raw_vectors.append((2*re,-2*im,x-y))
        goal_vectors.append((2*kr/b,-2*ki/b,p-d))
    invariant.extend(dot(raw_vectors[i],raw_vectors[j])-dot(goal_vectors[i],goal_vectors[j])
                     for i,j in ((0,0),(1,1),(0,1)))
    return invariant


def value_jac(params,data,target,wide):
    values=equations(params,data,target,wide)
    columns=[]
    for j in range(len(params)):
        jet=[arb_series([x,int(i==j)]) for i,x in enumerate(params)]
        columns.append([z[1] for z in equations(jet,data,target,wide)])
    return arb_mat([[x] for x in values]),arb_mat(columns).transpose()


def cross(a,b):
    return [a[(i+1)%3]*b[(i+2)%3]-a[(i+2)%3]*b[(i+1)%3] for i in range(3)]


def triad(vectors):
    a,b=vectors
    length=dot(a,a).sqrt()
    e=[z/length for z in a]
    coefficient=dot(e,b)
    other=[z-coefficient*w for z,w in zip(b,e,strict=True)]
    length=dot(other,other).sqrt()
    f=[z/length for z in other]
    return arb_mat([e,f,cross(e,f)]).transpose()


def reproduce(root,box,data,targets,wide,mass,model):
    moments,u,v=native_moments(box,data,wide)
    a=(u*v).sqrt()
    b=(1-u*v).sqrt()
    native_vectors=[]
    target_vectors=[]
    raw=[]
    for (i,j,re,im,minus,plus),(p,s,T,kr,ki,w) in zip(moments,targets,strict=True):
        raw.append(acb_mat([[i,acb(re,im),0,0],[acb(re,-im),j,0,0],
                            [0,0,plus,0],[0,0,0,minus]]))
        d=(T-u*v*s-(u-v)*w)/(1-u*v)
        native_vectors.append((2*re,-2*im,i-j))
        target_vectors.append((2*kr/b,-2*ki/b,p-d))
    rotation=triad(native_vectors)*triad(target_vectors).transpose()
    w=(1+sum(rotation[i,i] for i in range(3))).sqrt()/2
    assert w>0
    x,y,z=[(rotation[(i+2)%3,(i+1)%3]-rotation[(i+1)%3,(i+2)%3])/(4*w) for i in range(3)]
    U=acb_mat([[acb(w,-z),acb(-y,-x),0,0],[acb(y,-x),acb(w,z),0,0],[0,0,1,0],[0,0,0,1]])
    f,g=(v/(u+v)).sqrt(),(u/(u+v)).sqrt()
    E=U*acb_mat([[1,0,0],[0,0,b],[0,g,-a*f],[0,-f,-a*g]])
    adj=model["adj"]
    assert (adj(E)*E-model["diag"]([1]*3)).contains(acb_mat(3,3))
    flavor=run_path(str(root/"ch39-three-gram-source-flavor.py"))
    ps=[]
    for four,target,spec in zip(raw,targets,mass,strict=True):
        p,s,T,kr,ki,w=target
        exact=acb_mat([[p,0,acb(kr,ki)],[0,s,w],[acb(kr,-ki),w,T]])
        effective=adj(E)*four*E
        assert (effective-exact).contains(acb_mat(3,3))
        levels,proj=flavor["projections"](effective*effective)
        assert all(z.overlaps(m*m) for z,m in zip(levels,spec,strict=True))
        ps.append(proj)
    ref=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    V=acb_mat(ref["standard_matrix"](*[arb(i).sin() for i in (".22704",".04275",".00377")],arb("1.139")))
    pu,pd=ps
    assert all((pu[i]*pd[j]).trace().real.overlaps(abs(V[i,j])**2) for i in range(3) for j in range(3))
    j=(pu[0]*pd[0]*pu[1]*pd[1]).trace().imag
    ref_j=(V[0,0]*V[1,1]*V[0,1].conjugate()*V[1,0].conjugate()).imag
    assert j>0 and j.overlaps(ref_j)
    print("完整重建的有符号 J",j.str(12),flush=True)


def verify(root,case):
    n,t,picks,wide,start=case
    native,model,A,B,mass,nodes,coords=setup(root,n)
    ctx.cap=2
    basis=native["target_basis"](A,arb(t)*B)*model["diag"]([1,1,acb(0,-1)])
    target=[]
    for matrix,spec in zip((A,B),mass,strict=True):
        m=model["adj"](basis)*matrix*basis*spec[2]
        assert m[0,1].contains(0) and m[1,2].imag.contains(0)
        target.append((m[0,0].real,m[1,1].real,m[2,2].real,m[0,2].real,m[0,2].imag,m[1,2].real))
    data=frequency_data(nodes,coords,picks)
    point=list(map(arb,start))
    assert len(point)==9
    for _ in range(5):
        values,j=value_jac(point,data,target,wide)
        step=j.solve(values).entries()
        point=[(x-dx).mid() for x,dx in zip(point,step,strict=True)]
    radius=arb("1e-18")
    box=[x+arb(0,radius.upper()) for x in point]
    values,j=value_jac(point,data,target,wide)
    pre=j.inv()
    _,jac=value_jac(box,data,target,wide)
    defect=arb_mat([[int(i==j) for i in range(9)] for j in range(9)])-pre*jac
    gain=max(sum(abs(x).upper() for x in row) for row in defect.tolist())
    shift=max(abs(x).upper() for x in (pre*values).entries())
    assert shift+gain*radius<radius,(gain,shift)
    count=len(picks)
    assert all(0<x<1 for x in box[:1+count])
    u,v=box[count+1:count+3]
    assert u>0 and v>0 and u*v<1
    assert all(x>0 for x in box[count+3:])
    reproduce(root,box,data,target,wide,mass,model)
    print("编号、次数上界、收缩常数、根中心",n,2*len(picks),gain.str(7),*[x.str(13) for x in box],flush=True)


if __name__=="__main__":
    for case in CASES:
        verify(Path(__file__).resolve().parent,case)
