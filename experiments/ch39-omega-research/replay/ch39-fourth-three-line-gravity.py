"""由两条源端平凡线和一条目的端非平凡线同时校准第四壳的两项热迹与弱差。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat


def adj(matrix):
    return matrix.conjugate().transpose()


def column(matrix, j):
    phase = matrix[0,j]
    assert abs(phase)>0
    return [matrix[i,j]*phase.conjugate()/abs(phase) for i in range(matrix.nrows())]


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    _,_,incoming,readout,(low,h,k),weights,(moments,heat,fields) = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    shell = lib("three-ray-chiral-boundary")
    n = lib("rank-two-generator-up-quark")["N"]
    mid = lib("adjacent-pvm-transport")["spectrum"](shell,6,n+3)[1]
    target = adj(readout)
    sites = shell["geo"]["sites"]
    cycle = [sites.index((y,z,x)) for x,y,z in sites]
    groups = lib("charge-six-projector-algebra-transport")["groups"]
    assert {0,15}<=set(groups(mid,cycle)[0]) and 18 not in groups(target,cycle)[0]

    lines = acb_mat([column(mid,0),column(mid,15),column(target,18)]).transpose()
    eye = lambda k: acb_mat([[int(i==j) for j in range(k)] for i in range(k)])
    assert (adj(lines)*lines).contains(eye(3))
    high = acb_mat([column(incoming,0),column(incoming,1)]).transpose()*acb_mat([[(h-low).sqrt(),0],[0,(k-low).sqrt()]])
    row16 = column(target,15)
    for j in range(19):
        readout[15,j] = row16[j].conjugate()
    coords = readout*lines
    a = adj(lines)*high
    fixed = readout*high-coords*a
    u = acb_mat([[coords[15,j] for j in range(3)]])
    r = (u*adj(u))[0,0].real.sqrt()
    v = adj(u)*(1/r)
    # 第三条线属于非平凡字符，用于 u 的正交补。另一条由第一坐标正交化得到。
    assert u[0,2].contains(0)
    first = acb_mat([[1],[0],[0]])-v*v[0,0].conjugate()
    first *= 1/(adj(first)*first)[0,0].real.sqrt()
    basis = acb_mat([[v[j,0],first[j,0],int(j==2)] for j in range(3)])
    assert (adj(basis)*basis).contains(eye(3))
    perp = coords*acb_mat([[basis[j,1],basis[j,2]] for j in range(3)])
    old16 = acb_mat([[fixed[15,0],fixed[15,1]]])
    mz,mp = arb("91.1876"),arb("1.22089e19")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    weak_low = lib("native-rh-ray-pvm")["weak_series"](56)
    weak_high = lib("full-neutral-free-cutoff")["weak_column"]
    def scan(s):
        return weights[15]*weak_low(s)+(1-weights[15])*arb("-.074256577194663365545")
    left,right = arb("1.28"),arb(2)
    assert scan(left)>0>scan(right)
    for _ in range(42):
        point = (left+right)/2
        if scan(point)>0:
            left = point
        else:
            assert scan(point)<0
            right = point
    left -= arb("1e-11")
    right += arb("1e-11")
    times = left.union(right)

    def system(angles,s):
        cutoff = mz*s.sqrt()/mp
        row = (acb_mat([[(cutoff-low).sqrt(),0]])-old16)*(1/r)
        gram = adj(a)*a-adj(row)*row
        d0 = gram[0,0].real.sqrt()
        upper = gram[0,1]/d0
        schur = gram[1,1].real-abs(upper)**2
        assert d0>0 and schur>0
        chol = acb_mat([[d0,upper],[0,schur.sqrt()]])
        assert (adj(chol)*chol).contains(gram) or all(
            z.contains(0) for z in (adj(chol)*chol-gram).entries())
        phase = acb(0,angles[0]).exp()
        sn,cs = angles[1].sin_cos()
        rot = acb_mat([[phase*cs,phase*sn],[-sn,cs]])
        deriv = (acb_mat([[acb(0,1)*phase*cs,acb(0,1)*phase*sn],[0,0]]),
                 acb_mat([[-phase*sn,phase*cs],[-cs,-sn]]))
        out = fixed+coords*v*row+perp*rot*chol
        change = [perp*d*chol for d in deriv]
        amps = [low+abs(out[i,0])**2+abs(out[i,1])**2 for i in range(19)]
        assert amps[15].overlaps(cutoff)
        amps[15] = cutoff
        # 固定第16层在两个角度上的偏导数均为零，低截止热核由绝对值外包。
        errors = [abs(weights[15])*cutoff**2*sum(abs(a) for _,a,_ in fields)/(2*arb.pi()),
                  abs(weights[15])*cutoff**4*sum(abs(g) for g,_,_ in fields)/(8*arb.pi())]
        thermal = [arb(0,err.upper())-goal for err,goal in zip(errors,(arb(1),arb("1e-122")),strict=True)]
        jac = arb_mat(2,2)
        for i,f in enumerate(amps):
            if i==15:
                continue
            assert f>arb("1e-6")
            for j,value in enumerate(moments(fields,heat,f*f)):
                thermal[j] += weights[i]*value
            slopes = [2*sum((out[i,j].conjugate()*d[i,j]).real for j in (0,1)) for d in change]
            decay = [(-m/f**2).exp() for _,_,m in fields]
            dk = f*sum(a*e for (_,a,_),e in zip(fields,decay,strict=True))/arb.pi()
            dw = -f**3*sum(g*e for (g,_,_),e in zip(fields,decay,strict=True))/(2*arb.pi())
            for j,rate in enumerate((dk,dw)):
                for p,df in enumerate(slopes):
                    jac[j,p] += weights[i]*rate*df
        return arb_mat([[z] for z in thermal]),jac,amps,schur

    center = [arb("-.30239"),arb("3.19627")]
    for _ in range(4):
        val,jac,_,_ = system(center,times.mid())
        step = jac.mid().solve(val.mid())
        center = [(p-d).mid() for p,d in zip(center,step.entries(),strict=True)]
    width = arb("1e-10")
    box = [v+arb(0,width.upper()) for v in center]
    val,jac,_,_ = system(center,times)
    pre = jac.mid().inv()
    assert not pre.det().contains(0)
    _,wide,amps,schur = system(box,times)
    defect = arb_mat([[1,0],[0,1]])-pre*wide
    q = max(sum(abs(v).upper() for v in row) for row in defect.tolist())
    shift = max(abs(v).upper() for v in (pre*val).entries())/width
    assert q+shift<1,(q,shift)
    other = sum(weights[i]*weak_high(f*f,p0,p1) for i,f in enumerate(amps) if i!=15)
    ends = [weights[15]*weak_low(s)+other for s in (left,right)]
    assert ends[0]>0>ends[1]
    print("两个角度、Schur 补",*[v.str(14) for v in (*center,schur)],flush=True)
    print("第16层截止 GeV、第18和19层振幅",*[v.str(13) for v in (mz*times.sqrt(),amps[17],amps[18])],flush=True)
    print("位移占比、导数行和、弱差端点",*[v.str(13) for v in (shift,q,*ends)],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
