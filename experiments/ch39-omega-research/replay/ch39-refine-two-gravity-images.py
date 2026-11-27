"""在已证荷六三读数的参数盒中，用相同热方程收紧两条高源像。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def high_images(root):
    joint = run_path(str(root/"ch39-one-generator-up-and-gravity.py"))
    rank = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    shells = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    pivots = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    spectrum = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    grav = run_path(str(root/"ch39-three-index-free-atoms.py"))
    orbits = run_path(str(root/"ch39-up-quark-arbitrary-site-phase-floor.py"))["ORBITS"]
    _, box, gev, (low,h,k) = joint["verify"](root,return_parameters=True)
    n = rank["N"]+2
    moments,heat,fields = grav["broken_fields"](root,[n])
    field = fields[0]
    _,_,cross,_ = next(x for x in pivots["shell_weights"](root,n) if x[0]==6)
    c = [z.real for z in cross[0]]
    ctx.prec = 384
    src,target = (spectrum(shells,6,j)[1] for j in (n-1,n))
    E = acb_mat([[1/arb(len(o)).sqrt() if i in o else 0 for o in orbits] for i in range(19)])
    P = rank["projector"](orbits)
    col,adj,inner = joint["col"],rank["adj"],joint["inner"]
    fixed,off,small,big,tangent,rho = joint["mode_parameters"](P,E,src,target)
    e = ((gev/arb("1.22089e19")-low)/(h-low)).sqrt()
    trivial = [j for j in range(19) if inner(adj(E)*col(target,j),adj(E)*col(target,j)).real>arb(".5") and j!=15]
    def image(a):
        u = off+small*e+big*a
        t = (rho-inner(u,u).real).sqrt()
        return fixed+u+tangent*t, big-tangent*(inner(u,big).real/t)
    def equations(ab):
        a,b = ab
        y,yp = image(a)
        amps = {}
        for j in trivial:
            v = col(target,j)
            f = inner(v,y)
            amps[j] = (low+(h-low)*abs(f)**2,2*(h-low)*(f.conjugate()*inner(v,yp)).real,0)
        amps[2] = (low+(k-low)*(1-b*b),0,-2*(k-low)*b)
        amps[18] = (low+(k-low)*b*b,0,2*(k-low)*b)
        quiet = sum(c[j] for j in range(19) if j not in amps and j!=15)
        small_x = max(low.upper(),(gev/arb("1.22089e19")).upper())
        mass = abs(quiet)+abs(c[15])
        errors = [mass*small_x**2*sum(abs(v) for _,v,_ in field)/(2*arb.pi()),
                  mass*small_x**4*sum(abs(g) for g,_,_ in field)/(8*arb.pi())]
        values = [arb(0,r.upper())-v for r,v in zip(errors,(1,arb("1e-122")),strict=True)]
        jac = arb_mat(2,2)
        for j,(y,ya,yb) in amps.items():
            x = y**2
            thermal = moments(field,heat,x)
            decay = [(-z/x).exp() for _,_,z in field]
            dy = [y*sum(v*q for (_,v,_),q in zip(field,decay,strict=True))/arb.pi(),
                  -y**3*sum(g*q for (g,_,_),q in zip(field,decay,strict=True))/(2*arb.pi())]
            for i in range(2):
                values[i] += c[j]*thermal[i]
                jac[i,0] += c[j]*dy[i]*ya
                jac[i,1] += c[j]*dy[i]*yb
        return arb_mat([[v] for v in values]),jac

    for _ in range(3):
        point = [v.mid() for v in box]
        val,jac = equations(point)
        pre = jac.mid().inv()
        _,wide = equations(box)
        remainder = arb_mat([[1,0],[0,1]])-pre*wide
        residual = pre*val
        box = [v.intersection(x-residual[i,0]+sum(remainder[i,j]*(box[j]-point[j]) for j in range(2)))
               for i,(v,x) in enumerate(zip(box,point,strict=True))]
    a,b = box
    y,_ = image(a)
    q = col(src,18)
    t3,t19 = [joint["align"](col(target,j),q) for j in (2,18)]
    z = t3*(1-b*b).sqrt()+t19*b
    print("收紧后的第18、19层参数",a.str(24),b.str(24),flush=True)
    return col(src,17),q,y,z


if __name__=="__main__":
    high_images(Path(__file__).resolve().parent)
