"""Fit three mixing-angle centers on one branch, then bound its CP phase."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat

root = Path(__file__).parent
fixed = run_path(str(root/"ch39-fixed-local-flavor.py"))
single = run_path(str(root/"ch39-single-contact-ckm-box.py"))
identity, basis, radial = (fixed[k] for k in ("identity","basis","radial"))
shell = fixed["shell"]
s,t,u,delta0 = shell["centers"]
targets = (s*s*(1-u*u),t*t*(1-u*u),u*u)
indices = ((0,1),(1,2),(0,2))
q32,q152 = fixed["arcs"][32],fixed["arcs"][152]
k32,k152 = (1-arb(d)/3468 for d in (32,152))
native = fixed["native"]["plane_response"](acb.zeta_zero(11431).imag)
contact = ((1+3*k32)*q32+(1+3*k152)*q152)/4


def mixing_matrix(params, derivatives=False, response=native):
    gain,height,phase = params
    action = gain*response+height*contact+phase*identity
    mean = action.trace().real/2
    d = action-mean*identity
    gap = (d[0,0].real**2 + abs(d[0,1])**2).sqrt()
    sine,cosine = gap.sin_cos()
    sinc = sine/gap
    ephase = shell["cis"](mean)
    iu = acb(0,1)
    plane = ephase*(cosine*identity+iu*sinc*d)
    w = radial+basis*plane*basis.transpose()
    if not derivatives:
        return w
    derivatives_w = []
    for direction in (response,contact,identity):
        dm = direction.trace().real/2
        dep = direction-dm*identity
        dg = (d[0,0].real*dep[0,0].real + (d[0,1].conjugate()*dep[0,1]).real)/gap
        ds = dg*(cosine-sinc)/gap
        derivative = iu*dm*plane+ephase*(-sine*dg*identity + iu*(ds*d+sinc*dep))
        dw = basis*derivative*basis.transpose()
        derivatives_w.append(dw)
    return w, derivatives_w


def observables(params, derivatives=False, response=native):
    if not derivatives:
        return mixing_matrix(params, response=response)
    w, derivatives_w = mixing_matrix(params, True, response)
    cols = [[2*(w[i,j].conjugate()*dw[i,j]).real for i,j in indices] for dw in derivatives_w]
    residual = arb_mat([[abs(w[i,j])**2-z] for (i,j),z in zip(indices,targets)])
    return residual,arb_mat(cols).transpose()


def center_fit(seed, response=native):
    point = list(map(arb,seed))
    for _ in range(3):
        f,jac = observables(point,True,response)
        step = jac.solve(f)
        point = [(v-step[i,0]).mid() for i,v in enumerate(point)]
    r = arb("1e-20")
    box = [v+arb(0,r.upper()) for v in point]
    f,jac = observables(point,True,response)
    a = jac.inv()
    _,vary = observables(box,True,response)
    rem = arb_mat([[int(i==j) for j in range(3)] for i in range(3)])-a*vary
    lip = max(sum(abs(rem[i,j]) for j in range(3)).upper() for i in range(3))
    resid = a*f
    assert lip<1 and max(abs(resid[i,0]).upper() for i in range(3))+lip*r<r
    return box


def phase_readout(params,response=native):
    w = observables(params,response=response)
    ycp = (w[0,0]*w[1,0].conjugate()*w[1,1]*w[0,1].conjugate()).imag
    a,b,c = ((1-z*z).sqrt() for z in (s,t,u))
    sd = ycp/(a*b*c*c*s*t*u)
    cd = (abs(w[1,0])**2-(s*b)**2-(a*t*u)**2)/(2*a*b*s*t*u)
    return acb(cd,sd).arg(),ycp


def verify():
    for copies in range(1,5):
        for through in range(copies+1):
            if copies == 4 and through == 3:
                continue
            bypass = copies-through
            term = (bypass+through*k32)*q32+(bypass+through*k152)*q152
            single["exclusion_count"](11431,[term])
    params = center_fit((".7210227965417486","2.557578307275763","3.086954654848417"))
    delta,ycp = phase_readout(params)
    assert abs(delta-delta0) < shell["phase_radius"]
    print("minimal equal channels 4, one bypass and three spin comparisons;",
          "three-angle-fit phase",delta.str(15),"J",ycp.str(15))
    print("locally fixed parameters",*(p.str(15) for p in params))
    opposite = center_fit(("7.33802792808139","15.684754005933137",".8501645181600139"))
    neg,jneg = phase_readout(opposite)
    assert arb("-1.205") < neg < arb("-1.204") and jneg < 0
    print("same three angles on another positive branch:",neg.str(15),jneg.str(15))


if __name__ == "__main__":
    verify()
