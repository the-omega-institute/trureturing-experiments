"""由最短时钟的 CKM 读数辨识四个实味系数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat

ordered = run_path(str(Path(__file__).with_name("ch39-ordered-flavor-unitary.py")))
flavor = ordered["flavor"]


def step_with_derivatives(terms, duration, coeff):
    ident = flavor["identity"]
    action = sum((c*v for c,v in zip(coeff,terms)),acb_mat(2,2))
    mean = action.trace().real/2
    reduced = action-mean*ident
    gap = (reduced[0,0].real**2+abs(reduced[0,1])**2).sqrt()
    assert gap>0
    angle = duration*gap
    phase = acb(0,duration*mean).exp()
    sine,cosine = angle.sin(),angle.cos()
    amplitude = sine/gap
    wave = phase*(cosine*ident+acb(0,amplitude)*reduced)
    derivatives = []
    for term in terms:
        trace = term.trace().real/2
        change = term-trace*ident
        speed = (reduced[0,0].real*change[0,0].real+
                 (reduced[0,1].conjugate()*change[0,1]).real)/gap
        deriv = -duration*sine*speed*ident + acb(0,amplitude)*change
        deriv += acb(0,(duration*cosine-amplitude)*speed/gap)*reduced
        derivatives.append(acb(0,duration*trace)*wave+phase*deriv)
    return wave,derivatives


def observations(inputs,coeff,time_slopes=None):
    wave = flavor["identity"]
    variations = [acb_mat(2,2) for _ in range(4)]
    if time_slopes is not None:
        variations.append(acb_mat(2,2))
    for i,(terms,duration) in enumerate(inputs):
        step,ds = step_with_derivatives(terms,duration,coeff)
        if time_slopes is not None:
            action=sum((c*v for c,v in zip(coeff,terms)),acb_mat(2,2))
            ds.append(acb(0,time_slopes[i])*action*step)
        variations = [step*v+d*wave for v,d in zip(variations,ds)]
        wave = step*wave
    basis = flavor["basis"]
    v = flavor["radial"]+basis*wave*basis.transpose()
    pairs = ((0,1),(1,2),(0,2),(1,1))
    def factors(matrix):
        a,b,u,c = [matrix[i,j] for i,j in pairs]
        return [a,b,u.conjugate(),c.conjugate()]
    def product(row):
        a,b,c,d = row
        return a*b*c*d
    amplitudes = factors(v)
    values = [abs(z)**2 for z in amplitudes[:3]]+[product(amplitudes).imag]
    jac = []
    for d in variations:
        change = factors(basis*d*basis.transpose())
        jdot = sum(product(amplitudes[:i]+[z]+amplitudes[i+1:]).imag for i,z in enumerate(change))
        jac.append([2*(z.conjugate()*dz).real for z,dz in zip(amplitudes[:3],change[:3])]+[jdot])
    return values,arb_mat(jac).transpose()


def verify():
    ordered["previous"]["setup"]()
    inputs = []
    for n,k in ((11410,10985),(11431,676974),(11432,327502)):
        t = ordered["study"]["native"]["plane_response"](acb.zeta_zero(n).imag)
        terms = (t,flavor["arcs"][32],flavor["arcs"][152],flavor["identity"])
        inputs.append((terms,arb(k)/1015461))
    centers = ordered["previous"]["box"]["centers"]
    side = arb(2)**-28
    base,central_jac = observations(inputs,centers)
    pre = arb_mat([[central_jac[i,j].mid() for j in range(4)] for i in range(4)]).inv()
    cube = [c+arb(0,side) for c in centers]
    ranges,jac = observations(inputs,cube)
    identity = arb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    remainder = identity-pre*jac
    eta = max(sum(abs(remainder[i,j]) for j in range(4)).upper() for i in range(4))
    assert eta<arb("0.000782"),eta.str(12)
    normalized = ((1-eta)*side).lower()
    print("Sup row norm of I - preconditioned derivative",eta.str(14))
    print("Guaranteed centered observation radius, in preconditioned coordinates",normalized.str(15))
    print("Original determinant",jac.det().str(15))
    assert 0 not in jac.det()
    # 书中的四个中央目标给出精确的 A,B,U,J。
    s,t,v,angle = flavor["shell"]["centers"]
    u = v*v
    flavor_matrix=flavor["identity"]
    for terms,duration in inputs:
        step,_=step_with_derivatives(terms,duration,cube)
        flavor_matrix=step*flavor_matrix
    basis=flavor["basis"]
    matrix=flavor["radial"]+basis*flavor_matrix*basis.transpose()
    a,b,z,j=ranges
    quartet=matrix[0,1]*matrix[1,2]*matrix[0,2].conjugate()*matrix[1,1].conjugate()
    obs=[(a/(1-z)).sqrt(),(b/(1-z)).sqrt(),z.sqrt(),(quartet+a*b*z/(1-z)).arg()]
    widths=[arb(".000055")]*3+[arb(".006")]
    errors=[abs(v-w)/tolerance for v,w,tolerance in zip(obs,(s,t,v,angle),widths)]
    assert errors[1]>max(errors[i] for i in (0,2,3)) and obs[1]<t
    signs=[]
    for k in range(4):
        deriv=jac[1,k]*(1-z)+b*jac[2,k]
        assert 0 not in deriv
        signs.append(1 if deriv>0 else -1)
    optimal=[c+sgn*side for c,sgn in zip(centers,signs)]
    _,bb,zz,_=observations(inputs,optimal)[0]
    miss=(t-(bb/(1-zz)).sqrt())/widths[1]
    print("Best relative error within calibration box",miss.str(19),"corner",signs)
    assert arb(".83")<miss<arb(".84")
    target = [s*s*(1-u),t*t*(1-u),u,
              s*t*v*((1-s*s)*(1-t*t)).sqrt()*(1-u)*angle.sin()]
    guess = list(centers)
    for _ in range(6):
        vals,jac0 = observations(inputs,guess)
        change = jac0.solve(arb_mat([[a-b] for a,b in zip(vals,target)]))
        guess = [(v-change[k,0]).mid() for k,v in enumerate(guess)]
    vals,jac0 = observations(inputs,guess)
    pre = arb_mat([[jac0[i,j].mid() for j in range(4)] for i in range(4)]).inv()
    width = arb("1e-12")
    region = [v+arb(0,width) for v in guess]
    _,jac_region = observations(inputs,region)
    rem = identity-pre*jac_region
    rate = max(sum(abs(rem[i,j]) for j in range(4)).upper() for i in range(4))
    residual = pre*arb_mat([[a-b] for a,b in zip(vals,target)])
    assert max(abs(residual[i,0]) for i in range(4))+rate*width<width
    print("Coefficients for central targets",*[v.str(17) for v in region])
    print("Changes from prescribed centers",*[(v-c).str(14) for v,c in zip(region,centers)])
    # 区分同一 J 对应的相位分支。由正实部取第一象限内的解。
    wave = flavor["identity"]
    for terms,duration in inputs:
        step,_ = step_with_derivatives(terms,duration,region)
        wave = step*wave
    basis = flavor["basis"]
    v = flavor["radial"]+basis*wave*basis.transpose()
    quartet = v[0,1]*v[1,2]*v[0,2].conjugate()*v[1,1].conjugate()
    assert (quartet+target[0]*target[1]*u/(1-u)).real>0



if __name__=="__main__":
    verify()
