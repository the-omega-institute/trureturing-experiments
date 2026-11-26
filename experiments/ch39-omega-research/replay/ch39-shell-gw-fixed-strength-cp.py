"""The whole CKM phase circle at fixed strength, six site norms, and scalar spectrum."""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb, acb_mat, arb

site=run_path(str(Path(__file__).with_name("ch39-shell-gw-site-isospectral.py")))
common,model,adj,hs,rows=[site[k] for k in ("common","model","adj","hs","rows")]
strength=arb(25)/2


def trial(P,target):
    A=model["array"](P)
    T=A.conj().T@np.linalg.inv(A@A.conj().T)
    N=np.eye(54)-T@A
    x=model["array"](site["numeric_trial"](P,target)).ravel()
    rng=np.random.default_rng(18)
    trial=N@(rng.normal(size=54)+1j*rng.normal(size=54))
    x+=trial/np.linalg.norm(trial)
    for _ in range(25):
        blocks=x.reshape(6,9)
        residual=np.sum(abs(blocks)**2,axis=1)-float(strength)**2/6
        if max(abs(residual))<1e-12:
            break
        D=np.zeros((6,54),complex)
        for j in range(6):
            D[j,9*j:9*(j+1)]=2*blocks[j].conj()
        D=D@N
        x-=N@D.conj().T@np.linalg.solve((D@D.conj().T).real,residual)
    return acb_mat([[acb(float(z.real),float(z.imag))] for z in x])


def check(index, operator_error=arb(0), target_error=None):
    F=common["old"]["F"] if index==43 else common["operators"]()[1]
    t=common["native"]["relative_phase"](arb(676))
    P=acb_mat(rows(F)+rows(common["conjugate"](F,t)))
    target=common["arc"]["central_target"]()
    z=trial(P,target)
    T=adj(P)*(P*adj(P)).inv()
    x=z+T*(target-P*z)
    N=model["eye"](54)-T*P
    D=acb_mat([[2*x[j,0].conjugate() if j//9==i else 0 for j in range(54)] for i in range(6)])
    V=N*adj(D)
    gram=D*V
    R=acb_mat([[gram[i,j].real for j in range(6)] for i in range(6)])
    C=V*R.inv()
    masses=[sum(abs(x[9*i+j,0])**2 for j in range(9)) for i in range(6)]
    residual=sum((mass-strength**2/6)**2 for mass in masses).sqrt()
    ca=common["count"]
    assert ca(rows(P*adj(P)),arb(1)/400)==0
    # Only the two off-diagonal entries of the standard 1-3 rotation vary.
    # Right multiplication by R12 Dd weights their rows by the down masses.
    d,s,b=map(arb,(".0000135",".000268",".01401"))
    angle=arb(".22704")
    change=2*arb(".00377").sin()*(b*b+(d*angle.cos())**2+(s*angle.sin())**2).sqrt()
    if target_error is not None:
        change=target_error
    assert 20*operator_error<1
    inverse=20/(1-20*operator_error)
    shift=inverse*operator_error
    error=inverse*(change+operator_error*hs(x,x).real.sqrt())
    eta=residual+2*max(masses).sqrt()*error+error**2
    norm=hs(C,C).real.sqrt()
    radius=2*eta
    slope=2*max(masses).sqrt()*shift*norm + 2*error*(1+shift)*norm
    contraction=slope+2*((1+shift)*norm)**2*radius
    assert contraction<arb(".5")
    assert strength<4*arb.pi()
    print("zero",index,"fixed coupling",strength.str(4),"all phases; contraction <",contraction.upper().str(7))


if __name__=="__main__":
    for n in (43,195):
        check(n)
