"""Construct a local isospectral rotation with two quark-inert heavy doublets."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

flow = run_path(str(Path(__file__).with_name("ch39-shell-gw-first-heavy-flow.py")))
fit = flow["fit"]
common = fit["common"]
adj, hs, rows, count = (fit[key] for key in ("adj","hs","rows","count"))
array = common["flavor"]["array"]
g = arb("7.9")


def exact_array(z):
    return acb_mat([[acb(float(v.real),float(v.imag)) for v in row] for row in z])


def simple_mode(M,index):
    v,V = np.linalg.eigh(array(M))
    e = arb(float(v[index]))
    gap = arb(float(min(abs(v[index]-v[j]) for j in range(len(v)) if j!=index)))/2
    assert [count(rows(M),e+s*gap) for s in (-1,1)] == [index,index+1]
    t = exact_array(V[:,index:index+1])
    t /= hs(t,t).real.sqrt()
    residual = M*t-e*t
    error = arb(0,(2*hs(residual,residual).real.sqrt()/gap).upper())
    return acb_mat([[z+acb(error,error)] for z in t.entries()])


def fingerprints(index):
    model = flow["full"]["model"]
    if index == 43:
        S,F,_,_ = model["m"]["transported_spectra"]()
    else:
        S,F,_ = model["scalar_models"]()
    assert count(rows(S), arb(".1")) == 1 and count(rows(F), arb(".1")) == 2
    l = simple_mode(S,0)
    h,k = (simple_mode(F,j) for j in (0,1))
    assert hs(l,S*l).real < hs(h,F*h).real < hs(k,F*k).real
    target = []
    for j in range(6):
        z = acb_mat([[l[j,0],0,0]]+[[0,h[8*j+i,0],k[8*j+i,0]] for i in range(8)])
        target.append(adj(z)*z)
    return target


def axes():
    result = []
    for i in range(3):
        A = acb_mat(3,3)
        A[i,i] = 1
        result.append(A)
        for j in range(i):
            for v in (acb(1),acb(0,1)):
                A = acb_mat(3,3)
                A[i,j] = v/arb(2).sqrt()
                A[j,i] = A[i,j].conjugate()
                result.append(A)
    return result


def trial(P,T,K,target):
    a,t,k = map(array,(P,T,K))
    pinv = a.conj().T@np.linalg.inv(a@a.conj().T)
    rng = np.random.default_rng(88)
    Z = pinv@t+k@(rng.normal(size=(54,3))+1j*rng.normal(size=(54,3)))*.055
    H = np.array([array(b) for b in axes()])
    G = np.array([array(b) for b in target])
    for _ in range(60):
        X = Z.reshape(6,9,3)
        delta = X.conj().transpose(0,2,1)@X-G
        err = np.einsum("kij,bji->bk",H,delta).real.ravel()
        if np.linalg.norm(err)<1e-12:
            break
        gradients = np.zeros((54,54,3),complex)
        for j in range(6):
            gradients[9*j:9*j+9,9*j:9*j+9,:] = 2*(X[j]@H)
        jac = (k@gradients).reshape(54,-1)
        step = (jac.T@np.linalg.solve((jac.conj()@jac.T).real,err)).reshape(54,3)
        for i in range(15):
            zz = (Z-(.5**i)*step).reshape(6,9,3)
            if np.linalg.norm(zz.conj().transpose(0,2,1)@zz-G)<np.linalg.norm(err):
                Z = zz.reshape(54,3)
                break
        else:
            break
    return exact_array(Z)


def certify(index):
    heavy,up,down,lepton,F,D,*_ = flow["at_threshold"](index)
    P = acb_mat(rows(F)+rows(D))
    T = adj(P)*(P*adj(P)).inv()
    K = common["old"]["m"]["eye"](54)-T*P
    U,A,B = [fit["orbit"]["exact_frame"](np.array(z).reshape(3,3)) for z in common["old"]["units"]]
    target = acb_mat([[z/g,0,0] for z in (U*up*A).entries()+(U*down*B).conjugate().entries()])
    grams = fingerprints(index)
    trial_z = trial(P,target,K,grams)
    Z = trial_z+T*(target-P*trial_z)
    basis = axes()
    residual2 = arb(0)
    gradients = []
    for j in range(6):
        X = acb_mat(rows(Z)[9*j:9*j+9])
        residual = adj(X)*X-grams[j]
        residual2 += sum(abs(z).abs_upper()**2 for z in residual.entries())
        for A in basis:
            h = (2*X*A).table()
            direction = K*acb_mat([h[i%9] if i//9==j else [0,0,0] for i in range(54)])
            gradients.append(direction.entries())
    V = acb_mat(gradients)
    R = V*adj(V)
    real = acb_mat([[z.real for z in row] for row in R.table()])
    inverse_bound = arb(40)
    assert count(rows(real),1/inverse_bound**2) == 0
    radius = 2*residual2.sqrt().abs_upper()
    assert 2*inverse_bound**2*radius < arb(".5")
    # Both heavy lines are pointwise orthogonal to the light line. Flip both
    # signs at site zero to keep every local Gram and the fitted light column.
    F0 = acb_mat([r[:9] for r in rows(F)])
    singular = arb(".0374")
    assert count(rows(adj(F0)*F0),singular**2)==0
    active = 2*g*singular*grams[0][1,1].real.sqrt()
    assert active > arb(".241")
    print("zero",index,"two inert heavy quark vertices; alternative first-heavy up norm above",
          active.lower().str(9),"contraction below",(2*inverse_bound**2*radius).upper().str(6))
    return heavy, [up,down,lepton]



def continue_inert(index):
    heavy, data = certify(index)
    times = [(m/1000).log()/flow["loop"] for m in heavy[:3]]
    gauge = flow["gauge"](times[0],(flow["start"],flow["beta"]))
    for k in (1,2):
        slopes = [b+k*d for b,d in zip(flow["beta"],(arb(1)/10,arb(1)/6,0))]
        interval = times[k]-times[k-1]
        assert interval > 0
        params = (gauge,slopes)
        data = flow["evolve"](interval,data,params)
        gauge = flow["gauge"](interval,params)
    norms = [hs(M,M).real.sqrt() for M in data]
    assert arb(".50") < norms[0] < arb(".52")
    print("third heavy GeV",heavy[2].str(10),"light Yukawa norms immediately below",
          *[v.str(10) for v in norms])


if __name__ == "__main__":
    for n in (43,195):
        continue_inert(n)
