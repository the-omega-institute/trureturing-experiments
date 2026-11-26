"""Minimum real, complex and quaternionic dimensions for three real shell probes."""
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

fixed=run_path(str(Path(__file__).with_name("ch39-rh-fixed-clock-crossings.py")))
quat=run_path(str(Path(__file__).with_name("ch39-rh-heap-finite-shell-metric.py")))


def main():
    ctx.prec,ctx.cap=1536,3
    z,w,_=fixed["shell"]["quad"](21)
    assert len(w)==21 and min(w)>0 and max(w)<arb(1)/4
    # Each entry carries the projection of a current orthonormal coordinate
    # onto the first four axes. Work with the rank-four projection, diagonal 4w.
    eig=[]
    for j in range(21):
        eig.append((arb(int(j<4)),[arb(int(j==k)) for k in range(4)]))
    frame={}
    order=sorted(range(21),key=lambda j:float(w[j].mid()),reverse=True)
    for j in order[:-1]:
        d=4*w[j]
        eig.sort(key=lambda v:float(v[0].mid()),reverse=True)
        k=next(i for i in range(len(eig)-1) if eig[i][0]>d>eig[i+1][0])
        a,x=eig.pop(k)
        b,y=eig.pop(k)
        c=((d-b)/(a-b)).sqrt()
        s=((a-d)/(a-b)).sqrt()
        frame[j]=[(c*px+s*py)/d.sqrt() for px,py in zip(x,y)]
        eig.append((a+b-d,[-s*px+c*py for px,py in zip(x,y)]))
    j=order[-1]
    (d,x),=eig
    assert (d-4*w[j]).contains(0)
    frame[j]=[px/(4*w[j]).sqrt() for px in x]
    tol=arb(2)**(-200)
    assert all(abs(sum(q*q for q in u)-1)<tol for u in frame.values())
    for k in range(4):
        for l in range(4):
            assert abs(sum(w[j]*frame[j][k]*frame[j][l] for j in order)-arb(k==l)/4)<tol
    mul,conj=quat["mul"],quat["conj"]
    basis=((1,0,0,0),(0,1,0,0),(0,0,1,0))
    columns={j:[mul(mul(frame[j],v),conj(frame[j])) for v in basis] for j in order}
    for k in range(3):
        for l in range(3):
            sums=[arb(0) for _ in range(4)]
            for j in order:
                gram=mul(conj(columns[j][k]),columns[j][l])
                assert abs(gram[0]-int(k==l))<tol
                for a in range(4):
                    sums[a]+=w[j]*gram[a]
            assert all(abs(v-int(a==0 and k==l))<tol for a,v in enumerate(sums))
    print("Largest Gaussian weight",max(w).str(12),"verified quaternionic full zero-momentum isometry.")
    print("Three purely real angular inputs, exact 42 readings: real dimension 63, complex 42, quaternionic 21.")


if __name__=="__main__":
    main()
