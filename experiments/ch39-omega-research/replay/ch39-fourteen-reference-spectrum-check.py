"""Cross-check the 14-reference exposure with direct symmetric matrix diagonalization."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb,ctx

directory=Path(__file__).parent
uv=run_path(str(directory/"ch39-kinetic-rh-uv-running.py"))


def main():
    ctx.prec=384
    base=uv["shell"]["base"]
    a,b,c=[base["gram"](i,j) for i,j in ((1,0),(1,2),(0,2))]
    p=(1+2*a*b*c-a*a-b*b-c*c)/(1-c*c)
    z,w=uv["finite_model"](14)
    poles=np.array([float(sum(z)/(14*v)) for v in z])
    residues=np.array([float((1-p)*v) for v in w])
    d=8*np.pi**2/np.square([.467453,.63811,1.0583])
    K=np.dot(d,[5,-3,-2])
    F=109/15*(d[1]-d[2])-23/6*(d[0]-d[1])
    Q=1000*np.exp(K/44)
    modes=[(6,7),(4,6),(3+np.sqrt(5),3),(3-np.sqrt(5),3),(3+np.sqrt(3),2),(3-np.sqrt(3),2)]
    for mf,cycles,expected in [(6.52316395e13,True,119),(5.17571685e13,False,95)]:
        values=list(np.repeat(poles,25)) if cycles else []
        for lam,n in modes:
            M=np.diag(np.r_[lam+(125/mf)**2,poles])
            M[0,1:]=M[1:,0]=np.sqrt(lam*poles*residues)
            values.extend(np.repeat(np.linalg.eigvalsh(M),n))
        exposed=np.maximum(0,np.log(Q/(mf*np.sqrt(values))))
        assert sum(exposed>0)==expected and abs(22*sum(exposed)-15*F)<.001
        print("Direct numerical eigensystem:",len(values),"scalar masses,",expected,"open at Q.")


if __name__=="__main__":
    main()
