"""Allow coherent correlations between the twelve clock-edge perturbations."""
from pathlib import Path
import runpy
from flint import acb,arb,acb_mat

edge=runpy.run_path(str(Path(__file__).with_name("ch39-clock-edge-source.py")))
s=edge["s"]
src=edge["src"]


def inertia(matrix,value):
    a=acb_mat(matrix)
    negatives=0
    for i in range(a.nrows()):
        a[i,i]-=value
    for j in range(a.nrows()):
        pivot=a[j,j].real
        assert not pivot.contains(0)
        negatives+=int(pivot<0)
        for i in range(j+1,a.nrows()):
            for k in range(j+1,a.nrows()):
                a[i,k]-=a[i,j]*a[j,k]/pivot
    return negatives


if __name__=="__main__":
    p,q=edge["p"],edge["q"]
    m=s["matrix"](0,{724:acb(0,p).exp(),804:acb(0,q).exp()})
    H=acb_mat(9,12)
    for col,(i,j,d) in enumerate(edge["edges"]):
        H[i,col]=1
        H[j,col]=m[i,j].conjugate()
    assert inertia(H*H.conjugate().transpose(),arb(".34"))==0
    G=src["block"](m.inv(),range(9),range(9))
    assert inertia(G,arb(".0992547"))==0 and inertia(G,arb(".0992548"))==1
    assert inertia(G,arb(".11092"))==1
    print("physical trace-budget minimum eigenvalue in (.0992547,.0992548), gap > .01166")
