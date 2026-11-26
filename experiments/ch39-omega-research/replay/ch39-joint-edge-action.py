"""Joint zero-height and single-edge selection at any fixed total stiffness."""
from pathlib import Path
import runpy
from flint import acb
import numpy as np

joint=runpy.run_path(str(Path(__file__).with_name("ch39-joint-clock-action.py")))
fit=joint["fit"]
block,J=fit["block"],fit["J"]
adj=joint["adj"]


def records(p,q):
    m,orbits=fit["model"](p,q)
    det=m.det().real/joint["reference"]
    result=[]
    for orbit in orbits:
        h=orbit[0]
        resistance=(adj(h)*m.solve(h))[0,0].real
        scores=[det*((1-x)+4*x*resistance) for x in joint["ray"]["nodes"]]

        result.append(joint["controls"](scores)+[acb(0)]*8)
    return result


def face_filter(orbit):
    edge=joint["edge"]
    data=edge["coefficients"](orbit)
    m,orbits=fit["model"](0,0)
    hj=block(orbits[orbit][0],J,(0,))
    g=(adj(hj)*block(m,J,J).solve(hj))[0,0].real
    # Convert the affine numerator to x=S/(4+S), preserving each phase.
    base=data[0].transpose(2,0,1)
    correction=data[1].transpose(2,0,1)
    l=edge["times"](base,edge["rounded"](acb(4*g-1)))
    r=edge["times"](correction,edge["rounded"](acb(4*g)))
    data[1]=np.stack([a+b for a,b in zip(l,r)],axis=-1)
    def reject(level,idx,better):
        result=np.zeros(len(idx),dtype=bool)
        if better.all():
            return result
        a,b,r=edge["field"](data,level,idx[~better])
        excluded=np.zeros(len(a),dtype=bool)
        for f in edge["top_faces"]:
            loop=(edge["Q"],0,0)
            for u,w in zip(f,f[1:]+f[:1]):
                j=9+edge["top_edges"].index((u,w))
                loop=edge["times"](loop,(a[:,j],b[:,j],r[:,j]))
            re,im,err=loop
            assert max(int(np.max(abs(v))) for v in loop)<2**57
            excluded|=(im+err<0)|(re-err>0)|(9*(im-err)+re-err>0)
        result[~better]=excluded
        return result
    return reject


if __name__=="__main__":
    branches,witnesses=joint["inputs"](records)
    for i,branch in enumerate(branches):
        print("edge orbit, retained harmonics, boxes, level",i,len(branch),
              *joint["cover"](branch,witnesses,face_filter(i)),flush=True)
