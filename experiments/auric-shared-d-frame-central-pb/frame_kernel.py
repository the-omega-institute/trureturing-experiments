"""Exact modular lift for shared-D frame diagonality; no PB numerical sampling."""
import hashlib, json, sys, math
import numpy as np
from flint import nmod_mat
P=1009; II=469
assert all(P % k for k in range(2,math.isqrt(P)+1))
assert II*II%P==P-1
inv2=pow(2,-1,P); inv4=pow(4,-1,P)
eye=np.eye(4,dtype=np.int64)
pairs=[(i,j) for i in range(4) for j in range(i+1,4)]
B=[]
for i in range(3):
    a=np.zeros((4,4),dtype=np.int64); a[i,i]=1; a[3,3]=-1; B.append(a%P)
for i,j in pairs:
    for z in (1,II):
        a=np.zeros((4,4),dtype=np.int64); a[i,j]=z; a[j,i]=z if z==1 else -z; B.append(a%P)
B=np.array(B)
def coord(a):
    tr=int(np.trace(a))%P
    ans=[(int(a[i,i])-tr*inv4)%P for i in range(3)]
    for i,j in pairs:
        ans += [(int(a[i,j])+int(a[j,i]))*inv2%P,
                (int(a[i,j])-int(a[j,i]))*pow(2*II,-1,P)%P]
    return np.array(ans,dtype=np.int64)
for j,a in enumerate(B):
    assert np.array_equal(coord(a),np.eye(15,dtype=np.int64)[j])
T=[]
for i in (0,1):
    for z in (1,II):
        a=np.zeros((4,4),dtype=np.int64); a[i,3]=z; a[3,i]=z if z==1 else -z; T.append(a%P)
c=np.diag([1,2,3,0]).astype(np.int64)
wedges=[(i,j) for i in range(15) for j in range(i+1,15)]
seed=1; denominators=[]; vectors=[]
def frame():
    global seed
    u=eye.copy(); us=eye.copy()
    for h in range(4):
        v=[]
        for k in range(8):
            seed=(1664525*seed+1013904223)%2**32
            v.append(seed%5-2)
        vr=np.array(v[:4],dtype=np.int64); vi=np.array(v[4:],dtype=np.int64)
        d=int(vr@vr+vi@vi); assert d>0 and d%P
        denominators.append(d); vectors.append(v)
        w=(vr+II*vi)%P; ws=(vr-II*vi)%P
        hmat=(eye-2*pow(d,-1,P)*np.outer(w,ws))%P
        assert np.array_equal(hmat@hmat%P,eye)
        u=u@hmat%P; us=hmat@us%P
    assert np.array_equal(us@u%P,eye)
    return u,us
N=int(sys.argv[1]) if len(sys.argv)>1 else 48
rows=[]; inputs=[]
for f in range(N):
    u,us=frame()
    a=coord(u@c%P@us%P)
    out=np.array([coord(us@b%P@u%P)[3:] for b in B],dtype=np.int64).T
    assert out.shape==(12,15)
    for t in T:
        b=coord(u@t%P@us%P)
        w=np.array([(a[i]*b[j]-a[j]*b[i])%P for i,j in wedges],dtype=np.int64)
        inputs.append(w)
        for d in out:
            rows.append(np.outer(w,d).ravel()%P)
mat=np.array(rows,dtype=np.int64)
md=hashlib.sha256(mat.astype('<u2').tobytes()).hexdigest()
print(json.dumps({'frames':N,'shape':list(mat.shape),'prime':P,'root_i':II,
 'matrix_sha256_u16le':md,'denominator_minmax':[min(denominators),max(denominators)]}),flush=True)
mm=nmod_mat(mat.tolist(),P)
r=mm.rank()
rr,rt=mm.transpose().rref()
pivots=[]
for k in range(rt):
    pivots.append(next(j for j in range(rr.ncols()) if int(rr[k,j])))
minor=nmod_mat(mat[pivots,:].tolist(),P)
det=int(minor.det())
wrank=nmod_mat(np.array(inputs).tolist(),P).rank()
print(json.dumps({'rank':r,'nullity':1575-r,'input_wedge_rank':wrank,
 'minor_determinant_mod_p':det,
 'pivot_row_indices_sha256_u16le':hashlib.sha256(np.array(pivots,dtype='<u2').tobytes()).hexdigest(),
 'denominator_sequence_sha256_u16le':hashlib.sha256(np.array(denominators,dtype='<u2').tobytes()).hexdigest()}),flush=True)
assert r==1575 and rt==1575 and det!=0 and wrank==105
