# Publication repair (2026-10-23): require an explicit, hash-pinned theorem input.
"""Replay all three rational-box joint certificates and verify original matrix equations.
No floating point or root finder is used to verify the stored contractions.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse,hashlib,json,re
import sympy as s
import exact_intervals as iv
from certify_intersections import circuit,mfun,overlap,disjoint
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--theory',type=Path,required=True,help='read-only approved FIB_STRICT_REAL_UNEQUAL_DIAGONAL_ACTUAL_M_PAIR.md')
args=parser.parse_args()
raw=args.theory.read_bytes()
assert hashlib.sha256(raw).hexdigest()=='f39a62ba731197f79ba97288739282ba69af7f4602e4b18692beb9f8c309253d', 'theory source hash mismatch'
doc=raw.decode('utf-8')
for name,fun,n,kind in [('existing_M_recheck',mfun,3,'M'),('negative_source_joint_intersection',overlap,6,'overlap'),('disjoint_guard_failure_intersection',disjoint,8,'disjoint')]:
    cert=json.loads((HERE/(name+'.json')).read_text());c=list(map(Q,cert['center']));C=[list(map(Q,row)) for row in cert['preconditioner']];rho=Q(cert['rho']);box,data=iv.contraction(fun,c,C,rho)
    assert data['eta']==cert['eta'] and data['q']==cert['q']
    if kind=='M':
        ww=[circuit(Q(2,5)/box[0],box[0]),circuit(box[1]/box[2],box[2])]
        assert ww[0][1]['b'].lo>0 and ww[1][1]['b'].lo>0 and 0<ww[1][1]['r'].lo<ww[1][1]['r'].hi<iv.GRID
        assert ww[1][1]['r'].lo>iv.I.enc(Q(2,5)).hi
        assert (ww[0][0][2]-ww[1][0][2]).hi<0
        # Locate this root inside the exact published defining box using its rational centre.
        start=doc.index('c_1&=');stop=doc.index('and let $X=',start)
        pairs=re.findall(r'\\frac\{(-?\d+)\}\{(\d+)\}',doc[start:stop]);pub=[Q(int(a),int(b)) for a,b in pairs];assert len(pub)==3
        for i in range(3):assert abs(c[i]-pub[i])+rho<Q(1,10**38)
        # Verify uniqueness in the entire published box without accepting its printed arithmetic on trust.
        pubbox=[iv.I.box(q,Q(1,10**38)) for q in pub]
        iv.NVAR=3;f0=fun([iv.I.enc(q) for q in pub]);dfs=fun(iv.independent(pubbox))
        eta=max(sum((C[i][k]*f0[k] for k in range(3)),iv.I.enc(0)).mag() for i in range(3))
        qq=max(sum((iv.I.enc(int(i==j))-sum((C[i][k]*dfs[k].d[j] for k in range(3)),iv.I.enc(0))).mag() for j in range(3)) for i in range(3))
        assert eta+qq*Q(1,10**38)<Q(1,10**38) and qq<1
    else:
        ww=[circuit(box[2*i],box[2*i+1]) for i in range(n//2)]
        failures=[]
        for i,(_,g) in enumerate(ww):
            if g['r'].hi<0:failures.append([i,'r<0'])
            if g['r'].lo>iv.GRID:failures.append([i,'r>1'])
            if g['b'].hi<0:failures.append([i,'b<0'])
        assert failures==cert['source_guard_failures']
    assert all(w[1]['D'].nz() and w[1]['b'].nz() and w[1]['j'].nz() and w[1]['n'].nz() for w in ww)
    assert all((w[0][0]-24).nz() for w in ww)
    print(name,'RATIONAL_REPLAY_PASS',flush=True)
# Exact literal matrices and one entire gain, independently of response equality tests.
t,b,z=s.symbols('t b z');r=t*b;vals,gdata=circuit(t,b);d=gdata['d'];samples=gdata['samples']
g=sum(samples[i]*s.prod((z-k)/s.Integer(i-k) for k in range(4) if k!=i) for i in range(4));assert s.Poly(g,z).degree()<=3
J=s.Matrix([[0,1],[1,0]]);M=s.Matrix([[0,1],[1,1]]);P=s.Matrix([[r,b],[r*(1-r)/b,1-r]])
Qmat=s.eye(2);targets=[0,24,s.Rational(147,5),24,18,vals[0],vals[1]]
for i in range(7):
    if i:Qmat=[J,M,J,M,M,M][i-1]*Qmat
    assert s.cancel(d+g.subs(z,i)*s.trace(Qmat*P)-targets[i])==0,i
assert s.cancel(d+g.subs(z,7)*s.trace(J*Qmat*P)-vals[2])==0
assert s.cancel(d+g.subs(z,7)*s.trace(M*Qmat*P)-vals[3])==0
assert s.cancel(P.det())==0 and s.trace(P)==1
print('LITERAL_CLOCK_0_THROUGH_7_SYMBOLIC_PASS',flush=True)
print('ONE_CUBIC_OR_LOWER_AND_FIXED_DRIFT_PASS; WHOLE_GAIN_NONZERO_FROM_0_NE_24',flush=True)
