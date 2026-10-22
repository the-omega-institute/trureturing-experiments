"""Exact rational contraction of joint intersections and explicit source-guard failures.
Certificates here are NOT solutions in the original domain.
"""
from pathlib import Path
from fractions import Fraction as Q
import json,mpmath as mp
import exact_intervals as iv
HERE=Path(__file__).resolve().parent

def circuit(t,b):
    # Reused Support VI Lemma2.2 / XI 6.2, as implemented by negative-j-seed.py.
    r=t*b;j=t+b*(1-t*t);n=j+r;T=1+b;U=1+2*b
    D=5*((b*b+8*b+4)*j*n-2*T*U*(j+n));S=(8*b-1)*(2+5*b)
    C0=-24*D+3*S*j*n;C1=3*S*n;C2=3*(49*b+34)*j*n-54*U*(j+n);C3=3*S*j
    B5=-4*C0+15*C1-20*C2+10*C3;B6=-10*C0+36*C1-45*C2+20*C3;B7=-20*C0+70*C1-84*C2+35*C3
    vals=[(-C0+(U+n)*B5)/D,(-C0+(2*U+n)*B6)/D,(-C0+(T+2*j+4*r)*B7)/D,(-C0+(3*U+2*n)*B7)/D]
    return vals,{'r':r,'b':b,'D':D,'j':j,'n':n,'d':-C0/D,'samples':[C0/D,C1/D,C2/D,C3/D]}
def overlap(z):
    w=[circuit(z[2*i],z[2*i+1])[0] for i in range(3)]
    return [w[0][k]-w[1][k] for k in [0,1,2]]+[w[0][k]-w[2][k] for k in [0,1,3]]
def disjoint(z):
    w=[circuit(z[2*i],z[2*i+1])[0] for i in range(4)]
    return [w[0][k]-w[1][k] for k in [0,1,2]]+[w[2][k]-w[3][k] for k in [0,1,3]]+[w[0][k]-w[2][k] for k in [0,1]]
def mfun(z):
    w0=circuit((mp.mpf(2)/5)/z[0],z[0])[0] if not isinstance(z[0],(iv.I,iv.AD)) else circuit(Q(2,5)/z[0],z[0])[0]
    w1=circuit(z[1]/z[2],z[2])[0]
    return [w0[k]-w1[k] for k in [0,1,3]]
def mpq(q):return mp.mpf(q.numerator)/q.denominator
def certify(name,fun,seed,dimension,kind):
    mp.mp.dps=100
    f=lambda *z:fun(z)
    z=mp.findroot(f,[mp.mpf(str(q)) for q in seed],tol=mp.mpf('1e-90'),maxsteps=100)
    c=[Q(mp.nstr(q,85)) for q in z]
    D=mp.matrix([[mp.diff(lambda u:fun(list(z[:j])+[u]+list(z[j+1:]))[i],z[j]) for j in range(dimension)] for i in range(dimension)])
    Cmp=D**-1;C=[[Q(mp.nstr(Cmp[i,j],85)) for j in range(dimension)] for i in range(dimension)]
    rho=Q(1,10**40) if kind=='existing_M' else Q(1,10**35)
    box,cert=iv.contraction(fun,c,C,rho)
    cert.update({'name':name,'equation':kind,'coordinate_order':'(t0,b0,t1,b1,...)' if kind!='existing_M' else '(b0,r1,b1), r0=2/5','center':[str(a) for a in c],'preconditioner':[[str(a) for a in row] for row in C]})
    vals=[]
    if kind=='existing_M':circuits=[circuit(Q(2,5)/box[0],box[0]),circuit(box[1]/box[2],box[2])]
    else:circuits=[circuit(box[2*i],box[2*i+1]) for i in range(dimension//2)]
    for w,gd in circuits:
        entry={k:v.rounded(15) for k,v in gd.items() if k!='samples'};entry['samples']=[q.rounded(15) for q in gd['samples']];entry['response']=[q.rounded(15) for q in w];vals.append(entry)
        assert gd['D'].nz() and gd['b'].nz() and gd['j'].nz() and gd['n'].nz()
    cert['whole_box_bounds']=vals
    if kind=='existing_M':
        assert all(0<g['r'].lo<g['r'].hi<iv.GRID and g['b'].lo>0 for _,g in circuits[1:]);assert circuits[0][1]['b'].lo>0
        assert circuits[1][1]['r'].lo>iv.I.enc(Q(2,5)).hi
        cert['source_guard_result']='both strict original sources; distinct r; p!=24; terminal J difference strict negative'
        diff=circuits[0][0][2]-circuits[1][0][2];assert diff.hi<0;cert['J_difference']=diff.rounded(15)
    else:
        failures=[]
        for i,(w,g) in enumerate(circuits):
            if g['r'].hi<0:failures.append([i,'r<0'])
            if g['r'].lo>iv.GRID:failures.append([i,'r>1'])
            if g['b'].hi<0:failures.append([i,'b<0'])
        assert failures;cert['source_guard_failures']=failures
        cert['source_guard_result']='inadmissible; this exact joint intersection is not an original-world certificate'
    (HERE/(name+'.json')).write_text(json.dumps(cert,indent=2)+'\n');print(name,'EXACT_CONTRACTION_PASS','eta',float(Q(cert['eta'])),'q',float(Q(cert['q'])),'guard',cert['source_guard_result'],flush=True)
    return cert
if __name__=='__main__':
    certify('existing_M_recheck',mfun,[.040782959048828384,.6375033978214739,.0409169842106061],3,'existing_M')
    certify('negative_source_joint_intersection',overlap,[2.07682067127868,-.00443620843457,2.0925312631715,.32480988594141,48.66650356743,.02008133747385],6,'overlap')
    certify('disjoint_guard_failure_intersection',disjoint,[2.1926609464569244,.13086515643019678,2.4678572244469272,.3131724887686096,3.6638080123154135,.6654736134138814,-1.1704009163556084,.7274310011837438],8,'disjoint')
