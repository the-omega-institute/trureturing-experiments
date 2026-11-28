# Publication repair (2026-10-23): resolve original suppliers within this checkout.
"""Exact full-child inventory at the existing M-pair base.
Reuses the published normalized carrier and the private negative-child certificate's
Gaussian-rational translation/contraction machinery. No search miss is an exclusion.
"""
from pathlib import Path
from fractions import Fraction as Fr
from math import comb
from collections import defaultdict
import json,ast,mpmath as mp,sympy as sp
import exact_intervals as iv
from certify_intersections import circuit
HERE=Path(__file__).resolve().parent
oldpath=(Path(__file__).resolve().parents[1]/'fib-self-calibration/fixed-b-and-negative-child/negative-child-exhaustion.py')
code=ast.parse(oldpath.read_text());reused=['G','asG','base_translate','centre_translate']
for node in code.body:
    if isinstance(node,(ast.ClassDef,ast.FunctionDef)) and node.name in reused:exec(compile(ast.Module(body=[node],type_ignores=[]),str(oldpath),'exec'))
coeff=json.loads((HERE/'inverse_coefficients.json').read_text());x,p,v=sp.symbols('x p v')
Ps={k:sp.Poly.from_dict({tuple(e):int(c) for e,c in terms},(x,p,v)) for k,terms in coeff.items()}
oldcert=json.loads((HERE/'existing_M_recheck.json').read_text());cs=[Fr(q) for q in oldcert['center']];rhoM=Fr(oldcert['rho'])
# Bound the exact existing-M base from its source certificate, before choosing this base rectangle.
boxM=[iv.I.box(c,rhoM) for c in cs]
w0,gd=circuit(Fr(2,5)/boxM[0],boxM[0]);pc=Fr(w0[0].lo+w0[0].hi,2*iv.GRID);vc=Fr(w0[1].lo+w0[1].hi,2*iv.GRID)
pc=Fr(round(pc*10**40),10**40);vc=Fr(round(vc*10**40),10**40)
delta=Fr(1,10**24);rho=Fr(1,10**12)
assert w0[0].lo>iv.I.enc(pc-delta).lo and w0[0].hi<iv.I.enc(pc+delta).hi
assert w0[1].lo>iv.I.enc(vc-delta).lo and w0[1].hi<iv.I.enc(vc+delta).hi
R=Ps['R'];Rbase=base_translate(R)
A=Ps['alpha'];Z=Ps['zeta'];Abase=base_translate(A);Zbase=base_translate(Z)
spc=sp.Rational(pc.numerator,pc.denominator);svc=sp.Rational(vc.numerator,vc.denominator)
print('NROOTS_START',flush=True)
roots=sp.nroots(R.as_expr().subs({p:spc,v:svc}),n=75,maxsteps=500)
centres=[G(Fr(str(sp.re(z).evalf(65))),Fr(str(sp.im(z).evalf(65)))) for z in roots]
print('NROOTS_DONE',len(centres),flush=True)
def rnd(q,n=100):return Fr((q*10**n+Fr(1,2)).__floor__(),10**n)
def intervalpoly(poly,c,sigma):
    tr=centre_translate(poly,c);mid=tr.get((0,0,0),G());assert mid.b==0
    rad=sum((a.norm()*sigma**i*delta**(j+k) for (i,j,k),a in tr.items() if (i,j,k)!=(0,0,0)),Fr())
    return iv.I.box(mid.a,rad)
def horner_terms(terms,arguments):
    if not arguments:return iv.I.enc(sum(terms.values()))
    groups={}
    for e,c in terms.items():groups.setdefault(e[0],{})[e[1:]]=c
    z=iv.I.enc(0)
    for i in range(max(groups),-1,-1):z=z*arguments[0]+(horner_terms(groups[i],arguments[1:]) if i in groups else 0)
    return z

def centred_poly(P,c,sigma):
    assert c.b==0
    values=[c.a,pc,vc];radii=[sigma,delta,delta]
    pt=[iv.I.enc(z) for z in values];box=[iv.I.box(z,r) for z,r in zip(values,radii)]
    terms={tuple(e):int(a) for e,a in P.terms()}
    ans=horner_terms(terms,pt)
    for k in range(3):
        derivative={tuple(e[j]-int(j==k) for j in range(3)):a*e[k] for e,a in terms.items() if e[k]}
        if derivative:ans=ans+horner_terms(derivative,box)*iv.I.box(Fr(0),radii[k])
    return ans

cert={'status':'exact finite-base certificate, not global','source_base_certificate':'existing_M_recheck.json','carrier_supplier':str(oldpath),'carrier_normalization':'IX 4.7 / X12.12, SHA256 ff6b65b00b5d251082949e2318c2a53fb3f583a86c94b3d22c3e15c9ebe038e2','base_center':[str(pc),str(vc)],'base_radius':str(delta),'root_radius':str(rho),'roots':[]}
# Leading coefficient has fixed positive degree 18 throughout the base rectangle.
leadterms={(0,e[1],e[2]):c for e,c in R.terms() if e[0]==18}
lead=sp.Poly.from_dict(leadterms,(x,p,v));leadbox=intervalpoly(base_translate(lead),G(),Fr(0));assert leadbox.nz();cert['leading_coefficient']=leadbox.rounded(12)
positive=[];legal=[]
for i,c in enumerate(centres):
    tr=centre_translate(Rbase,c);inv=tr[1,0,0].inverse();gamma=G(rnd(inv.a),rnd(inv.b));eta=sum(((gamma*a).norm()*delta**(j+k) for (l,j,k),a in tr.items() if l==0),Fr());rad=rho if c.b==0 else 2*rho
    q=(G(1)-gamma*tr[1,0,0]).norm()+sum(((gamma*a*l).norm()*rad**(l-1)*delta**(j+k) for (l,j,k),a in tr.items() if l>=1 and (l,j,k)!=(1,0,0)),Fr())
    assert q<Fr(1,10) and eta+q*rho<rho,(i,float(eta),float(q))
    # A bound on all possible root displacement, sharpened from contraction.
    sigma=2*eta;assert eta/(1-q)<sigma<rho
    row={'id':i,'center':[str(c.a),str(c.b)],'gamma':[str(gamma.a),str(gamma.b)],'eta':str(eta),'q':str(q),'root_displacement_radius':str(sigma)}
    if c.b!=0:
        assert abs(c.b)>rho;row['filter']='nonreal x'
    elif c.a+rho<0:row['filter']='negative x'
    else:
        assert c.a-rho>0;positive.append(i)
        aa=centred_poly(A,c,sigma);zz=centred_poly(Z,c,sigma);assert aa.nz(),i
        d=-zz/aa;xx=iv.I.box(c.a,sigma);pp=iv.I.box(pc,delta);vv=iv.I.box(vc,delta)
        h3=(24-d)/xx;h4=18-d;h5=(pp-d)/(xx+1);h6=(vv-d)/(xx+2)
        H0=20*h3-45*h4+36*h5-10*h6;H1=10*h3-20*h4+15*h5-4*h6;H2=4*h3-6*h4+4*h5-h6;H7=-h3+4*h4-6*h5+4*h6
        B0=Fr(294,5)-2*d-H2;chi=d*d+B0*B0;assert chi.lo>0 and H1.nz()
        U=(-d*H0+B0*H2)/chi;b=(U-1)/2;j=U*(24-d)/H1;r=xx*U-j
        m=d+(3+2*xx)*H7;a=d+((1+1/U)/2+2*xx+2*r/U)*H7
        row['inverse']={key:val.rounded(12) for key,val in [('x',xx),('d',d),('U',U),('b',b),('r',r),('a',a),('m',m)]}
        if b.hi<0:row['filter']='b<0'
        elif r.hi<0:row['filter']='r<0'
        elif r.lo>iv.GRID:row['filter']='r>1'
        else:
            assert b.lo>0 and 0<r.lo<r.hi<iv.GRID,(i,row['inverse']);row['filter']='strict possible original companion';legal.append((i,r,b,a,m))
    cert['roots'].append(row);print('ROOT',i,row['filter'],row.get('inverse',{}),flush=True)
# Disjoint complex rectangles: degree 18 then proves exhaustive complex root coverage.
for i,c in enumerate(centres):
    for e in centres[:i]:assert abs(c.a-e.a)>2*rho or abs(c.b-e.b)>2*rho
# The full resultant is 256 D^8 F0^2 R. Retain F0 as well, including zero drift.
F0=sp.Poly(405*x**3+(975-18*p+5*v)*x**2+(90-36*p+5*v)*x-480,x,p,v)
Fbase=base_translate(F0);froots=sp.nroots(F0.as_expr().subs({p:spc,v:svc}),n=75,maxsteps=500)
cert['F0_roots']=[]
for ii,root in enumerate(froots):
    c=G(Fr(str(sp.re(root).evalf(65))),Fr(str(sp.im(root).evalf(65))))
    tr=centre_translate(Fbase,c);inv=tr[1,0,0].inverse();gamma=G(rnd(inv.a),rnd(inv.b));eta=sum(((gamma*a).norm()*delta**(j+k) for (l,j,k),a in tr.items() if l==0),Fr());rad=rho if c.b==0 else 2*rho
    q=(G(1)-gamma*tr[1,0,0]).norm()+sum(((gamma*a*l).norm()*rad**(l-1)*delta**(j+k) for (l,j,k),a in tr.items() if l>=1 and (l,j,k)!=(1,0,0)),Fr())
    assert q<Fr(1,10) and eta+q*rho<rho
    row={'id':ii,'center':[str(c.a),str(c.b)],'gamma':[str(gamma.a),str(gamma.b)],'eta':str(eta),'q':str(q)}
    if c.b!=0:assert abs(c.b)>rho;row['filter']='nonreal x'
    elif c.a+rho<0:row['filter']='negative x'
    else:
        sigma=2*eta;aa=centred_poly(A,c,sigma);zz=centred_poly(Z,c,sigma);assert aa.nz();d=-zz/aa
        xx=iv.I.box(c.a,sigma);pp=iv.I.box(pc,delta);vv=iv.I.box(vc,delta)
        h3=(24-d)/xx;h4=18-d;h5=(pp-d)/(xx+1);h6=(vv-d)/(xx+2)
        H0=20*h3-45*h4+36*h5-10*h6;H1=10*h3-20*h4+15*h5-4*h6;H2=4*h3-6*h4+4*h5-h6;B0=Fr(294,5)-2*d-H2
        chi=d*d+B0*B0;assert chi.lo>0 and H1.nz();U=(-d*H0+B0*H2)/chi;b=(U-1)/2;j=U*(24-d)/H1;r=xx*U-j
        rank=b*j-b*b-r*(1-r)
        row['inverse']={key:val.rounded(10) for key,val in [('x',xx),('d',d),('U',U),('b',b),('r',r),('rank_residual',rank)]}
        if b.hi<0:row['filter']='b<0'
        elif r.hi<0:row['filter']='r<0'
        elif r.lo>iv.GRID:row['filter']='r>1'
        else:assert rank.nz();row['filter']='rank residual !=0'
    cert['F0_roots'].append(row);print('F0_ROOT',ii,row['filter'],row.get('inverse',{}),flush=True)
assert len(cert['F0_roots'])==3
assert len(centres)==18 and len(legal)==2
legal.sort(key=lambda q:q[1].lo)
# Match the two possible original companions to the already existing exact M endpoints by r.
assert legal[0][1].lo<iv.I.enc(Fr(2,5)).lo<legal[0][1].hi
assert legal[1][1].lo<boxM[1].lo<boxM[1].hi<legal[1][1].hi
Jdiff=legal[0][3]-legal[1][3];assert Jdiff.nz();cert['actual_J_difference']=Jdiff.rounded(10);cert['legal_root_ids']=[z[0] for z in legal]
# Singular regular-chart sources: parent equation forces b=1/8 and Lambda=0.
r0=sp.Symbol('r');Phi=438272*r0**4-931328*r0**3+540800*r0**2-50728*r0-853
sing=[]
for (lo,hi),multiplicity in sp.Poly(Phi,r0).intervals(eps=sp.Rational(1,10**22)):
    if hi<=0 or lo>=1:continue
    rr=iv.I.box(Fr(int((lo+hi).p),int((lo+hi).q))/2,Fr(int((hi-lo).p),int((hi-lo).q))/2);bb=Fr(1,8);jj=bb+rr*(1-rr)/bb;nn=jj+rr;T=1+bb;U=1+2*bb
    p0=24+(U+nn)*(96-Fr(108,1)/T);Ap=1+(U+nn)*(4-15/jj+20/T-10/nn)
    v0=24+(2*U+nn)*(240-Fr(243,1)/T);Av=1+(2*U+nn)*(10-36/jj+45/T-20/nn)
    sep=(iv.I.box(pc,delta)-p0)*Av-(iv.I.box(vc,delta)-v0)*Ap
    assert sep.nz();sing.append({'r':[str(lo),str(hi)],'response_line_residual':sep.rounded(10)})
assert len(sing)==2,len(sing)
cert['singular_sources']=sing
cert['conclusion']='The entire original child at the existing exact M-pair base consists of exactly those two strict worlds. Their actual J answers differ, so that one base is not a solution of the independent-pair target. Both p and v remain free in the original unresolved problem.'
(HERE/'existing_M_base_entire_child.json').write_text(json.dumps(cert,indent=2)+'\n')
print('FULL_CHILD_EXACT_PASS',len(legal),'legal roots','Jdiff',cert['actual_J_difference'],'singular lines',len(sing),flush=True)
