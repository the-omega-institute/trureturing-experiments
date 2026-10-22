# Publication repair (2026-10-23): resolve original suppliers within this checkout.
# Source comparison also requires an explicit, hash-pinned theorem input.
"""Verify the compact printed entire-child certificate with rational arithmetic only.
The 21 auxiliary univariate rectangles support one fixed-base certification probe.
No root solver or floating-point operation contributes to the proof.
"""
from pathlib import Path
from fractions import Fraction as Fr
from math import comb
from collections import defaultdict
import argparse,hashlib,json,ast,re
import sympy as sp
import exact_intervals as iv
from certify_intersections import circuit
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--theory',type=Path,required=True,help='read-only approved FIB_STRICT_REAL_UNEQUAL_DIAGONAL_ACTUAL_M_PAIR.md')
args=parser.parse_args()
raw=args.theory.read_bytes()
assert hashlib.sha256(raw).hexdigest()=='f39a62ba731197f79ba97288739282ba69af7f4602e4b18692beb9f8c309253d', 'theory source hash mismatch'
doc=raw.decode('utf-8')
supplier=(Path(__file__).resolve().parents[1]/'fib-self-calibration/fixed-b-and-negative-child/negative-child-exhaustion.py')
for node in ast.parse(supplier.read_text()).body:
 if isinstance(node,(ast.FunctionDef,ast.ClassDef)) and node.name in ['G','asG','base_translate','centre_translate']:exec(compile(ast.Module(body=[node],type_ignores=[]),str(supplier),'exec'))
for node in ast.parse((HERE/'certify_M_base_child.py').read_text()).body:
 if isinstance(node,ast.FunctionDef) and node.name in ['horner_terms','centred_poly']:exec(compile(ast.Module(body=[node],type_ignores=[]),'certify_M_base_child.py','exec'))
data=json.loads((HERE/'existing_M_base_entire_child.json').read_text());pc,vc=map(Fr,data['base_center']);delta=Fr(1,10**24);rho=Fr(1,10**12);sigma=Fr(1,10**24)
co=json.loads((HERE/'inverse_coefficients.json').read_text());x,p,v=sp.symbols('x p v');Ps={k:sp.Poly.from_dict({tuple(e):int(a) for e,a in terms},(x,p,v)) for k,terms in co.items()};R=Ps['R'];A=Ps['alpha'];Z=Ps['zeta'];F=sp.Poly(405*x**3+(975-18*p+5*v)*x**2+(90-36*p+5*v)*x-480,x,p,v)
def rnd(q,n):return Fr((q*10**n+Fr(1,2)).__floor__(),10**n)
ans={'base_center':[str(pc),str(vc)],'base_radius':str(delta),'root_radius':str(rho),'sharp_radius':str(sigma),'eta_bound':'1/5000000000000000000000000','q_bound':'3/1000000000','rows':[],'pole_margins':{}}
margin=[];legal=[]
for name,P,records in [('R',R,data['roots']),('F0',F,data['F0_roots'])]:
    Pbase=base_translate(P);centres=[]
    for record in records:
        c=G(*(rnd(Fr(q),40) for q in record['center']));centres.append(c);tr=centre_translate(Pbase,c);inv=tr[1,0,0].inverse();gamma=G(rnd(inv.a,100),rnd(inv.b,100))
        eta=sum(((gamma*a).norm()*delta**(j+k) for (l,j,k),a in tr.items() if l==0),Fr());rad=rho if c.b==0 else 2*rho
        q=(G(1)-gamma*tr[1,0,0]).norm()+sum(((gamma*a*l).norm()*rad**(l-1)*delta**(j+k) for (l,j,k),a in tr.items() if l>=1 and (l,j,k)!=(1,0,0)),Fr())
        assert eta<Fr(2,10**25) and q<Fr(3,10**9) and eta+q*rho<rho and eta/(1-q)<sigma
        row={'polynomial':name,'id':record['id'],'center':[str(c.a),str(c.b)],'eta':str(eta),'q':str(q)}
        if c.b:assert abs(c.b)>rho;row['filter']='nonreal'
        elif c.a+rho<0:row['filter']='negative x'
        else:
            assert c.a-rho>0;aa=centred_poly(A,c,sigma);zz=centred_poly(Z,c,sigma);assert aa.nz();d=-zz/aa
            xx=iv.I.box(c.a,sigma);pp=iv.I.box(pc,delta);vv=iv.I.box(vc,delta)
            h3=(24-d)/xx;h4=18-d;h5=(pp-d)/(xx+1);h6=(vv-d)/(xx+2)
            H0=20*h3-45*h4+36*h5-10*h6;H1=10*h3-20*h4+15*h5-4*h6;H2=4*h3-6*h4+4*h5-h6;H7=-h3+4*h4-6*h5+4*h6
            B0=Fr(294,5)-2*d-H2;chi=d*d+B0*B0;assert chi.lo>0 and H1.nz();U=(-d*H0+B0*H2)/chi;assert U.nz();b=(U-1)/2;j=U*(24-d)/H1;r=xx*U-j;a=d+((1+1/U)/2+2*xx+2*r/U)*H7
            low=lambda a:Fr(min(abs(a.lo),abs(a.hi)),iv.GRID) if a.nz() else Fr(0)
            margin.append([low(aa),Fr(chi.lo,iv.GRID),low(H1),low(U)])
            row['bounds']={k:a.rounded(8) for k,a in [('x',xx),('b',b),('r',r),('a',a)]}
            if b.hi<0:row['filter']='b<0'
            elif r.hi<0:row['filter']='r<0'
            elif r.lo>iv.GRID:row['filter']='r>1'
            else:assert b.lo>0 and 0<r.lo<r.hi<iv.GRID;row['filter']='strict possible';legal.append((name,record['id'],r,b,a))
        ans['rows'].append(row)
    for i,c in enumerate(centres):
        for z in centres[:i]:assert abs(c.a-z.a)>2*rho or abs(c.b-z.b)>2*rho
    assert len(centres)==P.degree(x)
lead=sp.Poly.from_dict({(0,e[1],e[2]):a for e,a in R.terms() if e[0]==18},(x,p,v));leadbox=centred_poly(lead,G(),Fr(0));assert leadbox.lo>0;ans['leading_R']=leadbox.rounded(0)
legal.sort(key=lambda row:row[2].lo);assert len(legal)==2
Jdiff=legal[0][4]-legal[1][4];assert Jdiff.hi<0;ans['J_difference']=Jdiff.rounded(5);ans['legal_root_ids']=[[name,i] for name,i,*_ in legal]
ans['pole_margins']={k:str(min(row[i] for row in margin)) for i,k in enumerate(['abs_alpha','chi','abs_H1','abs_U'])}
assert all(Fr(ans['pole_margins'][k])>bound for k,bound in [('abs_alpha',Fr(10**14)),('chi',Fr(40)),('abs_H1',Fr(1,25)),('abs_U',Fr(1,10))])
assert list(map(Fr,ans['leading_R']))==[Fr(1873044541936855184250),Fr(1873044541936855184251)]
assert Fr(-3076889,10**5)<=Fr(Jdiff.lo,iv.GRID)<=Fr(Jdiff.hi,iv.GRID)<=Fr(-3076888,10**5)
# The whole originally printed M certificate box maps into this base rectangle.
start=doc.index('c_1&=');stop=doc.index('and let $X=',start);pub=[Fr(int(a),int(b)) for a,b in re.findall(r'\\frac\{(-?\d+)\}\{(\d+)\}',doc[start:stop])]
bb=iv.I.box(pub[0],Fr(1,10**38));w,g=circuit(Fr(2,5)/bb,bb)
assert (w[0]-pc).mag()<Fr(1,10**30) and (w[1]-vc).mag()<Fr(1,10**30)
ans['published_box_base_deviation_bound']='1/1000000000000000000000000000000'
# Compare the exact printed values with the certificate, including outward source filters.
section=doc[doc.index('## 9. The entire child'):]
assert Fr(re.search(r'p_c&=([0-9.]+)',section).group(1))==pc
assert Fr(re.search(r'v_c&=([0-9.]+)',section).group(1))==vc
printed_centres=[];printed_inverse={}
for line in section.splitlines():
    cells=[q.strip() for q in line.strip().strip('|').split('|')]
    if cells[0] not in ['$R$','$F0$']:continue
    name=cells[0].strip('$');i=int(cells[1])
    if len(cells)==4:
        printed_centres.append((name,i,[Fr(q.strip('$')) for q in cells[2:]]))
    elif len(cells)==5:
        printed_inverse[name,i]=[list(map(Fr,re.findall(r'-?[0-9]+\.[0-9]+',q))) for q in cells[2:4]]
assert printed_centres==[(row['polynomial'],row['id'],list(map(Fr,row['center']))) for row in ans['rows']]
for row in ans['rows']:
    if 'bounds' not in row:continue
    printed=printed_inverse.pop((row['polynomial'],row['id']))
    for key,limits in zip(['b','r'],printed):
        assert len(limits)==2 and limits==list(map(Fr,row['bounds'][key]))
assert not printed_inverse
for lo,i,hi in re.findall(r'([0-9.]+)&\\le a_\{(12|13)\}\\le([0-9.]+)',section):
    row=next(row for row in ans['rows'] if row['polynomial']=='R' and row['id']==int(i))
    assert [Fr(lo),Fr(hi)]==list(map(Fr,row['bounds']['a']))
(HERE/'paper_certificate.json').write_text(json.dumps(ans,indent=2)+'\n')
print('PRINTED_CERTIFICATE_RATIONAL_PASS; 18 R ROOTS + 3 F0 ROOTS; 2 POSSIBLE WHOLE WORLDS; J SEPARATED',flush=True)
print('PRINTED_SOURCE_CENTRES_FILTERS_AND_BOUNDS_MATCH',flush=True)
print('POLE_LOWER_BOUNDS',ans['pole_margins'],flush=True)
print('J_DIFFERENCE',ans['J_difference'],flush=True)
