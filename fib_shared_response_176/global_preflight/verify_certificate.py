import json,time,sys
from pathlib import Path
import sympy as s
from flint import fmpq_mpoly_ctx,fmpq
P=Path(__file__).parent;c=fmpq_mpoly_ctx.get(('u','b','e'),'lex');fu,fb,fe=c.gens();t=time.perf_counter();results={}
def rd(n):
 d=json.loads((P/(n+'.json')).read_text());return c.from_dict({tuple(k):fmpq(v) for k,v in d['terms']})
def conv(ex):return c.from_dict({m:fmpq(int(co.p),int(co.q)) for m,co in s.Poly(s.expand(ex),u,b,e,domain=s.QQ).terms()})
u,b,e,n,a=s.symbols('u b e n a');h=u*(u+1);ell=u*u+3*u-4;al=3*h-10;be=8*u*u+18*u-30
phi=al*a+be*n-23*h;delta=n*h+a*n*ell-4*a*h
T=(s.Rational(3,7)*n-a)*phi+s.Rational(175,29)*delta
Q=a*a-a*n-2*(u+1)*a-n*n+2*(2*u+1)*n-u*u+u+1
A=(10-3*b)*u*u+(8-3*b)*u+10*b-32;B=(10-8*b)*u*u+(19-18*b)*u+30*b-36;C=(23*b-29)*h
ps6=(n+2*a)*(10*u*u+19*u-36)-10*a*ell-29*h
ps7=4*(2*n+3*a)*(20*u*u+34*u-70)-105*a*ell-381*h
results['child_line_identity']=s.expand(ps6-b*phi-(A*a+B*n+C))==0
xsub={a:-B*n-C,n:A*n} # simultaneous homogeneous chart substitution uses homogenization explicitly
QA=s.Poly(s.expand(A*A*Q.subs(a,(-B*n-C)/A)),n)
TA=s.Poly(s.cancel(A*A*T.subs(a,(-B*n-C)/A)),n)
q2,q1,q0=QA.all_coeffs();t2,t1,t0=TA.all_coeffs()
SYL=s.expand((q2*t0-q0*t2)**2-(q2*t1-q1*t2)*(q1*t0-q0*t1))
results['R_Sylvester_identity']=conv(SYL)==conv(A)**4*rd('R')
Nh=s.cancel(-203*(t2*q0-q2*t0)/A);Dh=s.cancel(-203*(q2*t1-t2*q1)/A)
results['D_chart_identity']=conv(Dh)==rd('D')
G=rd('cancelled');nn=A*Nh;aa=-(B*Nh+C*Dh);zz=A*Dh
ph=conv(s.expand(al*aa+be*nn-23*h*zz));et=conv(s.expand(4*(2*nn+3*aa)*(20*u*u+34*u-70)-105*aa*ell-381*h*zz))
results['E_cancellation_identity']=rd('E')*G==fe*ph-et
results['G_is_A_unit']=G==-conv(A)/3
R=rd('R');E=rd('E');Hr=rd('H_raw');H=rd('H_factor_1');SUP=rd('collision_support');disc=rd('H_main_discriminant')
results['image_resultant_identity']=R.resultant(E,0)==Hr
hf=Hr.factor(); prod=c.constant(hf[0])
for p,m in hf[1]: prod=prod*p**m
results['image_factorization_identity']=prod==Hr and Hr%H==0
results['discriminant_identity']=H.resultant(H.derivative(2),2)==disc
results['support_nonzero_squarefree']=bool(SUP) and SUP.gcd(SUP.derivative(1)).total_degree()==0
for name in ['A','D','phi']:
 p=rd(name) if name!='phi' else rd('phi_reconstructed')
 results['exception_resultant_'+name]=R.resultant(p,0)==rd('res_R_'+name)
 results['support_covers_'+name]=SUP%((rd('res_R_'+name)/rd('res_R_'+name).gcd(rd('res_R_'+name).derivative(1))))==0
results['support_covers_discriminant']=SUP%(disc/disc.gcd(disc.derivative(1)))==0
content=Hr/H;results['support_covers_content']=SUP%(content/content.gcd(content.derivative(1)))==0
results['elapsed_seconds']=time.perf_counter()-t
(P/'verification_results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results),flush=True)
assert all(v for k,v in results.items() if k!='elapsed_seconds')
