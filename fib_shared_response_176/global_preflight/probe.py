import sys,time,json,resource
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import scan_fibers as v
import sympy as s
from flint import fmpq_mpoly_ctx, fmpq
ROOT=Path(__file__).parent
start=time.perf_counter(); events=[]
def emit(label,**kw):
    row=dict(label=label,seconds=time.perf_counter()-start,**kw); events.append(row)
    print(json.dumps(row,default=str),flush=True)
    (ROOT/'probe_results.json').write_text(json.dumps(events,indent=2,default=str)+'\n')
def save(name,p):
    (ROOT/(name+'.json')).write_text(json.dumps({'variables':list(ctx.names()),'terms':[[[int(z) for z in k],str(v)] for k,v in p.to_dict().items()]})+'\n')
def conv(ex):
    return ctx.from_dict({m:fmpq(int(c.p),int(c.q)) for m,c in s.Poly(s.expand(ex),u,b,e,domain=s.QQ).terms()})
u,b=v.u,v.b; e=s.symbols('e'); ctx=fmpq_mpoly_ctx.get(('u','b','e'),'lex'); fu,fb,fe=ctx.gens()
R=conv(v.R.as_expr()); N=v.nhat; D=v.dhat; A=v.A; B=v.B; C=v.C
# A*D*a=-(B*N+C*D); D*n=N. No division used below.
aa=-(B*N+C*D); nn=A*N; zz=A*D
ph=s.expand((3*v.h-10)*aa+(8*u*u+18*u-30)*nn-23*v.h*zz)
et=s.expand(4*(2*nn+3*aa)*(20*u*u+34*u-70)-105*aa*(u*u+3*u-4)-381*v.h*zz)
ph=conv(ph); et=conv(et); G=ph.gcd(et); ph=ph/G; et=et/G
E=fe*ph-et
emit('inputs',R_degrees=R.degrees(),R_total=R.total_degree(),R_terms=len(R.to_dict()),N_degrees=s.Poly(N,u,b).degree_list(),D_degrees=s.Poly(D,u,b).degree_list(),A_degrees=s.Poly(A,u,b).degree_list(),E_degrees=E.degrees(),cancelled= str(G),phi_degrees=ph.degrees(),eta_degrees=et.degrees())
for name,p in [('R',R),('E',E),('A',conv(A)),('D',conv(D)),('cancelled',G),('phi_reconstructed',ph)]: save(name,p)
# Exceptional chart support remains in the necessary certificate, not discarded.
for name,p in [('A',conv(A)),('D',conv(D)),('cancelled',G),('phi',ph)]:
    t=time.perf_counter(); z=R.resultant(p,0); save('res_R_'+name,z)
    emit('exception_resultant',name=name,degrees=z.degrees(),terms=len(z.to_dict()),nonzero=bool(z),factorization=str(z.factor()) if z and len(z.to_dict())<100 else None,cost=time.perf_counter()-t)
t=time.perf_counter(); H=R.resultant(E,0);save('H_raw',H)
emit('image_resultant',degrees=H.degrees(),terms=len(H.to_dict()),nonzero=bool(H),cost=time.perf_counter()-t,maxrss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
# Only b,e remain. Remove b-content and repeated algebraic factors explicitly.
t=time.perf_counter(); fac=H.factor(); emit('image_factorization',constant=str(fac[0]),factors=[dict(degrees=p.degrees(),terms=len(p.to_dict()),multiplicity=k) for p,k in fac[1]],cost=time.perf_counter()-t)
for i,(p,k) in enumerate(fac[1]):save('H_factor_'+str(i),p)
# The b-content is independent of e. Compute its raw resultant by homogeneity,
# avoiding a redundant high-degree Sylvester calculation.
main=next(p for p,m in fac[1] if p.degrees()[2]>0)
content=H/main
main_disc=main.resultant(main.derivative(2),2)
raw_disc=content**(2*int(main.degrees()[2])-1)*main_disc
save('H_disc_raw',raw_disc)
emit('raw_discriminant_from_content_identity',nonzero=bool(raw_disc),degrees=raw_disc.degrees(),terms=len(raw_disc.to_dict()),identity='Res_e(c*H,c*H_e)=c^31*Res_e(H,H_e)')
