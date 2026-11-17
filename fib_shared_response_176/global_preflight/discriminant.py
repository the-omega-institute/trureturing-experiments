import json,time,resource
from pathlib import Path
from flint import fmpq_mpoly_ctx,fmpq,fmpq_poly
P=Path(__file__).parent
c=fmpq_mpoly_ctx.get(('u','b','e'),'lex');u,b,e=c.gens();t=time.perf_counter(); rows=[]
def rd(n):
 d=json.loads((P/(n+'.json')).read_text());return c.from_dict({tuple(k):fmpq(v) for k,v in d['terms']})
def save(n,p):
 (P/(n+'.json')).write_text(json.dumps({'variables':list(c.names()),'terms':[[[int(z) for z in k],str(v)] for k,v in p.to_dict().items()]})+'\n')
def out(label,**kw):
 row=dict(label=label,seconds=time.perf_counter()-t,**kw);rows.append(row);print(json.dumps(row,default=str),flush=True);(P/'discriminant_results.json').write_text(json.dumps(rows,indent=2,default=str)+'\n')
H=rd('H_factor_1'); ht=H.to_dict();sp0=fmpq_poly([ht.get((0,0,j),0) for j in range(17)]);g=sp0.gcd(sp0.derivative())
out('b_zero_squarefree_certificate',degree=sp0.degree(),gcd=str(g),nonzero_discriminant=bool(g.degree()==0),discriminant=str(sp0.resultant(sp0.derivative())))
start=time.perf_counter();D=H.resultant(H.derivative(2),2);save('H_main_discriminant',D)
out('main_discriminant',nonzero=bool(D),degrees=D.degrees(),terms=len(D.to_dict()),cost=time.perf_counter()-start,maxrss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
start=time.perf_counter();f=D.factor();out('main_discriminant_factorization',constant=str(f[0]),factors=[dict(degrees=p.degrees(),multiplicity=m,terms=len(p.to_dict())) for p,m in f[1]],cost=time.perf_counter()-start)
for i,(p,m) in enumerate(f[1]):save('disc_factor_'+str(i),p)
