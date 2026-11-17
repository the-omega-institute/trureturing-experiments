import json,time
from pathlib import Path
from flint import fmpq_mpoly_ctx,fmpq,fmpq_poly
P=Path(__file__).parent;c=fmpq_mpoly_ctx.get(('u','b','e'),'lex');u,b,e=c.gens();t=time.perf_counter();rows=[]
def rd(n):
 d=json.loads((P/(n+'.json')).read_text());return c.from_dict({tuple(k):fmpq(v) for k,v in d['terms']})
def save(n,p):
 (P/(n+'.json')).write_text(json.dumps({'variables':list(c.names()),'terms':[[[int(z) for z in k],str(v)] for k,v in p.to_dict().items()]})+'\n')
def out(label,**kw):
 row=dict(label=label,seconds=time.perf_counter()-t,**kw);rows.append(row);print(json.dumps(row,default=str),flush=True);(P/'support_results.json').write_text(json.dumps(rows,indent=2,default=str)+'\n')
R=rd('R');H=rd('H_factor_1');Hr=rd('H_raw')
# Exact root-support certificate: all denominators and leading-coefficient exceptions included.
lcR=c.from_dict({(0,k[1],0):v for k,v in R.to_dict().items() if k[0]==16})
lcH=c.from_dict({(0,k[1],0):v for k,v in H.to_dict().items() if k[2]==16})
content=Hr/H
parts=[lcR,lcH,content,rd('res_R_A'),rd('res_R_D'),rd('res_R_phi'),rd('H_main_discriminant')]
cert=c.constant(1)
for p in parts:
 sq=p/p.gcd(p.derivative(1));cert=cert*sq/cert.gcd(sq)
# normalize monic
lead=max(cert.to_dict(),key=lambda k:k[1]);cert=cert/c.constant(cert.to_dict()[lead]);save('collision_support',cert)
out('support',nonzero=bool(cert),degrees=cert.degrees(),terms=len(cert.to_dict()),lc_R=str(lcR),lc_H=str(lcH),b_content_degrees=content.degrees())
f=cert.factor();out('support_factors',factors=[dict(degrees=p.degrees(),multiplicity=m) for p,m in f[1]])
# Exact Sturm real-root counts, no numeric sampling. Implement using rational FLINT remainders.
def univ(p):
 deg=int(p.degrees()[1]);return fmpq_poly([p.to_dict().get((0,j,0),0) for j in range(deg+1)])
def sturm_count(poly):
 st=[poly,poly.derivative()]
 while st[-1]:
  r=-(st[-2]%st[-1])
  if not r:break
  # positive normalization controls coefficient swell, preserves signs
  r=r/abs(r.leading_coefficient());st.append(r)
 def variation(sgns):return sum(a!=b for a,b in zip(sgns,sgns[1:]))
 signs_p=[1 if x.leading_coefficient()>0 else -1 for x in st]
 signs_m=[sg*((-1)**x.degree()) for sg,x in zip(signs_p,st)]
 return variation(signs_m)-variation(signs_p),len(st)
for i,(p,m) in enumerate(f[1]):
 save('support_factor_'+str(i),p);ss=time.perf_counter();poly=univ(p);count,sl=sturm_count(poly)
 out('factor_real_root_count',index=i,degree=poly.degree(),count=count,sturm_length=sl,cost=time.perf_counter()-ss,polynomial=str(p) if poly.degree()<5 else None)
out('complete',total_cost=time.perf_counter()-t)
