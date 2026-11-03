"""Bounded falsifiers for the actual-source evaluator; no all-cache certificate.

Failure criteria: a capacity/inversion disagreement, non-normalized split,
changed replay mass, false null probability, wrong mixed prediction, or adaptive
joint mass differing from the product of its conditional updates.
"""
from fractions import Fraction as F
from pathlib import Path
from random import Random
import hashlib
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import exact_checks as c
import exact_evaluator as e
assert sys.version_info >= (3,12)
rng=Random(62914)
records=[]
# Irregular masks with one whole endpoint anchored, so direct source evaluation
# remains small while the other endpoint has arbitrary non-prefix masks.
for anchor in ('left','right'):
 for trial in range(5):
  rows=[]
  for j in range(4):
   if anchor=='left':
    fixed=list(c.complete(19,0)[:7]); positions=range(7,15)
   else:
    fixed=list(c.complete(0,73)[7:]); positions=range(7)
   fixed += [(h,rng.randrange(2)) for h in positions if rng.randrange(3)!=0]
   rows.append(tuple(sorted(fixed)))
  D=tuple(rows)
  direct=e.normalizer(D); inv=c.evidence(D)
  assert direct==inv
  records.append({'anchor':anchor,'cache':[dict(row) for row in D],'Z':list(map(str,direct))})
# Redundant identical rows and complete replay use every original index once.
repeat=(c.complete(19,73),)*4
zr=e.normalizer(repeat)
assert zr==c.evidence(repeat)==(F(1,128*256*c.M**3),)*2
for h,bit in repeat[0]:
 assert c.conjoin(repeat,0,h,bit)==repeat
 assert c.conjoin(repeat,0,h,1-bit) is None
# A contrary response is D conjunction E, not a new random sample.
null=(tuple(sorted(repeat[0]+((0,1-dict(repeat[0])[0]),))),)+repeat[1:]
assert e.normalizer(null)==c.evidence(null)==(0,0)
# Empty evidence and both exact theta endpoints; the source mixture remains
# symbolic affine in theta and requires no exact-real numerical approximation.
assert c.evidence(((),)*4)==(1,1)
base=(c.complete(0,0),c.complete(0,1),c.complete(1,0),())
zbase=c.evidence(base)
closed=base[:3]+(c.complete(1,1),)
zclosed=c.evidence(closed)
mixed={}
for theta in (F(0),F(1,7),F(1,3),F(1,2),F(6,7),F(1)):
 p=c.mix(zclosed,theta)/c.mix(zbase,theta)
 assert p==F(37+54*theta,91*c.M)
 mixed[str(theta)]=str(p)
# A bounded two-invocation adaptive controller reads an old partial slot, then
# picks a different sample/coordinate based on that response. Constants and
# replay are identity updates; the finite cylinder ratios must telescope.
D=(((0,0),(7,1)),(),(),())
z=c.evidence(D)
leaves=[]
for bit in (0,1):
 first=c.conjoin(D,0,6,bit)
 zfirst=c.evidence(first)
 assert all(x>0 for x in zfirst)
 j,h=((1,14),(2,4))[bit]
 subtotal=[F(0),F(0)]
 for nextbit in (0,1):
  leaf=c.conjoin(first,j,h,nextbit)
  zl=c.evidence(leaf)
  for i in (0,1):
   conditional=zfirst[i]/z[i]*zl[i]/zfirst[i]
   assert conditional==zl[i]/z[i]
   subtotal[i]+=zl[i]
  leaves.append({'first_bit':bit,'next_bit':nextbit,'next_coordinate':[j,h],'mass_by_source':[str(zl[i]/z[i]) for i in (0,1)]})
 assert tuple(subtotal)==zfirst
for i in (0,1):
 assert sum(F(x['mass_by_source'][i]) for x in leaves)==1
# Verify the advertised run count separately from the distinct argument count.
args=[]
for width in (7,8):
 for r in (1,2,3,4):
  for trial in range(20):
   masks=tuple(tuple((h,(trial+b+h)%2) for h in range(width) if (trial+3*b+5*h)%4==0) for b in range(r))
   args.append((width,masks))
assert len(args)==160 and len(set(args))==32
# Additional genuinely distinct non-prefix label masks, including empty allowed
# sets, conflicting singleton choices, absent restrictions, and no blocks.
labelargs=set()
for width in (7,8):
 for r in range(5):
  for trial in range(7):
   masks=tuple(tuple((h,rng.randrange(2)) for h in range(width) if rng.randrange(3)==0) for b in range(r))
   labelargs.add((width,masks))
labelargs.update({(7,(None,)),(8,(None,)),(7,(((0,0),)*0,)),(7,(((0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(6,0)),)*2)})
for width,masks in labelargs:
 assert c.label_dp(masks,width)==c.label_inclusion_exclusion(masks,width)
result={'status':'passed','python':sys.version,'direct_vs_inversion_irregular_caches':records,'repeated_full_label_Z':list(map(str,zr)),'null_conjunction_Z':['0','0'],'theta_mixed_prediction':mixed,'bounded_adaptive_leaves':leaves,'supplied_label_comparisons':len(args),'supplied_distinct_label_arguments':len(set(args)),'additional_distinct_label_arguments':len(labelargs),'scope':'Finite falsifiers only; neither universal all-cache coverage nor Lean verification.'}
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('direct_vs_inversion_irregular_caches','bounded_adaptive_leaves')},sort_keys=True))
