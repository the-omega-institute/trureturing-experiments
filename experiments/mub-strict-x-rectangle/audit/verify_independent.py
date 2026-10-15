#!/usr/bin/env python3
"""Portable independent root, center-symmetry, graph, and replay-report audit.

Standard Python only. Acceptance arithmetic uses arbitrary-precision integers
and exact Fractions. Supplied preconditioners are untrusted proposals: every use
is checked through a fresh contraction/inclusion calculation. This is an
external computational audit, not Lean admission or an independent cover trace.
"""
from dataclasses import dataclass
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
import argparse,hashlib,json
S=1<<40
RADIUS=1<<28
OLD_RADIUS=1<<24
@dataclass(frozen=True)
class D:
 lo:int
 hi:int
 def __post_init__(self):assert self.lo<=self.hi
 @classmethod
 def number(cls,n):return cls(n*S,n*S)
 def __add__(self,o):
  o=o if isinstance(o,D)else D.number(o)
  return D(self.lo+o.lo,self.hi+o.hi)
 __radd__=__add__
 def __neg__(self):return D(-self.hi,-self.lo)
 def __sub__(self,o):return self+-o
 def __rsub__(self,o):return -self+o
 def __mul__(self,o):
  o=o if isinstance(o,D)else D.number(o)
  products=[a*b for a in(self.lo,self.hi)for b in(o.lo,o.hi)]
  return D(min(products)//S,-((-max(products))//S))
 __rmul__=__mul__
 def reciprocal(self):
  assert not self.lo<=0<=self.hi
  return D(S*S//self.hi,-((-S*S)//self.lo))
 def __truediv__(self,o):
  o=o if isinstance(o,D)else D.number(o)
  return self*o.reciprocal()
 def sq(self):
  m=0 if self.lo<=0<=self.hi else min(self.lo*self.lo,self.hi*self.hi)
  M=max(self.lo*self.lo,self.hi*self.hi)
  return D(m//S,-((-M)//S))
 def absmax(self):return max(-self.lo,self.hi)
@dataclass(frozen=True)
class Z:
 r:D
 i:D
 @classmethod
 def number(cls,n):return cls(D.number(n),D.number(0))
 def __add__(self,o):
  o=o if isinstance(o,Z)else Z.number(o)
  return Z(self.r+o.r,self.i+o.i)
 __radd__=__add__
 def __mul__(self,o):
  o=o if isinstance(o,Z)else Z.number(o)
  return Z(self.r*o.r-self.i*o.i,self.r*o.i+self.i*o.r)
 def conj(self):return Z(self.r,-self.i)
 def normsq(self):return self.r.sq()+self.i.sq()
 def __neg__(self):return Z(-self.r,-self.i)
 def __sub__(self,o):return self+-o
 def reciprocal(self):
  n=self.normsq();return Z(self.r/n,-self.i/n)
 def __truediv__(self,o):
  o=o if isinstance(o,Z)else Z.number(o)
  return self*o.reciprocal()

def require(ok,msg):
 if not ok:raise ValueError(msg)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_bounds(path):
 lines=path.read_text().splitlines();require(len(lines)==37 and lines[0]=='40','matrix bounds format')
 flat=[]
 for ln in lines[1:]:
  v=list(map(int,ln.split()));require(len(v)==4,'matrix bound entry')
  a,b,c,d=v;flat.append(Z(D(a,b),D(c,d)))
 return [flat[6*i:6*i+6]for i in range(6)]
def read_centres(path):
 lines=path.read_text().splitlines();require(len(lines)==61 and lines[0]=='60','center count')
 ans={}
 for ln in lines[1:]:
  v=list(map(int,ln.split()));require(len(v)==7,'center row')
  idx,mask,*coords=v;require(0<=idx<60 and idx not in ans and 0<=mask<32,'center labels')
  require(all(-S<=t<=S for t in coords),'compact centers')
  ans[idx]=(mask,coords)
 return ans
def read_proposals(path,centres):
 lines=path.read_text().splitlines();require(len(lines)==61 and lines[0]=='60','proposal count')
 ans={}
 for ln in lines[1:]:
  v=list(map(int,ln.split()));require(len(v)==27,'proposal matrix')
  idx,mask=v[:2];require(idx in centres and idx not in ans and mask==centres[idx][0],'proposal labels')
  ans[idx]=[[D(x,x)for x in v[2+5*i:7+5*i]]for i in range(5)]
 return ans
def phase(X,mask):
 ans=[Z.number(1)]
 for j,t in enumerate(X):
  sign=1-2*((mask>>j)&1);den=1+t.sq()
  ans.append(Z(sign*(1-t.sq())/den,2*sign*t/den))
 return ans
def evaluate(H,X,mask):
 u=phase(X,mask);du=[]
 for j,t in enumerate(X):
  sign=1-2*((mask>>j)&1);den=1+t.sq()
  du.append(Z((-4*sign)*t/den.sq(),(2*sign)*(1-t.sq())/den.sq()))
 f=[];J=[]
 for a in range(5):
  y=sum((H[n][a].conj()*u[n]for n in range(6)),Z.number(0))
  f.append(y.normsq()-6);row=[]
  for j in range(5):
   z=H[j+1][a].conj()*du[j];row.append(2*(y.r*z.r+y.i*z.i))
  J.append(row)
 return f,J
def krawczyk(H,X,mask,P):
 m=[D((t.lo+t.hi)//2,(t.lo+t.hi)//2)for t in X]
 f0,_=evaluate(H,m,mask);_,J=evaluate(H,X,mask)
 K=[];q=0
 for i in range(5):
  k=m[i]
  for n in range(5):k=k-P[i][n]*f0[n]
  row=0
  for j in range(5):
   e=D.number(int(i==j))
   for n in range(5):e=e-P[i][n]*J[n][j]
   row+=e.absmax();k=k+e*(X[j]-m[j])
  q=max(q,row);K.append(k)
 return K,q

def validate_roots(H,centres,proposals,radius):
 roots={};maxq=0;maxd=0
 for idx in range(60):
  mask,coords=centres[idx];X=[D(t-radius,t+radius)for t in coords]
  K,q=krawczyk(H,X,mask,proposals[idx]);d=max((k-D(t,t)).absmax()for k,t in zip(K,coords))
  require(q<S and all(a.lo<b.lo<=b.hi<a.hi for a,b in zip(X,K)),f'root inclusion {idx}')
  maxq=max(maxq,q);maxd=max(maxd,d);Y=K
  for _ in range(3):
   K,_=krawczyk(H,Y,mask,proposals[idx]);Y=[D(max(a.lo,b.lo),min(a.hi,b.hi))for a,b in zip(Y,K)]
  roots[idx]=(mask,Y)
 return roots,{'max_contraction_num':maxq,'max_displacement_num':maxd,'scale':S}

def explicit_center():
 # The exact algebraic field is Q(i,sqrt(21)); coefficients use basis
 # 1,i,sqrt(21),i*sqrt(21). Independently check the seed Gram and symmetry.
 zero=(F(0),)*4;one=(F(1),F(0),F(0),F(0))
 def add(a,b):return tuple(x+y for x,y in zip(a,b))
 def neg(a):return tuple(-x for x in a)
 def conj(a):return(a[0],-a[1],a[2],-a[3])
 def mul(a,b):
  v=[F(0)]*4
  for j in range(4):
   for k in range(4):
    i=j%2+k%2;r=j//2+k//2
    v[i%2+2*(r%2)]+=a[j]*b[k]*(-1)**(i//2)*21**(r//2)
  return tuple(v)
 b=(F(-3,5),F(4,5),F(0),F(0));e=(F(-2,5),F(0),F(0),F(1,5))
 A=[[b if i==j else one for j in range(3)]for i in range(3)]
 B=[[e if i==j else one for j in range(3)]for i in range(3)]
 H=[A[i]+B[i]for i in range(3)]+[[conj(B[j][i])for j in range(3)]+[neg(conj(A[j][i]))for j in range(3)]for i in range(3)]
 require(all(mul(z,conj(z))==one for row in H for z in row),'exact center unit entries')
 for i in range(6):
  for j in range(6):
   g=zero;h=zero
   for k in range(6):g=add(g,mul(H[i][k],conj(H[j][k])));h=add(h,mul(conj(H[k][i]),H[k][j]))
   target=(F(6),F(0),F(0),F(0))if i==j else zero
   require(g==h==target,'exact center Grams')
 shift=[1,2,0,4,5,3]
 require(all(H[shift[i]][j]==H[i][shift.index(j)]for i in range(6)for j in range(6)),'exact shift commutator')
 r=isqrt(21*S*S);require(r*r<21*S*S<(r+1)*(r+1),'sqrt enclosure')
 sq=D(r,r+1)
 bd=Z(D.number(-3)/5,D.number(4)/5);ed=Z(D.number(-2)/5,sq/5)
 A=[[bd if i==j else Z.number(1)for j in range(3)]for i in range(3)]
 B=[[ed if i==j else Z.number(1)for j in range(3)]for i in range(3)]
 return [A[i]+B[i]for i in range(3)]+[[B[j][i].conj()for j in range(3)]+[-A[j][i].conj()for j in range(3)]for i in range(3)]

def center_fixation(roots,centres,canonical):
 for idx in canonical:
  mask,Y=roots[idx];v=phase(Y,mask);w=[v[j]for j in(1,2,0,4,5,3)]
  require(w[0].normsq().lo>0,'center projective denominator')
  w=[z/w[0]for z in w];coords=centres[idx][1]
  for j in range(5):
   z=w[j+1]*Z.number(1-2*((mask>>j)&1))
   require((z+1).normsq().lo>0,'center inverse chart denominator')
   t=Z(D.number(0),D.number(-1))*(z-1)/(z+1)
   require(t.i.lo<=0<=t.i.hi,'real inverse chart enclosure')
   require(coords[j]-OLD_RADIUS<t.r.lo<=t.r.hi<coords[j]+OLD_RADIUS,f'center shift fixation {idx}')
 return True

def pair_graph(roots,data):
 rays={i:phase(X,mask)for i,(mask,X)in roots.items()}
 edges=set(map(tuple,data['allowed_edges']));can=set(data['canonical_vertices']);colors=data['noncanonical_colors']
 require(len(edges)==114 and len(can)==6 and len(colors)==60,'graph size')
 require(all(0<=i<j<60 for i,j in edges)and all(0<=i<60 for i in can),'graph labels')
 require(all(colors[i]in(0,1)for i in range(60)if i not in can),'graph colors')
 require(all((i in can and j in can)or(i not in can and j not in can and colors[i]!=colors[j])for i,j in edges),'graph partition')
 low=S;worst=None;high=0;count=0
 for i,j in combinations(range(60),2):
  ip=sum((a.conj()*b for a,b in zip(rays[i],rays[j])),Z.number(0));n=ip.normsq()/36
  require(n.hi<S,f'projective distinctness {i},{j}');high=max(high,n.hi)
  if(i,j)not in edges:
   require(n.lo*500000>S,f'nonedge threshold {i},{j}')
   if n.lo<low:low=n.lo;worst=[i,j]
   count+=1
 require(count==1656,'nonedge count')
 return {'distinct_pair_tests':1770,'nonedges':count,'minimum_nonedge_lower':str(F(low,S)),'worst_pair':worst,'maximum_distinctness_upper':str(F(high,S))}

def read_enclosures(path):
 lines=path.read_text().splitlines();require(len(lines)==61 and lines[0]=='60','root-enclosure count')
 ans={}
 for line in lines[1:]:
  v=list(map(int,line.split()));require(len(v)==12,'root-enclosure row');idx,mask=v[:2]
  require(0<=idx<60 and idx not in ans and 0<=mask<32,'root-enclosure label')
  ans[idx]=(mask,[D(v[2+2*j],v[3+2*j])for j in range(5)])
 return ans

def run(result,input_path,proposal_path,source):
 data=json.loads(input_path.read_text());report=json.loads((result/'verification.json').read_text())
 require(data['schema']=='x-uniform-supergraph-patch-v1','input schema')
 require(report['status']=='PASS'and report['parameter_center']==['-1/5','0'],'result status/center')
 require(F(report['parameter_radius'])==F(1,1<<24)and F(report['uniqueness_radius'])==F(1,1<<12),'result radii')
 source_names={'driver':'verify_enlarged_patch.py','historical_core':'check_real_x_global_cover.cpp','enlarged_core':'cover_larger_patch.cpp','historical_driver':'check_real_x_global_cover.py','parameter_generator':'check_real_x_supergraph_patch.py','fraction_kernel':'check_strict_x_counterexample.py'}
 for key,name in source_names.items():require(sha(source/name)==report['hashes'][key],f'source hash {key}')
 require(sha(input_path)==report['hashes']['input'],'input hash')
 oldtext=(source/'check_real_x_global_cover.cpp').read_text()
 require(oldtext.count('1LL<<24')==2 and(source/'cover_larger_patch.cpp').read_text()==oldtext.replace('1LL<<24','1LL<<28'),'minimal source change')
 centres_path=result/'centre_preflight/centers.txt'
 for file,key in [(centres_path,'centres'),(result/'matrix_bounds.txt','matrix_bounds'),(result/'root_enclosures.txt','root_enclosures')]:require(sha(file)==report['hashes'][key],f'data hash {key}')
 charts=[]
 for k in range(32):
  c=json.loads((result/f'chart_{k:02d}.json').read_text())
  require(c==report['charts'][k]and c['chart']==k and c['status']=='COVERED'and c['pending']==c['unresolved']==0,f'chart {k}')
  require(c['nodes']==2*(c['excluded']+c['known_root_leaves'])+c['contracted']-1,f'tree count {k}')
  charts.append(c)
 require(sum(c['nodes']for c in charts)==report['nodes'],'aggregate nodes')
 centres=read_centres(centres_path);P=read_proposals(proposal_path,centres);H=read_bounds(result/'matrix_bounds.txt')
 newroots,margins=validate_roots(H,centres,P,RADIUS)
 require(margins['max_contraction_num']==385012914919 and margins['max_displacement_num']==206508646,'quoted margins')
 require(margins['max_contraction_num']*1000<351*S and margins['max_displacement_num']*100<77*RADIUS,'safety ratios')
 # Independent fixed-P refinements do not trust the C++ sharpened enclosures.
 independent_graph=pair_graph(newroots,data)
 reported_graph=pair_graph(read_enclosures(result/'root_enclosures.txt'),data)
 H0=explicit_center()
 require(all(a.r.lo<=b.r.lo<=b.r.hi<=a.r.hi and a.i.lo<=b.i.lo<=b.i.hi<=a.i.hi for ar,br in zip(H,H0)for a,b in zip(ar,br)),'exact center lies in new bounds')
 centerroots,center_margins=validate_roots(H0,centres,P,OLD_RADIUS)
 require(center_fixation(centerroots,centres,set(data['canonical_vertices'])),'fresh center fixation')
 require(OLD_RADIUS<RADIUS,'nested boxes')
 return {'status':'INDEPENDENT_PORTABLE_AUDIT_PASS','parameter_radius':'1/16777216','root_radius':'1/4096','root_count':60,'all_32_chart_reports_checked':True,'nodes':report['nodes'],'pending':0,'unresolved':0,'exact_center_Grams_and_shift_commutation_checked':True,'fresh_center_root_and_shift_fixation_checked':True,'old_boxes_strictly_contained_in_new':True,'enlarged_root_margins':margins,'center_root_margins':center_margins,'independently_refined_graph':independent_graph,'reported_enclosure_graph':reported_graph,'canonical_propagation':'Fresh exact-center fixation plus connected continuous unique-root branch lemma; see proof documentation','global_coverage_level':'Audited fresh reports and reviewed core, not an independent trace replay','arithmetic':'Independent unbounded-integer dyadic intervals; exact Fraction center-field identities','lean_kernel_verified':False,'sha256':{'independent_checker':sha(Path(__file__)),'preconditioner_proposals':sha(proposal_path),'combined_cover_report':sha(result/'verification.json'),'input':sha(input_path),'parameter_bounds':sha(result/'matrix_bounds.txt'),'centers':sha(centres_path)}}

if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--result',type=Path,required=True)
 ap.add_argument('--input',type=Path,required=True)
 ap.add_argument('--preconditioners',type=Path,required=True)
 ap.add_argument('--source',type=Path,help='Source directory; defaults to the input file parent')
 ap.add_argument('--output',type=Path,help='Optional JSON output; stdout always receives the report')
 a=ap.parse_args()
 require(__debug__,'run without Python -O')
 r=run(a.result,a.input,a.preconditioners,a.source or a.input.parent)
 text=json.dumps(r,indent=2)+'\n'
 if a.output:a.output.write_text(text)
 print(text,end='')
