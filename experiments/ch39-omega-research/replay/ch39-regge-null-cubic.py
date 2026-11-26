"""平直 Regge 三阶角亏导数与两个在壳外腿的纯规范插入。"""
from pathlib import Path

import sympy as sp
from itertools import combinations, permutations
from collections import defaultdict
from runpy import run_path
flat=run_path(str(Path(__file__).with_name("ch39-regge-flat-cell.py")))
edges=flat['EDGES']; E=range(10); q=[sp.Rational((j-i)*(5-j+i),5) for i,j in edges]
R=sp.Rational
def gram(v):
 def at(i,j):
  return 0 if i==j else v[edges.index(tuple(sorted((i,j))))]
 return sp.Matrix(4,4,lambda i,j:(at(0,i+1)+at(0,j+1)-at(i+1,j+1))/2)
grad=(-sp.ones(1,4)).col_join(sp.eye(4))
inv=gram(q).inv(); frame=grad*inv
G=[gram([sp.S.One if x==j else sp.S.Zero for x in E]) for j in E]
df=[-frame*m*frame.T for m in G]
ddf=[[frame*(G[b]*inv*G[c]+G[c]*inv*G[b])*frame.T for c in E] for b in E]
T=defaultdict(lambda:0)
for u,v in edges:
 face=list(combinations([i for i in range(5) if i not in (u,v)],2))
 cols=[edges.index(e) for e in face]
 cs=R(1,2) if v-u in (1,4) else 0
 A=sp.sqrt(R(3,20) if cs else R(1,5)); sine=sp.sqrt(1-cs**2)
 dx=[(sum(q[j] for j in cols)/2-q[a])/4 if a in cols else 0 for a in E]
 gg=[-(m[u,u]+m[v,v])/4 for m in df]
 dc=[-m[u,v]/2+cs*g for m,g in zip(df,gg)]
 for b in E:
  for c in E:
   m=ddf[b][c]
   g2=-(m[u,u]+m[v,v])/4+(df[b][u,u]*df[c][u,u]+df[b][v,v]*df[c][v,v])/8
   c2=-(m[u,v]+df[b][u,v]*gg[c]+df[c][u,v]*gg[b])/2+cs*(g2+gg[b]*gg[c])
   d_angle=c2/sine+cs*dc[b]*dc[c]/sine**3
   for a in cols:
    x2=R(1,8)*(1 if a!=b else -1) if b in cols else 0
    T[a,b,c]+=sp.expand((dx[a]*d_angle/(2*A)+2*(x2/(2*A)-dx[a]*dx[b]/(4*A**3))*dc[c]/sine)/sp.sqrt(5))
# 两条交叉偏振腿与一条精确顶点位移腿的三阶系数，先保存动量的齐次多项式。
z=sp.Symbol('z')
vector=[(int(i<2)-int(i>=2)) for i in range(4)]
def xy(v):
 return (sum(a for i,a in enumerate(vector) if v & (1<<i)),v.bit_count())
cross={5:2,6:-2,9:-2,10:2}
acc=[defaultdict(lambda:0),defaultdict(lambda:0)]
for order in permutations(range(4)):
 vert=[sum(1<<i for i in order[:j]) for j in range(5)]
 data=[]
 for i,j in edges:
  a,b=vert[i],vert[j];lx,lt=xy(a);rx,rt=xy(b)
  data.append((rx-lx,rt-lt,lx+rx,lt+rt,cross.get(b-a,0)))
 for i,(x,t,m,n,_) in enumerate(data):
  for j,(_,_,m1,n1,c1) in enumerate(data):
   if not c1:continue
   for k,(_,_,m2,n2,c2) in enumerate(data):
    if not c2:continue
    cx,ct=m1+m2-2*m,n1+n2-2*n
    for idx,xi in enumerate([t,x]):
     acc[idx][x,t,cx,ct]+=T[i,j,k]*4*c1*c2*xi
for direction, sums in enumerate(acc):
 poly=sp.factor(sum(val*(x+t*z)*((a+b*z)**2/8+(x+t*z)**2/6) for (x,t,a,b),val in sums.items()))
 expected = 64*z*(5*z**2+1) if direction==0 else 64*(5*z**2+1)/5
 assert sp.expand(poly-expected)==0
 assert sp.simplify(poly.subs(z,sp.I/sp.sqrt(5)))==0
 print('纯 Regge 的二导数 Ward 系数',poly)
