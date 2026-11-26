"""单次两面更新在时间与空间顶点位移上的三阶 Ward 障碍。"""
from pathlib import Path
from runpy import run_path
from itertools import permutations,combinations
from collections import defaultdict
import sympy as S
F=S.Rational
flat=run_path(str(Path(__file__).with_name("ch39-regge-flat-cell.py")));edges=flat['EDGES'];parts=flat['canonical']()
q=[F((j-i)*(5-j+i),5) for i,j in edges]
verts={m:S.Matrix([(m>>i)&1 for i in range(4)]) for m in range(16)}
e=S.Matrix([1,-1,0,0]); f=S.Matrix([0,0,1,-1]); time=S.ones(4,1); space=S.Matrix([1,1,-1,-1])
cross=e*f.T+f*e.T
wave=space+S.I*time/S.sqrt(5)
# 两条交叉腿取 k，顶点位移腿取 −2k，分别使用时间和空间协向量。
gauge_time=-2*(wave*time.T+time*wave.T)
gauge_space=-2*(wave*space.T+space*wave.T)
Q=[[S.expand((v.T*h*v)[0]) for v in verts.values()] for h in [cross,gauge_time,gauge_space]]
wavepos=[S.expand((v.T*wave)[0]) for v in verts.values()]
local={}
for face,c,m in parts:
 A=S.sqrt(F(1,5) if c==0 else F(3,20)); ie=[edges.index(e) for e in combinations(face,2)]
 area=[(sum(q[a] for a in ie)/2-q[k])/(8*A) if k in ie else 0 for k in range(10)]
 ang=[S.expand(sum(q[i]*m[i,k] for i in range(10))*S.sqrt(5)/A) for k in range(10)]
 local[face]=(A,area,ang)
DD=defaultdict(lambda:[0,0,0]);DA={};partners=defaultdict(lambda:[[],[]])
for perm in permutations(range(4)):
 v=[sum(1<<i for i in perm[:j]) for j in range(5)]
 dr=[v[j]-v[i] for i,j in edges]
 mid=[(wavepos[v[i]]+wavepos[v[j]])/2 for i,j in edges]
 triangles={}
 for face,(A,area,ang) in local.items():
  i,j,k=face; key=v[j]-v[i],v[k]-v[j]
  deriv=tuple(S.expand(sum(a*Q[l][m] for a,m in zip(area,dr))) for l in range(3))
  if key in DA:assert DA[key]==(A,deriv)
  DA[key]=(A,deriv)
  triangles[face]=(key,A)
  for l in (0,):
   DD[key][l]+= S.expand(-sum(c*Q[l][m]*(z-wavepos[v[i]])**2/2 for c,m,z in zip(ang,dr,mid)))
 for face,(key,A) in triangles.items():
  for g,(other,B) in triangles.items():
   if A!=B:partners[key][len(set(face)&set(g))-1].append(other)
res=[[0,0],[0,0]]
for key,(a,ds) in DA.items():
 dx,dg,dp=ds
 for group,items in enumerate(partners[key]):
  values=[0,0]
  for other in items:
   b,(ex,et,es)=DA[other]
   for j,pair in enumerate((dg+et,dp+es)):
    values[j]+=2*(dx+ex)*pair*DD[key][0]
  for j,value in enumerate(values):
   res[j][group]-=b*b/a*value/len(items)

a=S.sqrt(3)
expected=[[4*S.I*(-171+50*a)/45,4*S.I*(86*a+183)/225],
          [-304*S.sqrt(15)/1125-4*S.sqrt(5)/25,84*S.sqrt(5)/125+872*S.sqrt(15)/1125]]
required=[]
for values, row in zip(res,expected):
 assert all(S.simplify(v-e)==0 for v,e in zip(values,row))
 vertex,edge=row
 required.append(S.simplify(edge/(edge-vertex)))
assert 0<required[0]<S.Rational(45,100)<S.Rational(76,100)<required[1]<1
print("两个独立规范插入各自要求的路由概率",*[v.evalf(13) for v in required])
