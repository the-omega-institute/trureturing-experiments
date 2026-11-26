"""给出同球面、同原壳图、同尺度读数的两参非相似有理变形。"""
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path
import sys


SEED=(7,17,23)
POLE=(17,17,17)
E1=(17,-7,0)
E2=(23,0,-7)


def dot(x,y):
    return sum(a*b for a,b in zip(x,y,strict=True))


def distance(x,y):
    z=tuple(a-b for a,b in zip(x,y,strict=True))
    return dot(z,z)


def deform(sites,a,b):
    t=tuple(a*i+b*j for i,j in zip(E1,E2,strict=True))
    v=dot(t,t)
    seed=tuple(((867-v)*x+2*867*y)/(867+v) for x,y in zip(SEED,t,strict=True))
    assert dot(seed,seed)==867
    moves=dict(zip(permutations(SEED),permutations(seed),strict=True))
    return [moves.get(x,x) for x in sites]


def verify(root):
    graph=run_path(str(root/"ch39-arc-relative-determinant.py"))
    read=run_path(str(root/"ch39-scale-free-geometric-pair.py"))["readout"]
    sites=graph["sites"]
    orbit=set(permutations(SEED))
    assert len(orbit)==6 and orbit<=set(sites)
    D={(i,j):distance(x,y) for (i,x),(j,y) in combinations(enumerate(sites),2)}
    moving={d for (i,j),d in D.items() if sites[i] in orbit or sites[j] in orbit}
    assert min(D.values())==12
    assert max(d for d in moving if d<=804)==632
    assert min(d for d in moving if d>804)==1064
    # 每站移动 < 1/2 时，每对距离移动 < 1。以下均为无根式的保守界。
    assert 632 < Fraction(51,2)**2 and Fraction(53,2)**2 < 724 < 804
    assert 1064>32**2 and 31**2>804
    near=set(permutations((13,13,23)))
    far=set(permutations((1,5,29)))
    assert {x for x in sites if distance(POLE,x)==68}==near
    assert {x for x in sites if distance(POLE,x)==544}==far
    assert {distance(POLE,x) for x in orbit}=={136}
    assert 68<9**2 and 11**2<136<12**2 and 13**2<544
    assert min(max(distance(x,z) for z in sites) for x in sites if x!=POLE)==804
    assert max(distance(POLE,z) for z in sites)==544
    assert 544<24**2 and 804>28**2
    assert dot(SEED,SEED)==867 and dot(SEED,E1)==dot(SEED,E2)==0
    assert dot(E1,E1)<19**2 and dot(E2,E2)<25**2 and Fraction(19+25,200)<Fraction(1,4)
    assert E1[2]==0 and E2[2]!=0 and E1[1]!=0
    pair,edges,sides,group=read(sites)
    assert pair==(32,724)
    changed=deform(sites,Fraction(1,200),Fraction(-1,200))
    assert sum(p!=q for p,q in zip(changed,sites,strict=True))==6 and read(changed)==(pair,edges,sides,group)
    i,j=sites.index(SEED),sites.index((13,13,23))
    # 基于极点的三角形，辐边几何相位为 1。cos² 去除方向与中心符号。
    def cosine_squared(points):
        x,y=points[i],points[j]
        real=867*(867+dot(POLE,x)+dot(x,y)+dot(y,POLE))
        imag=graph["det"](POLE,x,y)
        return Fraction(real**2,real**2+867*imag**2)
    pole=sites.index(POLE)
    assert all(tuple(sorted(e)) in edges for e in combinations((i,j,pole),2))
    c0,c1=cosine_squared(sites),cosine_squared(changed)
    change=abs(c1-c0)
    assert c1>c0 and change>Fraction(1,2500000)
    print("原活动轨道边的内/外平方弦长界",632,1064)
    print("有理参数、变形后的有序代表及尺度对",Fraction(1,200),Fraction(-1,200),changed[i],pair)
    print("固定三角面 cos² 的正改变量与下界",change,Fraction(1,2500000))


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
