"""由无标记点集恢复壳图、完整等距群及相似不变的整数尺度对。"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path
import sys


def volume(points):
    u,v,w=[[y-x for x,y in zip(points[0],p)] for p in points[1:]]
    return sum(u[i]*(v[(i+1)%3]*w[(i+2)%3]-v[(i+2)%3]*w[(i+1)%3]) for i in range(3))


def readout(points):
    vertices=range(len(points))
    # 点标签仅用于记录。所有选择只读点间距离。
    dist=[[sum((x-y)**2 for x,y in zip(p,q)) for q in points] for p in points]
    profiles=[tuple(sorted(row)) for row in dist]
    center,=[i for i in vertices if profiles.count(profiles[i])==1]
    groups={}
    for i,d in enumerate(dist[center]):
        if i!=center:
            groups.setdefault(d,[]).append(i)
    near=groups[min(groups)]
    far=groups[max(groups)]
    assert len(near)==3 and len(far)==6
    anchors=[center,*near]
    assert volume([points[i] for i in anchors])!=0
    signatures=[tuple(dist[i][j] for j in anchors) for i in vertices]
    index={s:i for i,s in enumerate(signatures)}
    assert len(index)==len(points)
    group=[]
    for perm in permutations(range(1,4)):
        # 唯一固定中心和最近三点决定整个等距变换。以四距离坐标解出逐点像。
        inv=[perm.index(i)+1 for i in range(1,4)]
        action=tuple(index[(s[0],*(s[j] for j in inv))] for s in signatures)
        assert action[center]==center and [action[j] for j in near]==[anchors[i] for i in perm]
        assert all(dist[i][j]==dist[action[i]][action[j]] for i,j in combinations(vertices,2))
        sign=1 if volume([points[action[i]] for i in anchors])*volume([points[i] for i in anchors])>0 else -1
        group.append((sign,action))
    assert len({a for s,a in group})==6 and sum(s==1 for s,a in group)==3

    base_length=min(dist[i][j] for i,j in combinations(far,2))
    bases=[(i,j) for i,j in combinations(far,2) if dist[i][j]==base_length]
    assert len(bases)==3 and {i for e in bases for i in e}==set(far)
    long_pairs=[]
    for i,j in bases:
        sides=[tuple(sorted((dist[i][k],dist[j][k]))) for k in near]
        maximal,=[pair for pair in set(sides) if all(
            pair[0]>=x and pair[1]>=y for x,y in sides)]
        assert sides.count(maximal)==2 and base_length<maximal[0]<maximal[1]
        long_pairs.append(maximal)
    shorter,longer=long_pairs[0]
    assert len(set(long_pairs))==1
    edges={(i,j) for i,j in combinations(vertices,2) if dist[i][j]<=longer}
    degree=Counter(i for edge in edges for i in edge)
    assert {i for i in vertices if degree[i]==min(degree.values())}==set(far)
    clocks=[]
    for d in (shorter,longer):
        by_length=[e for e in edges if dist[e[0]][e[1]]==d]
        assert len(by_length)==6 and all(len(set(e)&set(far))==1 and len(set(e)&set(near))==1 for e in by_length)
        clocks.append({(i,j) if i in far else (j,i) for i,j in by_length})
    for _,p in group:
        assert all({(p[i],p[j]) for i,j in clock}==clock for clock in clocks)
    # 环路净计数以相反遍历同时变号。符号正规化不影响非零谓词。
    frequencies={}
    for i,j,k in combinations(vertices,3):
        if {(i,j),(j,k),(i,k)}<=edges:
            cycle=((i,j),(j,k),(k,i))
            freq=tuple(sum(int((a,b) in c)-int((b,a) in c) for a,b in cycle) for c in clocks)
            frequencies[i,j,k]=max(freq,tuple(-x for x in freq))
    proper=Counter()
    full=Counter()
    rest=set(frequencies)
    orbit_sizes=Counter()
    while rest:
        tri=rest.pop()
        orbit={tuple(sorted(p[i] for i in tri)) for sg,p in group if sg>0}
        assert all(frequencies[t]==frequencies[tri] for t in orbit)
        proper[frequencies[tri]]+=1
        orbit_sizes[frequencies[tri],len(orbit)]+=1
        rest.difference_update(orbit)
    rest=set(frequencies)
    while rest:
        tri=rest.pop()
        orbit={tuple(sorted(p[i] for i in tri)) for _,p in group}
        assert all(frequencies[t]==frequencies[tri] for t in orbit)
        full[frequencies[tri]]+=1
        rest.difference_update(orbit)
    fixed=sorted((sign,sum(t==tuple(sorted(p[i] for i in t)) for t in frequencies),
                  sum(t==tuple(sorted(p[i] for i in t)) for t,f in frequencies.items() if f!=(0,0)))
                 for sign,p in group)
    multiplicities=Counter(frequencies.values())
    assert len(edges)==135 and len(frequencies)==529
    assert multiplicities=={(0,0):433,(1,0):42,(0,1):42,(1,-1):12}
    assert proper=={(0,0):147,(1,0):14,(0,1):14,(1,-1):4}
    assert full=={(0,0):83,(1,0):7,(0,1):7,(1,-1):2}
    assert fixed==[(-1,19,0)]*3+[(1,4,0)]*2+[(1,529,96)]
    assert all(n==3 for (f,n) in orbit_sizes if f!=(0,0))
    K=sum(n for k,n in proper.items() if k!=(0,0))
    D=K*Fraction(shorter,base_length)
    return (K,D),edges,(base_length,shorter,longer),fixed


def verify(root):
    graph=run_path(str(root/"ch39-arc-relative-determinant.py"))
    points=graph["sites"]
    pair,edges,sides,fixed=readout(points)
    assert edges=={(i,j) for i,j,d,z in graph["edges"]}
    assert sides==(32,724,804) and pair==(32,724)
    # 同时检验重编号、坐标反射、平移与有理缩放。任意相似度的证明由距离比及共轭群给出。
    altered=[tuple(Fraction(7,11)*z+b for z,b in zip((y,x,t),(3,-2,7))) for x,y,t in reversed(points)]
    altered_pair,_,altered_sides,_=readout(altered)
    assert altered_pair==pair and altered_sides==tuple(Fraction(49,121)*s for s in sides)
    print("完整等距群的定向符号、全部/活动三角固定数",fixed)
    print("重建的边数、极值三角平方边长及缩放后平方边长",len(edges),sides,altered_sides)
    print("保向活动三角轨道数、短侧平方长与底边平方长的比、尺度对",pair[0],Fraction(sides[1],sides[0]),(pair[0],int(pair[1])))


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
