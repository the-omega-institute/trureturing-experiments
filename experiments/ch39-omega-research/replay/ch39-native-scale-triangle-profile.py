"""十九点壳的钟弦交换自同构、最小度边的度量三角数据与面轨道数。"""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys


def triangle_lengths(edges, face):
    return tuple(sorted(edges[pair] for pair in combinations(face, 2)))


def verify(root):
    tomography = run_path(str(root / "ch39-mincut-real-triangle-tomography.py"))
    graph = tomography["graph"]
    sites = graph["sites"]
    edges = {(i,j):d for i,j,d,_ in graph["edges"]}
    degree = Counter(i for pair in edges for i in pair)
    assert len(sites)==19 and len(edges)==135
    assert min(degree.values()) == 10
    minimum = {i for i in range(len(sites)) if degree[i]==min(degree.values())}
    matched = [pair for pair in edges if set(pair)<=minimum]
    assert len(minimum)==6 and sorted(i for pair in matched for i in pair)==sorted(minimum)
    faces=tomography["enumerate_faces"]()
    assert len(faces)==529
    diagonals={edges[e] for e in matched}
    assert diagonals=={32}
    profiles=[]
    for edge in matched:
        profile=Counter(triangle_lengths(edges,f) for f in faces if set(edge)<=set(f))
        assert sum(profile.values())==9
        profiles.append(profile)
    expected={(32,244,244):1, (32,544,544):1, (32,152,152):1,
              (32,216,296):2, (32,504,632):2, (32,724,804):2}
    assert all(p==expected for p in profiles)
    # 共用底边时，在两条侧边长度的逐坐标偏序中找支配所有类型者。
    profile=profiles[0]
    dominating,=[s for s in profile if all(all(x>=y for x,y in zip(s,t,strict=True)) for t in profile)]
    assert dominating==(32,724,804)
    print("最低度顶点及其诱导匹配、支配三角与每条边上的重数",sorted(minimum),matched,dominating,profile[dominating])

    # 裸图的三对真双生顶点；同时交换后，将有向外点至内点的两钟互换。
    permutation=list(range(19))
    for i,j in matched:
        permutation[i],permutation[j]=j,i
        for k in range(19):
            if k not in (i,j):
                assert (tuple(sorted((i,k))) in edges)==(tuple(sorted((j,k))) in edges)
    assert {tuple(sorted((permutation[i],permutation[j]))) for i,j in edges}==set(edges)
    clock0,clock1=dominating[1:]
    assert {(permutation[i],permutation[j]) for (i,j),d in edges.items() if d==clock0}=={
           (i,j) for (i,j),d in edges.items() if d==clock1}
    print("交换两钟的图自同构",[(i,j) for i,j in matched])

    mixed=Counter()
    for f in faces:
        lens=triangle_lengths(edges,f)
        if clock0 in lens and clock1 in lens:
            mixed[lens]+=1
    assert mixed=={(32,724,804):6,(200,724,804):6}
    rotation=[sites.index((b,c,a)) for a,b,c in sites]
    rotate=lambda f:tuple(sorted(rotation[i] for i in f))
    moving={f:(m,z) for f,(m,z) in faces.items() if m!=(0,0)}
    labels=Counter(moving.values())
    assert len(moving)==96 and len(labels)==32 and set(labels.values())=={3}
    assert all(rotate(f)!=f and moving[rotate(f)]==label for f,label in moving.items())
    assert len({z for _,z in moving.values()})==32
    real_labels={(k,(z[0],abs(z[1]))) for k,z in moving.values()}
    assert len(real_labels)==16
    assert Counter(k for k,z in labels)=={(1,0):14,(0,1):14,(1,-1):4}
    print("混合三角型与重数、有钟面数、定向原生标签及合并共轭的标签数",mixed,len(moving),len(labels),len(real_labels))

    # sin²(arg z) 只依赖相位和其共轭类。原生面标签给出严格的有理值。
    moments=[]
    means=[]
    for k in ((1,0),(0,1)):
        pure=[z for m,z in moving.values() if m==k]
        assert len(pure)==42
        moments.append(sum(Fraction(867*b*b,a*a+867*b*b) for a,b in pure)/len(pure))
        means.append(sum(tomography["unit"](z).real for z in pure)/len(pure))
    assert moments[0]>moments[1] and means[0]<means[1]
    print("两钟纯面 sin² 的精确均值、平均余弦",moments,[x.str(14) for x in means])


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
