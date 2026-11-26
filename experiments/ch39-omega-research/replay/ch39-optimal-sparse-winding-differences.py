"""以第6枚星图及跨号运输边达到全5050项能差的锐误差因子。"""
import argparse
from hashlib import sha256
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx


def move_direct(residual, i, j):
    """将正负两端可配的质量移到一条有向边，返回其精确式的区间包围。"""
    if residual[i] < 0 < residual[j]:
        i, j = j, i
    if not (residual[i] > 0 > residual[j]):
        return arb(0)
    if residual[i] < -residual[j]:
        amount = residual[i]
        residual[i] = arb(0)
        residual[j] += amount
    else:
        assert residual[i] > -residual[j]
        amount = -residual[j]
        residual[j] = arb(0)
        residual[i] -= amount
    return amount


def optimal_tree(row):
    positive = [j for j,v in enumerate(row) if v > 0]
    negative = [j for j,v in enumerate(row) if v < 0]
    assert len(positive)+len(negative)==len(row)
    residual = list(row)
    tree = set()
    while positive and negative:
        i,j = positive[-1],negative[-1]
        tree.add(tuple(sorted((i,j))))
        if len(positive)==len(negative)==1:
            # 行和精确为零，最后两份流量相等。
            assert (residual[i]+residual[j]).contains(0)
            break
        move_direct(residual,i,j)
        if residual[i]==0:
            positive.pop()
        if residual[j]==0:
            negative.pop()
    assert len(tree)==100
    return tree


def direct_then_hub(row, hub, tree):
    residual=list(row)
    cost=arb(0)
    for i,j in sorted(tree):
        if hub not in (i,j):
            cost+=move_direct(residual,i,j)
    return cost,[abs(v) for v in residual]


def hub_distances(edges, hub, n):
    adj=[[] for _ in range(n)]
    for i,j in edges:
        adj[i].append(j)
        adj[j].append(i)
    distance=[-1]*n
    distance[hub]=0
    queue=[hub]
    for i in queue:
        for j in adj[i]:
            if distance[j]==-1:
                distance[j]=distance[i]+1
                queue.append(j)
    assert min(distance)==0
    return distance


def verify(replay, show_edges):
    # 被复用的前缀证书从第一个位置参数解析依赖目录。
    sys.argv=[sys.argv[0],str(replay)]
    prefix=run_path(str(replay/"ch39-winding-difference-prefix.py"))
    ctx.prec=320
    equality,full,_,_=prefix["rows"]()
    rows,gain=prefix["inverse_and_gain"](equality,full[:101],[[1]]*101)
    norms=[sum(abs(x) for x in row)/2 for row in rows]
    q=max(range(105),key=lambda i:float(norms[i]))
    target=norms[q]
    assert q+1==79 and all(target>v for i,v in enumerate(norms) if i!=q)
    assert (gain/2).overlaps(target)
    assert arb("8.889184950660")<target<arb("8.889184950662")
    hub=5
    star={tuple(sorted((hub,j))) for j in range(101) if j!=hub}
    tree=optimal_tree(rows[q])
    graph=star|tree
    print("最坏绕数、锐因子、直达运输树与星图的边数",q+1,target.str(17),len(graph))
    # 在其余各行先用运输树的现有直边抵消异号需求。
    # 剩余流沿至第六枚的最短路发送；据此逐根尝试删除星图辐边。
    routes=[direct_then_hub(row,hub,tree) for i,row in enumerate(rows) if i!=q]
    def worst_cost(edges):
        dist=hub_distances(edges,hub,101)
        return max((c+sum(d*t for d,t in zip(dist,r,strict=True))).upper()
                   for c,r in routes)

    assert worst_cost(graph)<target
    for edge in sorted(star-tree):
        trial=graph-{edge}
        if worst_cost(trial)<target:
            graph=trial
    remaining=worst_cost(graph)
    assert remaining<target
    print("其余104行的最大成本包围",remaining.str(14))
    payload="".join(f"{i+1} {j+1}\n" for i,j in sorted(graph))
    digest=sha256(payload.encode("ascii")).hexdigest()
    assert len(graph)==123 and remaining<arb("8.88544")
    assert digest=="e588e7b5c25571752b8df790e82c471e65f4d1ca4ece32b14337dc47724dd02b"
    print("最终边数及逐行名单 SHA-256",len(graph),digest)
    if show_edges:
        print(payload,end="")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("replay",nargs="?",type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument("--edges",action="store_true",help="按升序输出一对一行的全部编号")
    args=parser.parse_args()
    verify(args.replay,args.edges)
