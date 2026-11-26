"""由右弱壳的四边环相位构造真实的三代 CP 不变量。"""
from collections import Counter
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import fmpq_mat
from sympy import I, Matrix, conjugate, simplify, sqrt

shell = run_path(str(Path(__file__).with_name("ch39-shell-vector-gap.py")))


def cycles(lap, edges):
    nodes = [[] for _ in range(lap.nrows())]
    for a,b in edges:
        nodes[a].append(b)
        nodes[b].append(a)
    result = []

    def extend(path):
        for other in nodes[path[-1]]:
            if other == path[0] and len(path) >= 3 and path[1] < path[-1]:
                result.append(path)
            elif other > path[0] and other not in path:
                extend(path+[other])

    for a in range(len(nodes)):
        extend([a])
    return nodes, result


def loop_rank(loops, edges):
    rows = []
    for path in loops:
        row = [0]*len(edges)
        for a,b in zip(path,path[1:]+path[:1],strict=True):
            row[edges.index((max(a,b),min(a,b)))] = 1 if a>b else -1
        rows.append(row)
    return fmpq_mat(rows).rank()


def verify():
    _,_,colour,right = shell["neutral"]["graph_constants"](port_builder=shell["capture"])
    for graph, signature, local_rank in ((colour,{3:3,9:1},3),(right,{4:3},2)):
        lap,edges = graph
        neighbours,loops = cycles(lap,edges)
        assert Counter(len(c) for c in loops) == signature
        short = [c for c in loops if len(c)<=4]
        assert loop_rank(short,edges)==local_rank
        assert loop_rank(loops,edges)==len(edges)-lap.nrows()+1
        if graph is right:
            hubs = [(a,b,sorted(set(neighbours[a])&set(neighbours[b])))
                    for a,b in combinations(range(lap.nrows()),2)
                    if len(set(neighbours[a])&set(neighbours[b]))==3]
            assert len(hubs)==1
            print("hubs and three intermediates",hubs[0],flush=True)

    omega=(-1+I*sqrt(3))/2
    assert simplify(omega**2+omega+1)==0
    up = Matrix.diag(1,4,9)
    for phase in (omega,conjugate(omega)):
        down = Matrix([[1,0,0],[phase,1,0],[1,1,1]])
        h = down*down.H
        commutator=up*h-h*up
        value=simplify(commutator.det())
        assert simplify(value-120*(conjugate(phase)-phase))==0
        print("determinant of flavour commutator",value,flush=True)


if __name__=="__main__":
    verify()
