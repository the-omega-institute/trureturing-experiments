"""对边类型的两面有理读数、两种螺旋度的二阶归一及真实零点的三阶差。"""
from collections import defaultdict
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
import sympy as sp


ROOT = Path(__file__).parent


def verify_readout():
    a, b = sp.symbols("a b",positive=True)
    seed = (1+a*a-b*b)/2
    mean = 2*seed*a*a/(seed**2+a*a)
    assert sp.factor((a-mean)*(seed**2+a*a)-a*(seed-a)**2)==0
    for expression, expected in ((mean,a),(mean.diff(a),1),(mean.diff(b),0)):
        assert sp.simplify(expression.subs(b,1-a)-expected)==0


def cubic_coefficient():
    flat = run_path(str(ROOT/"ch39-regge-flat-cell.py"))
    edges = flat["EDGES"]
    q = sp.Matrix([sp.Rational((j-i)*(5-j+i),5) for i,j in edges])
    scale = 4*sp.sqrt(5)-2*sp.sqrt(15)
    parts = flat["canonical"]()
    triangles = defaultdict(list)
    partners = defaultdict(lambda: [[],[]])
    for order in permutations(range(4)):
        vertices = [sum(1<<k for k in order[:j]) for j in range(5)]
        dq = sp.Matrix([int(vertices[j]-vertices[i]==1) for i,j in edges])
        entries = {}
        for face, cosine, area_angle in parts:
            area = sp.sqrt(sp.Rational(1,5) if cosine==0 else sp.Rational(3,20))
            indices = [edges.index(e) for e in combinations(face,2)]
            total = sum(q[k] for k in indices)
            da = sp.expand(sum((total/2-q[k])*dq[k] for k in indices)/(8*area))
            angle = sp.expand((q.T*area_angle*dq)[0]*sp.sqrt(5)/area)
            i,j,k = face
            key = vertices[j]-vertices[i],vertices[k]-vertices[j]
            triangles[key].append((area,da,angle))
            entries[face] = key,area,da
        for face,(key,area,da) in entries.items():
            for other,(_,other_area,dg) in entries.items():
                if area != other_area:
                    index = len(set(face).intersection(other))-1
                    partners[key][index].append((da+dg)**2)
    result = 0
    for key, rows in triangles.items():
        area, derivative, _ = rows[0]
        assert all(a==area and d==derivative for a,d,_ in rows)
        other = sp.sqrt(sp.Rational(3,20)) if len(rows)==4 else sp.sqrt(sp.Rational(1,5))
        vertex, edge = partners[key]
        assert len(vertex)==len(rows) and len(edge)==4*len(rows)
        difference = sum(vertex)/len(vertex) - sum(edge)/len(edge)
        result -= 3*scale**3*(other**2/area)*difference*sum(v for _,_,v in rows)
    assert sp.simplify(result-705+1629*sp.sqrt(3)/4)==0
    return 705-1629*arb(3).sqrt()/4


def verify():
    verify_readout()
    source = run_path(str(ROOT/"ch39-regge-conserved-source.py"))
    source["verify"]()
    shell = run_path(str(ROOT/"ch39-arc-relative-determinant.py"))
    ctx.prec = 192
    weight = [(1+shell["phases"](acb.zeta_zero(n).imag)[724].real)/2 for n in (11366,11431)]
    difference = cubic_coefficient()*(weight[0]-weight[1])/(8*arb.pi())
    assert arb("-.00550956") < difference < arb("-.00550954")
    print("两枚实际零点给同一平方边长形变的实作用三阶导数差",difference.str(14))


if __name__=="__main__":
    verify()
