"""Eliminate both extreme common scales by the first heat-moment equation."""
from pathlib import Path
import runpy

m = runpy.run_path(str(Path(__file__).with_name("ch39-fixed-inertia-scale.py")))
q = m["q"]
minimum,maximum = q(390),q(429)
def exclude(a,b,sign,steps):
    margin = None
    for i in range(steps):
        rr = a+(b-a)*q(i).union(q(i+1))/steps
        for j in range(6):
            v = m["v0"]+(m["v1"]-m["v0"])*q(j).union(q(j+1))/6
            for k in range(4):
                u = m["u0"]+(m["u1"]-m["u0"])*q(k).union(q(k+1))/4
                _,res = m["chart"](rr,v,u)
                gap = (sign*res).lower()
                assert gap>0
                margin = gap if margin is None else min(margin,gap)
    print("excluded r interval and moment margin",a,b,margin.str(12))
if __name__=="__main__":
    exclude(m["r0"],minimum,1,340)
    exclude(maximum,m["r1"],-1,20)
