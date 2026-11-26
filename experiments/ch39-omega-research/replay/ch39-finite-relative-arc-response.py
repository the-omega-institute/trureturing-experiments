"""One rational Cayley rotor per measured overlap on the arc-selected boundary."""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from sympy import I, Matrix, Rational as Q, eye, zeros

folder = Path(__file__).resolve().parent
cells = run_path(str(folder/"ch39-finite-shell-relative-index.py"))
vertices, edges, _ = cells["shell_cells"]()
pairs = list(combinations(range(4),2))


def main():
    a,b = vertices[:2]
    assert a == (0,1,3,4) and b == (0,1,4,3)
    for v,w in edges:
        p,q = vertices[v],vertices[w]
        cost = sum((x-y)**2 for x,y in zip(p,q))//2
        crossed = sum(p.index(k) != q.index(k) for k in (0,4))
        assert crossed == (1 if cost==1 else 0)
    assert a.index(4)!=b.index(4)
    H = Matrix([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/2
    up,down = zeros(3),zeros(3)
    coeffs = {p:0 for p in pairs}
    for p in vertices:
        t = Q(sum(x*y for x,y in zip(p,a)),26)
        assert Q(3,13)<=t<=1
        z = ((1+I*t)/(1-I*t)).expand(complex=True)
        # Project onto the direction of (5,8), the midpoint of the two endpoint phases.
        assert 5*z.as_real_imag()[0]+8*z.as_real_imag()[1]>=8
        chi = (-1)**sum(p[i]>p[j] for i,j in pairs)
        i,j = sorted((p.index(3),p.index(4)))
        root = H[i,:]-H[j,:]
        up += z*root.T*root/16
        down += chi*z.conjugate()*root.T*root/16
        coeffs[i,j] += chi*z.conjugate()/16
    for (i,j),c in coeffs.items():
        k,l = sorted(set(range(4))-{i,j})
        bound = Q(abs((a[i]-a[j])*(a[k]-a[l])),4*26**2)
        assert (c*c.conjugate()).expand()<=bound**2
    Hu,Hd = [(Y*Y.H).applyfunc(lambda x:x.expand()) for Y in (up,down)]
    cp = ((Hu*Hd-Hd*Hu).applyfunc(lambda x:x.expand()).det()).expand()
    assert cp.as_real_imag()[0] == 0 and cp.as_real_imag()[1] != 0
    # The upper bound on the bottom/top ratio, independent of the common coupling.
    ratio = Q(1401,86160)
    assert Q(9,676)**2*89 < 64*ratio**2
    print("Anchor alone gives nonzero CP in the selected cover; up singular values in [8/sqrt(89),1] in common units.")
    print("The bottom/top ratio is at most 9 sqrt(89)/(676*8), below 0.01401/0.8616.")


if __name__=="__main__":
    main()
