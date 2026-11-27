"""在三条极短边共用的复门球面上，检验第七和第八壳同时恢复引力双矩。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def shell_data(root, fifth, columns, levels, number):
    geometry = run_path(str(root/"ch39-fifth-seventh-shortest-site-gates.py"))["geometry"]
    base, blocks, sine, basis, station, weights, heat, field = geometry(
        root, fifth, columns, levels, number=number)
    pairs = ((6,11),(7,12),(8,18))
    diagonal = sum(blocks[1:],blocks[0])
    fixed = base-diagonal
    real_z = sum(sine[1:],sine[0])
    imag_z = acb_mat([[acb(0,1)*sum(basis[b,i].conjugate()*station[a,j]
                                     +basis[a,i].conjugate()*station[b,j] for a,b in pairs)
                      for j in range(2)] for i in range(19)])
    low,heavy,top = levels
    assert top < arb(".75553")
    slopes = (sum(a for _,a,_ in field)/(2*arb.pi()),-sum(g for g,_,_ in field)/(8*arb.pi()))
    errors = run_path(str(root/"ch39-fifth-seventh-cutoff-separation.py"))["moment_errors"]
    normalized = errors(field,heat,weights,[".75553"]*19)
    return (fixed,diagonal,real_z,imag_z),weights,(low,heavy-low,top-low),slopes,normalized


def experiment(root, fifth, columns, levels):
    shells = [shell_data(root,fifth,columns,levels,n) for n in (2_306_888_100_092,2_306_888_100_093)]
    # θ/π∈[0,1]、φ/π∈[-1,1] 覆盖完整球面，极点包含在边界。
    stack = [(arb(0),arb(1),arb(-1),arb(1))]
    covered = 0
    minima = [arb(100)]*3
    thresholds = list(map(arb,(".20004",".40",".65")))
    while stack:
        a,b,l,r = stack.pop()
        sine, cosine = (arb.pi()*a.union(b)).sin_cos()
        sp,cp = (arb.pi()*l.union(r)).sin_cos()
        x,y = sine*cp,sine*sp
        residuals = []
        for (fixed,dg,rx,iy),weights,(low,v1,v2),slopes,targets in shells:
            proj = fixed+cosine*dg+x*rx+y*iy
            amps = [low+v1*abs(proj[i,0])**2+v2*abs(proj[i,1])**2 for i in range(19)]
            for power, slope, target in zip((2,4),slopes,targets,strict=True):
                error = slope*(sum(w*t**power for w,t in zip(weights,amps,strict=True))-target)
                residuals.append(error.abs_lower())
        gaps = (max(residuals[2:]),max(residuals[0],residuals[2]),max(residuals))
        if all(g > t for g,t in zip(gaps,thresholds,strict=True)):
            covered += 1
            minima = [min(m,g) for m,g in zip(minima,gaps,strict=True)]
        elif b-a >= r-l:
            assert b-a > arb("1e-8"),(a,b,l,r,gaps)
            mid = (a+b)/2
            stack.extend(((a,mid,l,r),(mid,b,l,r)))
        else:
            mid = (l+r)/2
            stack.extend(((a,b,l,mid),(a,b,mid,r)))
    print("闭盒数、第八壳双矩、两壳牛顿项、全部四项的最大绝对残差下界",
          covered,*[v.str(16) for v in minima],flush=True)
    return minima


def verify(root):
    joint = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*joint["verify"](root,return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
