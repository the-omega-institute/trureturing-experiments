"""第八壳全部大权层中只许第十层取高截止时的引力障碍。"""
from pathlib import Path
from runpy import run_path
from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number = 2_306_888_100_093
    lib('sixth-third-pivot-gravity-cutoff')['verify'](root, number=number,maximum='.811612')
    _, source = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root, refine_node=2_306_888_100_086)[2]
    _, _, pivots, _ = next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root, number) if row[0]==6)
    _, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root,[number])
    ctx.prec = 192
    c=[v.real for v in pivots[2]]
    fixed=(0,15,16,17,18)
    selected=(9,6,14)
    others=[i for i in range(19) if i not in (*fixed,*selected)]
    assert all(c[i]<0 for i in others)
    cap=arb(".75553")
    assert source[-1].sqrt()<cap
    trace=sum(v.sqrt() for v in source).upper()
    bounds=[".001" if i in fixed else ".75553" for i in range(19)]
    m2,m4=lib("fifth-seventh-cutoff-separation")["moment_errors"](field,heat,c,bounds)
    tv=sum(abs(c[i]) for i in fixed)
    k=m2.upper()+tv*arb("1e-6")
    w=-m4.lower()+tv*arb("1e-12")
    largest=max((-c[i]).upper() for i in others)
    def impossible(lo,hi):
        rem=trace-sum(lo)
        if rem<0:return True
        a=c[9]*lo[0]**2+c[6]*hi[1]**2+c[14]*hi[2]**2-k
        b=c[9]*hi[0]**4+c[6]*lo[1]**4+c[14]*lo[2]**4+w
        return b<0 or a>largest*rem**2 or (a>0 and a**3>b*(largest*rem)**2)
    stack=[((arb(0),)*3,(cap,)*3)]
    count=0
    while stack:
        lo,hi=stack.pop();count+=1
        if impossible(lo,hi):continue
        i=max(range(3),key=lambda j:float(hi[j]-lo[j]));assert hi[i]-lo[i]>arb(".0001"),(lo,hi)
        mid=(lo[i]+hi[i])/2
        stack.extend(((lo,hi[:i]+(mid,)+hi[i+1:]),(lo[:i]+(mid,)+lo[i+1:],hi)))
    print("只放开第十个大权层的三维引力排除盒数",count,flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
