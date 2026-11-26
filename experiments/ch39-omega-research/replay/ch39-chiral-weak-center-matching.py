"""11661 的三代边界和中心奇偶质量在两种一圈定义下的交汇实例。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent
chiral = run_path(str(ROOT / "ch39-chiral-common-gap-crossing.py"))
dot = lambda a,b: sum(x*y for x,y in zip(a,b))


def spectrum(model, full):
    read = chiral["source"]["spectral_supplier"](model)
    rows = read(acb.zeta_zero(11661).imag)
    layers = {q: chiral["ordered"]([h for h,z in rows if z==q])
              for q in (0,*model["charges"])}
    out = []
    for q, w in zip(model["charges"],model["weyl"]):
        assert layers[q][2]<9<layers[q][3]
        ratio = 1 if q%2 else 8
        out.extend((ratio*h,[2*x for x in w],"F") for h in layers[q][3:])
    out.extend((8*h,[0,-7,-arb(21)/2],"V") for h in layers[0])
    out.extend((h,model["higgs"],"S") for h in layers[3][0 if full else 3:])
    return out


def verify():
    model,_,_,initial,_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    beta = model["beta"]
    d = (beta[1]-beta[2], beta[2]-beta[0], beta[0]-beta[1])
    pi = arb.pi()
    base = [2*pi*v for v in initial]
    def electroweak(rhs):
        return (rhs[0]-rhs[1])/(beta[0]-beta[1])
    for full in (True,False):
        rows = spectrum(model,full)
        ctx.prec = 192
        def report(label,t,L,rhs):
            mu = model["mz"]*L.exp()
            kappa = mu*(-t).exp()
            assert arb("1e13")<mu<arb("1e17")
            assert all(kappa*h>model["mz"] for h,_,_ in rows)
            common = (rhs[0]-beta[0]*L)/(2*pi)
            assert common>0
            print("十九标量" if full else "十六标量",label,
                  "交点/GeV、质量单位/GeV、共同逆耦合",mu.str(11),kappa.str(11),common.str(10))
        sample = arb("4.437" if full else "4.528")
        active = [(h,w) for h,w,_ in rows if h.log()<sample]
        assert len(active)==(50 if full else 47)
        sums = [sum(w[i] for _,w in active) for i in range(3)]
        rhs = [base[i]+sum(h.log()*w[i] for h,w in active) for i in range(3)]
        t = dot(d,rhs)/dot(d,sums)
        assert all(h.log()<t if h.log()<sample else h.log()>t for h,_,_ in rows)
        rhs = [a-v*t for a,v in zip(rhs,sums)]
        report("尖阈值",t,electroweak(rhs),rhs)

        def smooth(t,L):
            result = list(base)
            high = (2*t).exp()
            low = (2*(t-L)).exp()
            for h,w,kind in rows:
                fn = model["scalar"] if kind=="S" else lambda m,p:model["sheet"](m,p,kind=="V")
                minimum = fn(h.upper()**2,high.lower())-fn(h.lower()**2,low.upper())
                maximum = fn(h.lower()**2,high.upper())-fn(h.upper()**2,low.lower())
                value = minimum.union(maximum)
                result = [a-x*value for a,x in zip(result,w)]
            return result
        lbox = arb(27).union(arb("27.3"))
        left,right = (arb("5.14"),arb("5.17")) if full else (arb("5.25"),arb("5.28"))
        sign = 1 if dot(d,smooth(left,lbox))<0 else -1
        assert sign*dot(d,smooth(left,lbox))<0<sign*dot(d,smooth(right,lbox))
        for _ in range(28):
            mid = (left+right)/2
            val = sign*dot(d,smooth(mid,lbox))
            if val<0:
                left=mid
            elif val>0:
                right=mid
            else:
                break
        assert right-left<arb("1e-8")
        t = left.union(right)
        rhs = smooth(t,lbox)
        L = electroweak(rhs)
        # 两个 t 面的投影异号，两个 L 面由此严格包含保证异号。
        assert lbox.contains_interior(L)
        report("有限减法",t,L,rhs)


if __name__=="__main__":
    verify()
