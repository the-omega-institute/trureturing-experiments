"""在全部弱差可达的三线支撑上，用条件 Gram 矩阵外包各层振幅和两个热迹。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    _,_,incoming,readout,(low,h,k),weights,(_,heat,fields) = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    shell = lib("three-ray-chiral-boundary")
    mid = lib("adjacent-pvm-transport")["spectrum"](shell,6,lib("rank-two-generator-up-quark")["N"]+3)[1]
    triple = lib("fourth-shell-three-line-weak")
    cases = triple["verify_data"](root,mid,incoming,readout,(low,h,k),weights,trivial_counts=(1,2,3))
    adj = triple["adj"]
    def norm(row):
        sizes = [abs(z) for z in row.entries()]
        lower = sum(max(arb(0),s.lower())**2 for s in sizes).sqrt()
        upper = sum(s.upper()**2 for s in sizes).sqrt()
        return lower.union(upper)
    target = adj(readout)
    high = incoming*acb_mat([[(h-low).sqrt(),0],[0,(k-low).sqrt()]])
    original = readout*high
    limit = 16*arb("91.1876")/arb("1.22089e19")
    small = (limit-low).sqrt()
    max_f = arb(".76")
    assert k<max_f
    ak = sum(a for _,a,_ in fields)/(2*arb.pi())
    aw = -sum(g for g,_,_ in fields)/(8*arb.pi())
    xmax = max_f**2
    err_k = xmax*sum(abs(a)*(1-heat(2,z/xmax)) for _,a,z in fields)/(2*arb.pi())
    err_w = xmax**2*sum(abs(g)*(arb(1)/2-heat(3,z/xmax)) for g,_,z in fields)/(4*arb.pi())
    assert err_k<arb(".000131") and err_w<arb(".000021")
    tv = sum(map(abs,weights))

    def band(intervals,slope):
        def value(square):
            return square*(ak+slope*aw*square)
        total = arb(0)
        for weight,(lower,upper) in zip(weights,intervals,strict=True):
            l,r = lower**2,upper**2
            col = value(l).union(value(r))
            if slope>0:
                crit = -ak/(2*slope*aw)
                if not (crit<l or r<crit):
                    col = col.union(value(crit))
            total += weight*col
        return total+arb(0,(tv*(err_k+abs(slope)*err_w)).upper())

    slopes = [arb(s) for s in ("0","-50","1.3","1.4",".3",".4",".5")]
    for ntriv,(_,allowed) in cases.items():
        open_cases = []
        power18,power19 = [],[]
        for (left,right),_,_ in allowed:
            e = acb_mat([[*[mid[i,j-1] for j in left],*[target[i,j-1] for j in right]] for i in range(19)])
            a = adj(e)*high
            co = readout*e
            fixed = original-co*a
            b = acb_mat([[fixed[15,0],fixed[15,1]]])
            u = acb_mat([[co[15,j] for j in range(3)]])
            r = norm(u)
            v = adj(u)*(1/r)
            row0 = -b*(1/r)
            gram = adj(a)*a-adj(row0)*row0
            delta = small/r
            err = delta*(2*norm(row0)+delta)
            diff = ((gram[0,0].real-gram[1,1].real)**2+4*abs(gram[0,1])**2).sqrt()
            eig0 = ((gram.trace().real-diff)/2-err).lower()
            eig1 = ((gram.trace().real+diff)/2+err).upper()
            s0,s1 = max(arb(0),eig0).sqrt(),eig1.sqrt()
            intervals = []
            for i in range(19):
                if i==15:
                    intervals.append((low.lower(),limit.upper()))
                    continue
                g = acb_mat([[co[i,j] for j in range(3)]])
                proj = (g*v)[0,0]
                perp = norm(g-u*(proj/r))
                b0 = norm(acb_mat([[fixed[i,0],fixed[i,1]]])+row0*proj)
                jitter = abs(proj)*delta
                bottom = max(arb(0),(b0-jitter-perp*s1).lower(),(perp*s0-b0-jitter).lower())
                top = (b0+jitter+perp*s1).upper()
                intervals.append(((low+bottom**2).lower(),min(k.upper(),(low+top**2).upper())))
            assert all(l.is_finite() and r.is_finite() for l,r in intervals)
            power18.append(intervals[17][1])
            power19.append(intervals[18][1])
            if all(band(intervals,t).overlaps(1+t*arb("1e-122")) for t in slopes):
                open_cases.append((left,right))
        print("平凡线、已排除、热迹带容许",ntriv,len(allowed)-len(open_cases),len(open_cases),open_cases,flush=True)
        print("该类所有弱差零点第18、19层振幅上界",max(power18).str(13),max(power19).str(13),flush=True)
        assert (len(allowed)-len(open_cases))=={1:81,2:81,3:12}[ntriv]
    print("单层牛顿和真空列的绝对近似误差",err_k.str(12),err_w.str(12),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
