"""七壳弱差可由八站消去时，极点必须保持不动。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    helpers=lib("seven-weak-seven-actual-sites")
    adj,restrict,norm,small=(helpers[k] for k in ("adj","restrict","norm","smaller_eigenvalue"))
    ctx.prec=256
    graph=lib("three-ray-chiral-boundary")
    frame=lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])[0][0]
    orbits=[[i for i in range(19) if abs(frame[i,k])>arb(".5")] for k in range(7)]
    assert orbits[-1]==[14] and sorted(sum(orbits,[]))==list(range(19))
    f,b,t=levels
    x=fifth*columns*acb_mat([[(b-f).sqrt(),0],[0,(t-f).sqrt()]])
    source=adj(frame)*x
    rows=[]
    for number in range(2_306_888_100_092,2_306_888_100_099):
        _,basis=lib("adjacent-pvm-transport")["spectrum"](graph,6,number)
        rows.append((adj(basis)*frame).tolist())
    first=[row[0] for row in rows]
    radius=(arb(3000)/arb("1.22089e19")-f).sqrt()
    choices=lib("seven-cyclic-weak-horizon")["CANDIDATES"]
    def supports(active):
        for sizes in product(*(range(1,len(orbits[k])+1) for k in active)):
            if sum(sizes)==8:
                scale=[(arb(len(orbits[k]))/s).sqrt() for k,s in zip(active,sizes,strict=True)]
                for groups in product(*(combinations(orbits[k],s) for k,s in zip(active,sizes,strict=True))):
                    local=restrict(x,[i for g in groups for i in g])
                    base=acb_mat([[sum(x[i,h] for i in g)/arb(len(g)).sqrt()-a*source[k,h] for h in (0,1)]
                                  for g,k,a in zip(groups,active,scale,strict=True)])
                    yield groups,scale,adj(local)*local,base
    for absent,drop in [(3,3),(4,2)]:
        active=[k for k in range(7) if k!=absent]
        lines=[row for j,row in enumerate(first) if j!=drop]
        inv=acb_mat([[row[k] for k in active] for row in lines]).inv()
        mid=-inv*acb_mat([[row[absent]] for row in lines])*restrict(source,[absent])
        rad=[radius*sum(abs(z) for z in row) for row in inv.tolist()]
        count,gap=0,arb(10)
        for groups,s,gram,base in supports(active):
            y=base+acb_mat([[a*z for z in row] for a,row in zip(s,mid.tolist(),strict=True)])
            err=sum((2*norm(row)+a*d)*a*d for row,a,d in zip(y.tolist(),s,rad,strict=True))
            v=(-small(gram-adj(y)*y)-err).lower()
            assert v.is_finite() and v>arb(".02"),(groups,v)
            count+=1;gap=min(gap,v)
        assert count==2835
        print("固定非极轨道、排除的八站支撑数、违反余量",orbits[absent],count,gap.str(12),flush=True)
    caps=[arb(0) for _ in range(7)]
    mixed=0
    for ids in product(*choices):
        if any(ids):
            inv=acb_mat([r[i] for r,i in zip(rows,ids,strict=True)]).inv()
            for k,row in enumerate(inv.tolist()):
                z=radius*sum(abs(a) for a in row)
                assert z.is_finite()
                caps[k]=max(caps[k],z.upper())
            mixed+=1
    assert mixed==971
    lines=[row for j,row in enumerate(first) if j!=2]
    inv=acb_mat([row[:6] for row in lines]).inv()
    coeff=(-inv*acb_mat([[row[6]] for row in lines])).entries()+[arb(1)]
    rad=[radius*sum(abs(z) for z in row) for row in inv.tolist()]+[arb(0)]
    # test only fixed directions; each valid directional test individually excludes the support.
    tests=[acb_mat([[1],[0]]),acb_mat([[0],[1]])]
    for size in (arb(".5"),arb(1),arb(2)):
        for k in range(8):
            angle=arb.pi()*k/4
            tests.append(acb_mat([[1],[size*acb(angle.cos(),angle.sin())]])/(1+size**2).sqrt())
    count,first_gap,mixed_gap=0,arb(10),arb(10)
    for groups,s,gram,base in supports(list(range(7))):
        a=acb_mat([[t*v] for t,v in zip(s,coeff,strict=True)])
        projected=base-a*adj(a)*base/(adj(a)*a)[0,0]
        d=norm([t*z for t,z in zip(s,rad,strict=True)])
        v=(-small(gram-adj(projected)*projected)-(2*norm(projected.entries())+d)*d).lower()
        assert v.is_finite() and v>arb(".06"),(groups,v)
        first_gap=min(first_gap,v)
        delta=adj(base)*base-gram
        # ||y+e||²>=||y||²-2||y||·||e|| gives a sharper one-sided test.
        best=arb(-10)
        for t in tests:
            err=2*sum(abs(z)*s0*d for z,s0,d in zip((base*t).entries(),s,caps,strict=True))
            v=((adj(t)*delta*t)[0,0].real-err).lower()
            assert v.is_finite(),groups
            best=max(best,v)
        assert best>arb(".04"),(groups,best)
        mixed_gap=min(mixed_gap,best);count+=1
    assert count==4374
    print("触及七轨的八站支撑数、全第一线余量、混合选线余量",count,first_gap.str(12),mixed_gap.str(12),flush=True)


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))
