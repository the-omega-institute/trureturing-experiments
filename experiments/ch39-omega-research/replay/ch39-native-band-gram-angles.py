"""五荷原壳最低三维谱带与四种算术三源空间的主夹角、交维数与投影秩。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx


def adj(m):
    return m.conjugate().transpose()


def restricted(lap,sites):
    vals,v=lap.eig(right=True,algorithm="rump")
    order=sorted(range(19),key=lambda i:vals[i].real.mid())
    assert all(z.imag.contains(0) for z in vals)
    assert all(vals[i].real<vals[j].real for i,j in zip(order[:3],order[1:4],strict=True))
    cols=[]
    for j in order[:3]:
        n=sum(abs(v[i,j])**2 for i in range(19)).sqrt()
        cols.append([v[i,j]/n for i in sites])
    return acb_mat(cols).transpose(),[vals[j].real for j in order[:4]]


def angles(block,n):
    cross=block-n*(adj(n)*block)
    gram=adj(cross)*cross
    ev=sorted((v.real for v in gram.eig(algorithm="rump")),key=lambda v:v.mid(),reverse=True)
    return ev


def blind_matrix(ray,cyclic,q,h):
    base,clocks=cyclic["clock_matrices"](ray,q)
    pole=ray["soft"]["pole"]
    delta=arb(23)/25*arb(2)**-20/(h*h+ray["soft"]["point"]["study"]["offset"])
    for site in ray["soft"]["star"]["orbit"]:
        base[pole,pole]-=delta
        base[site,site]-=delta
        base[pole,site]+=delta
        base[site,pole]+=delta
    bases=list(cyclic["rotation_bases"](ray))
    for index in (1,2):
        basis=bases[index]
        j,=[i for i in range(basis.ncols()) if not basis[6,i].is_zero()]
        other=[i for i in range(basis.ncols()) if i!=j]
        fixed,first,second=[adj(basis)*a*basis for a in (base,*clocks)]
        core=acb_mat([[fixed[i,k] for k in other] for i in other])
        pair=core.eig(right=True,algorithm="rump")
        vals,vec=pair
        ix=min(range(len(vals)),key=lambda i:vals[i].real.mid())
        assert all(vals[ix].real<vals[k].real for k in range(len(vals)) if k!=ix)
        co=[sum(vec[k,ix].conjugate()*a[i,j] for k,i in enumerate(other))
            for a in (fixed,first,second)]
        a,b,c=co
        cs=(abs(c)**2-abs(a)**2-abs(b)**2)/(2*abs(a)*abs(b))
        if cs>1 or cs<-1:
            continue
        assert -1<cs<1
        u=(a/abs(a))*(b.conjugate()/abs(b))*acb(cs,(1-cs**2).sqrt())
        v=-(a+b*u)/c
        assert abs(u).contains(1) and abs(v).contains(1)
        off=clocks[0]*u+clocks[1]*v
        actual=base+off+adj(off)
        return actual,vals[ix].real,index
    raise ValueError("两非平凡字符中未找到可闭合的内轨连接三角形")


def verify(root):
    ray=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    triple=run_path(str(root/"ch39-four-gram-triplet-ckm.py"))
    cyclic=run_path(str(root/"ch39-native-light-projector-pump.py"))
    ctx.prec=256
    height=acb.zeta_zero(2_306_888_100_085).imag
    gram=native["gram_cholesky"](root)
    normals=[acb_mat([[v] for v in triple["normal"](gram,o)]) for o in range(4)]
    sites=[ray["soft"]["pole"],*ray["boundary_sites"]()]
    allvalues=[]
    for q in (1,-4,2,-3,6):
        block,_=restricted(ray["laplacian"](q,height),sites)
        for o,n in enumerate(normals):
            es=angles(block,n)
            assert 0<es[2]<es[1]<es[0]<1
            vals=[x.sqrt().acos()*180/arb.pi() for x in es]
            print(q,o,"主角（度）",*[v.str(8) for v in vals],flush=True)
            allvalues+=es
        off,e,index=blind_matrix(ray,cyclic,q,height)
        block,lowest=restricted(off,sites)
        assert e<lowest[3] and sum(z.overlaps(e) for z in lowest[:3])==1
        for n in normals:
            ev=angles(block,n)
            assert ev[1]>0 and ev[2].contains(0)
        print("独立有效钟相位，全部三源秩二的字符",q,index,flush=True)
    assert min(z.lower() for z in allvalues)>arb(".0063743")**2
    assert max(z.upper() for z in allvalues)<arb(".411014")**2
    print("实际零点的全范围奇异值",min(z.lower() for z in allvalues).sqrt().str(10),
          max(z.upper() for z in allvalues).sqrt().str(10),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
