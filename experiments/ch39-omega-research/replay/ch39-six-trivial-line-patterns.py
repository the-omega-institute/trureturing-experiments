"""保持七维平凡字符的六壳共同解的必要低线个数，从第五壳源重建。"""
from itertools import product
from pathlib import Path
from runpy import run_path
from flint import arb, acb_mat

def certify(root,levels,blocks,model):
    sieve=run_path(str(root/'ch39-six-cyclic-line-sieve.py'))
    candidates=run_path(str(root/'ch39-seven-cyclic-weak-horizon.py'))['CANDIDATES'][:6]
    tzero=sieve['squared'](blocks[0])
    T=sum(sieve['squared'](b) for b in blocks[1:])
    survivors=[]
    for ids in product(*candidates):
        targets=sieve['bounds'](model,levels,tzero,T,ids)
        for s,(u,v) in zip(model,targets,strict=True):
            w=[-s['w'][i] for i in range(19) if i not in s['triv']]
            assert all(x>0 for x in w)
            # 二者的松量覆盖十二条非平凡线的有质量热核误差；
            # bounds 自带的余量继续保留，宁可扩大可行区。
            u+=arb(0,(arb('.04')/(s['ck']*T**2)).upper())
            v+=arb(0,(arb('.045')/(s['cw']*T**4)).upper())
            if max(arb(0),u.lower())**2 > sum(w)*v.upper():
                break
            if any(b.upper() < sum(x**(-arb(1)/(q-1)) for x in w)**(1-q)
                   for q,b in ((2,u),(4,v))):
                break
        else:
            survivors.append(ids)
    patterns=list(product((0,15),repeat=6))
    expected={i for i in patterns if i.count(15) in (2,5)}|{(15,15,0,0,0,15)}
    assert set(survivors)==expected and len(survivors)==22
    eps=(arb(3000)/arb("1.22089e19")-levels[0]).sqrt()
    eye=acb_mat([[int(i==j) for j in range(7)] for i in range(7)])
    bound=arb("1e19")
    for ids in survivors:
        B=acb_mat([sieve["row"](shell["zero"],i).tolist()[0] for shell,i in zip(model,ids,strict=True)])
        dual=sieve["adj"](B)*(B*sieve["adj"](B)).inv()
        null=eye-dual*B
        k=max(range(7),key=lambda j:null[j,j].real.lower())
        n=acb_mat([[null[j,k]] for j in range(7)])
        n=n/sieve["squared"](n).sqrt()
        cap=eps*sum(sum(abs(dual[k,j])**2 for k in range(7)).sqrt() for j in range(6))
        floors=[]
        for shell,i in zip(model,ids,strict=True):
            row=sieve["row"](shell["zero"],15-i)
            err=eps*sum(abs(z) for z in (row*dual).entries())
            lower=max(arb(0),(abs((row*n)[0,0])*(tzero-cap**2).sqrt()-err).lower())
            floors.append(((levels[0]+lower**2)*arb("1.22089e19")).lower())
        for j in range(1,6):
            if ids[j]!=ids[j-1]:
                bound=min(bound,floors[j-1],floors[j])
    assert bound>arb("4.2e7")
    print("每处换线两条较高截止的下界 GeV",bound.str(12),flush=True)
    print('六壳候选总数、排除数、余下组合中第16线的出现数分布',324,324-len(survivors),
          {k:sum(i.count(15)==k for i in survivors) for k in (2,3,5)},flush=True)

def experiment(root,fifth,columns,levels):
    sieve=run_path(str(root/'ch39-six-cyclic-line-sieve.py'))
    candidates=run_path(str(root/'ch39-seven-cyclic-weak-horizon.py'))['CANDIDATES'][:6]
    blocks,(_,_,model)=sieve['prepare_shells'](root,fifth,columns,levels,candidates,arb(0))
    fields=run_path(str(root/'ch39-three-index-free-atoms.py'))['broken_fields'](
        root,range(2_306_888_100_092,2_306_888_100_098))[2]
    for s,field in zip(model,fields,strict=True):
        ek,ew=sieve['thermal_defect'](field,levels[-1]);w=-sum(t for i,t in enumerate(s['w']) if i not in s['triv'])
        assert ek*w<arb('.04') and ew*w<arb('.045')
    certify(root,levels,blocks,model)

if __name__=='__main__':
    root=Path(__file__).resolve().parent
    source=run_path(str(root/'ch39-fifth-full-positive-joint.py'))
    experiment(root,*source['verify'](root,return_operator='spectral'))
