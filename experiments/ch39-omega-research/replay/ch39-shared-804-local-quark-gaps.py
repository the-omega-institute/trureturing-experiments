"""共同 804 钟的两段局部排除与总跨度中值见证。"""
from pathlib import Path
from runpy import run_path
from flint import acb,acb_mat,arb,ctx
root=Path(__file__).resolve().parent
def lib(name):
    return run_path(str(root/('ch39-'+name+'.py')))
graph,pump=lib('three-ray-chiral-boundary'),lib('native-light-projector-pump')
ctx.prec=160
adj=pump['adjoint']
W=lib('moving-preparation-quark-shape')['orthonormal_sources'](graph,adj)
sv=lib('first-zero-relative-quark-spectrum')['singular_values']
models=[]
for q in (-4,2):
    L,(A,B)=pump['clock_matrices'](graph,q)
    models.append([(adj(Q)*L*Q,adj(Q)*A*Q,adj(Q)*B*Q,adj(Q)*W) for Q in pump['rotation_bases'](graph)])
def readout(k,a,b,filter=True):
    a,b=acb(0,a).exp(),acb(0,b).exp()
    rows=[]
    min_gap=100
    for L,A,B,W in models[k]:
        C,D=A*a,B*b
        H=L+C+adj(C)+D+adj(D)
        v,gap=pump['first_line'](H)
        if gap<min_gap:min_gap=gap
        column=acb_mat([[z] for z in v])
        row=adj(column)
        energy=(row*H*column)[0,0].real
        rows.append(((row*W)/(1+(energy/9)**512) if filter else row*W).entries())
    return acb_mat(rows),min_gap
T=acb_mat(lib('shell-clock-ckm-blocks')['standard_matrix'](*(arb(z).sin() for z in ('.22704','.04275','.00377')),arb('1.139')))
for i,u in enumerate(('.00000615','.00311','.8616')):
    for j,d in enumerate(('.0000135','.000268','.01401')):
        T[i,j]*=arb(u)/arb(d)
t=sv(T)
goals=t[1]/t[0],t[2]/t[1]

def projector_center(k,a,b,ra,rb):
    center,gap=readout(k,a,b,False)
    displacement=arb(2).sqrt()*(ra+rb)+arb('6e-8')
    assert gap > 2*displacement
    error=arb(2).sqrt()*displacement/(gap-displacement)+arb('2e-8')
    return center,error
for turn in (1,2):
    beta=turn*arb.pi()-arb('.87833265')/2
    rb,ra=arb('.0005'),arb('.0012')
    U,eu=projector_center(0,arb('5.60000358'),-2*beta,ra,2*rb)
    cells=[(i,64) for i in range(64)]
    n=0
    worst=arb(0)
    least=arb(1)
    finest=0
    while cells:
        i,N=cells.pop()
        D,ed=projector_center(1,arb.pi()*(2*i+1)/N,beta,arb.pi()/N,rb)
        sd=sv(D)[0].lower()
        if sd>ed:
            j=sv(U*D.inv())
            lower=(j[1]-eu/sd)/(1+ed/sd)
            upper=(j[2]+eu/sd)/(1-ed/sd)
            if lower>0:
                bound=(arb(9)/7)**2*upper/lower
                if bound<(18 if turn==1 else 32):
                    n+=1
                    finest=max(finest,N)
                    worst=max(worst,bound.upper())
                    least=min(least,(sd-ed).lower())
                    continue
        assert N<8192,(turn,N,i)
        cells.extend(((2*i,2*N),(2*i+1,2*N)))
    assert worst < arb('17.878' if turn==1 else '31.821')
    print('lift',turn,'cells',n,'finest',finest,'bound',worst.str(15),'least',least.str(12))

beta=arb.pi()-arb('.87833265')/2
D,_=readout(1,arb(0),beta)
mid = (arb('5.60024461948') + arb('5.60024461961'))/2
center, error = projector_center(0,mid,-2*beta,arb('1e-10'),arb(0))
assert sv(center)[2] > error
for decimal in ('5.60024461948','5.60024461961'):
    U,_=readout(0,arb(decimal),-2*beta)
    a,b,c=sv(U*D.inv())
    diff=c/a-goals[0]*goals[1]
    assert (diff if decimal.endswith("48") else -diff) > 0
    print(decimal,(c/a-goals[0]*goals[1]).str(12),(c/b).str(12))
