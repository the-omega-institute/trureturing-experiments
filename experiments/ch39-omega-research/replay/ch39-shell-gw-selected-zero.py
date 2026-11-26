"""Global zero selection for all positively matched Wilson weights and the resulting chiral character."""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat,arb, ctx


folder=Path(__file__).resolve().parent
s=run_path(str(folder/"ch39-shell-gw-charge-conjugation.py"))
model=s["s"]
ctx.prec=96
angles=[(1-arb(d)/26).acos() for d in (3,5,7)]
f,g=angles[1]-angles[0],angles[2]-angles[0]
lam=arb(1)/51**5


def loss(h):
    return 1+(f*h).cos()*(g*h).cos()


def certify_height():
    h=acb.zeta_zero(43).imag
    u=loss(h)
    ends=[lam/30,30*lam]
    n=int((u/ends[0]+h*h).sqrt().zeta_nzeros().unique_fmpz())
    for start in range(1,n+1,128):
        for idx,z in enumerate(acb.zeta_zeros(start,min(128,n+1-start)),start=start):
            if idx!=43:
                dh=z.imag*z.imag-h*h
                du=loss(z.imag)-u
                assert all(du+t*dh>0 for t in ends)
    print("unique selected positive zero 43, both ends and tail, initial indices tested",n,
          "geometric excess",u.str(15))


def action_kernel(sign,p,q):
    base=s["scalar_kernel"](sign,False)
    matrix=model["equivariant_matrix"](base,(-p if sign<0 else p),q)
    for i in range(36):
        matrix[i][i]-=arb(1)/4 if i<18 else -arb(1)/4
    return base,matrix


def character_blocks(base,matrix):
    orbit=[2,3,4,5,0,1]
    rotate=[orbit[i%6]+6*(i//6) for i in range(12)]
    phases={0:acb(1)}
    queue=[0]
    for j in queue:
        for i in range(12):
            if i not in phases and abs(base[i][j])>arb(".01"):
                phases[i]=phases[j]*base[rotate[i]][rotate[j]]/base[i][j]
                queue.append(i)
    assert len(phases)==12
    for i in range(12):
        for j in range(12):
            assert (base[rotate[i]][rotate[j]]-phases[i]*base[i][j]*phases[j].conjugate()).contains(0)
    cube=phases[0]*phases[rotate[0]]*phases[rotate[rotate[0]]]
    fix=acb(0,-cube.arg()/3).exp()
    # Geometric and family coordinates undergo the same cycle.
    rot=[3*rotate[i//3]+(i%3+1)%3 for i in range(36)]
    coeff=[fix*phases[i//3] for i in range(36)]
    seeds=[i for i in range(36) if i==min(i,rot[i],rot[rot[i]])]
    assert len(seeds)==12
    sub=[]
    for k in range(3):
        factor=acb(0,-2*arb.pi()*k/3).exp()
        v=acb_mat(36,12)
        for j,start in enumerate(seeds):
            x=start
            c=acb(1)/arb(3).sqrt()
            for step in range(3):
                v[x,j]=c
                c*=coeff[x]*factor
                x=rot[x]
        comp=v.conjugate().transpose()*acb_mat(matrix)*v
        sub.append([[comp[i,j] for j in range(12)] for i in range(12)])
    return sub


if __name__=="__main__":
    certify_height()
    ctx.prec=240
    h=acb.zeta_zero(43).imag
    for sign in (1,-1,0):
        base,k=action_kernel(sign,f*h,g*h)
        gap=arb(1)/100
        assert model["model"]["count"](k,-gap)==model["model"]["count"](k,gap)==18-3*sign
        model["calc"]["opposite_modes"](k,3*sign,gap)
        indices=[6-model["model"]["count"](sub,arb(1)/1001) for sub in character_blocks(base,k)]
        assert indices==[sign]*3
        print("charge",sign,"rotation character indices",indices,"no opposite chirality pairs")
