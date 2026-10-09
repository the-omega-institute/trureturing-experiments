import numpy as np
def tree(s):
    l=s*(s+1); k=l+1; N=1+k*(l+1)
    L=np.zeros((N,N)); idx=1
    hubs=[]
    for j in range(k):
        h=idx; idx+=1; hubs.append(h); L[0,h]=L[h,0]=-1
        for t in range(l):
            f=idx; idx+=1; L[h,f]=L[f,h]=-1
    for i in range(N): L[i,i]=-L[i].sum()
    return L,hubs,N
def peak(L,u,v,tau,tol=1e-7):
    w,V=np.linalg.eigh(L)
    # group eigenvalues
    groups=[]; 
    for i,x in enumerate(w):
        if groups and abs(x-groups[-1][0])<tol: groups[-1][1].append(i)
        else: groups.append([x,[i]])
    B=0; Uvu=0
    for th,ids in groups:
        Evu=sum(V[v,i]*V[u,i] for i in ids)
        B+=abs(Evu); Uvu+=np.exp(1j*tau*th)*Evu
    return B,abs(Uvu),len(groups)
for s in [1,2,3,4,5]:
    L,hubs,N=tree(s)
    B,U,g=peak(L,0,hubs[0],np.pi)
    a=s*s+1;b=(s+1)**2+1
    pred=2*(s+1)**2/(((s+1)**2+1)*(2*s+1))
    print(f"s={s} N={N} distinct eigs={g} B={B:.12f} |U(pi)|={U:.12f} gap={B-U:.2e} predicted={pred:.12f}")
# also check other tau for even s: best over a grid
for s in [2,4]:
    L,hubs,N=tree(s); best=max(peak(L,0,hubs[0],t)[1] for t in np.linspace(0,2*np.pi,2001)); print('even s',s,'max|U| on grid',best, 'B',peak(L,0,hubs[0],0)[0])
