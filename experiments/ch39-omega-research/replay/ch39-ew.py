from fractions import Fraction
from flint import arb, ctx

ctx.prec = 160
def q(value):
    fraction = Fraction(value)
    return arb(fraction.numerator)/fraction.denominator

def matrix(rows):
    return [[q(i) for i in row] for row in rows]

beta = list(map(q, ['41/10','-19/6',-7]))
B = matrix([['199/50','27/10','44/5'],['9/10','35/6',12],['11/10','9/2',-26]])
D = matrix([['17/10','1/2','3/2'],['3/2','3/2','1/2'],[2,2,0]])
C = matrix([['17/20','9/4',8],['1/4','9/4',8],['9/4','9/4',0]])
# (4,2,1) decomposes as a colored doublet plus a colorless doublet.

table = [
    "6.15e-6 .14e-6 5.84e-6 .13e-6 5.54e-6 .13e-6 5.07e-6 .12e-6",
    ".00311 .00005 .00295 .00005 .00280 .00005 .00257 .00005",
    ".8616 .0043 .8242 .0043 .7880 .0044 .7309 .0046",
    "1.35e-5 .02e-5 1.28e-5 .02e-5 1.22e-5 .02e-5 1.12e-5 .02e-5",
    ".000268 .000003 .000254 .000003 .000242 .000003 .000222 .000003",
    ".01401 .00009 .01321 .00009 .01246 .00009 .01130 .00008",
    ".22704 .00082 .22705 .00082 .22705 .00082 .22705 .00082",
    ".04275 .00042 .04307 .00042 .04339 .00042 .04394 .00043",
    ".00377 .00008 .00380 .00008 .00383 .00008 .00388 .00008",
    "1.139 .023 1.139 .023 1.139 .023 1.139 .023",
]


from math import ceil
loop=16*arb.pi()**2
initial_g=[1/q(g)**2 for g in ('.467453','.63811','1.0583')]
initial_y=[2*arb(row.split()[0],row.split()[1] if i in (3,4) else 0).log() for i,row in enumerate(table[:6])]
initial_y.extend(2*q(v).log() for v in ['.0000028227','.00059465','.010101'])
low_input=initial_g+initial_y
def slopes(vec):
    g=[1/v for v in vec[:3]]
    y=[v.exp() for v in vec[3:12]]
    tr=[sum(y[k:k+3]) for k in (0,3,6)]
    T=3*tr[0]+3*tr[1]+tr[2]
    gdot=[-2*beta[i]/loop-2*(sum(B[i][j]*g[j]-D[i][j]*tr[j] for j in range(3)))/loop**2 for i in range(3)]
    ydot=[]
    for i in range(3):
        mixing=arb(0) if i==2 else arb(0).union(tr[1-i])
        G=sum(a*b for a,b in zip(C[i],g))
        ydot.extend((3*(v-mixing)+2*(T-G))/loop for v in y[3*i:3*i+3])
    if len(vec)==12:return gdot+ydot
    l=vec[12]
    hyper=3*g[0]/5
    y4=3*sum(v*v for v in y[:6])+sum(v*v for v in y[6:])
    dl=(24*l*l+(4*T-9*g[1]-3*hyper)*l-2*y4+3*(2*g[1]**2+(g[1]+hyper)**2)/8)/loop
    dm=(12*l+2*T-9*g[1]/2-3*hyper/2)/loop
    return gdot+ydot+[dl,dm]

def evolve(state,delta,keep_last=False):
    n=ceil(abs(float(delta))/.002)
    h=delta/n
    pads=[q('.0003')]*3+[q('.0006')]*9+[q('.00009'),q('.00015')]
    hull=None
    for _ in range(n):
        tube=[x+arb(0,r) for x,r in zip(state,pads)]
        deriv=slopes(tube)
        assert all(abs(h*v)<r for v,r in zip(deriv,pads))
        if keep_last:
            v=[x+arb(0).union(h)*d for x,d in zip(state,deriv)]
            hull=v if hull is None else [a.union(b) for a,b in zip(hull,v)]
        state=[x+h*v for x,v in zip(state,deriv)]
    return hull if keep_last else state


interval=(q(1000)/q('91.1876')).log()
mz=evolve(low_input,-interval)
lam=q('.5573').union(q('.5611'))/4
vev=q('248.372').union(q('248.436'))
one=evolve(mz+[lam,(lam*vev**2).log()],interval)
end=evolve(low_input+one[12:],(q('1171921831891.2').union(q('1171921831891.3')).union(q('1171921831891.1'))/1000).log())
for name,x in [('Higgs quartic',end[12]),('absolute mass squared',end[13].exp())]:
    print(name, x.lower().str(12),x.upper().str(12),flush=True)

