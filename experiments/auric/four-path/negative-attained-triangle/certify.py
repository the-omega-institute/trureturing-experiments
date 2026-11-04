"""Exact Arb enclosure of the original negative attained-triangle comparison.
Python >=3.12; python-flint==0.9.0. No binary floating point enters acceptance.
"""
import sys,json,time,hashlib,argparse
from fractions import Fraction as F
from flint import arb,ctx,fmpq
ctx.prec=128

def A(x):
    if isinstance(x,F): return arb(fmpq(x.numerator,x.denominator))
    return arb(x)
def ball(lo,hi): return arb(fmpq((lo+hi).numerator,(lo+hi).denominator)/2,fmpq((hi-lo).numerator,(hi-lo).denominator)/2)
class D:
    def __init__(self,v,g=None): self.v=v if isinstance(v,arb) else A(v);self.g=g or (A(0),)*3
    def __add__(self,o):
        o=asD(o);return D(self.v+o.v,tuple(a+b for a,b in zip(self.g,o.g)))
    __radd__=__add__
    def __neg__(self): return D(-self.v,tuple(-a for a in self.g))
    def __sub__(self,o): return self+-asD(o)
    def __rsub__(self,o): return asD(o)+-self
    def __mul__(self,o):
        o=asD(o);return D(self.v*o.v,tuple(a*o.v+self.v*b for a,b in zip(self.g,o.g)))
    __rmul__=__mul__
    def __pow__(self,n):
        assert n>=0
        if n==0:return D(1)
        v=self.v**n;return D(v,tuple(n*self.v**(n-1)*a for a in self.g))
def asD(x):return x if isinstance(x,D) else D(x)
h=A(F(1,2)).sqrt();s=((1-h)/2).sqrt();c=((1+h)/2).sqrt();k=h-s
e2=4*(h+s)**2;d2=4*(h+c)**2

def poly(t,a,b,mode):
    u=a;v=(1-a)*b;r=1-t*u;z=1-t*v
    hc=A(F(4,25))+A(F(2,3))*t
    H0=t*hc
    dd=A(F(37,150))+A(F(2,3))*t if mode==0 else A(F(17,60))+t*A(F(1,3))
    aa=e2-H0;bb=d2-H0-dd
    rdef=4*s*h-2*s*k*u+2*h*k*v-t*(s*h-s*k*u*u+h*k*v*v)
    W=4*rdef-hc
    cb=c*c*u*(2-t*u)+h*h*v*(2-t*v)
    P=W*bb-4*aa*cb
    N=(1-a)**2*a*b*(1-b)*(2-t*u)*(2-t*v)*(4-t*(1+u+v))*((2-t)**2-t*t*(u-v)**2)
    Q=r*r*z*z*P*P-16*aa*aa*c*c*h*h*t*N
    return W,bb,P,Q

def evaluate(box,mode):
    xx=[ball(lo,hi) for lo,hi in box]
    natural=poly(*xx,mode)
    if all(x>0 for x in natural):return True,None,[str(x.lower()) for x in natural]
    dx=[D(x,tuple(A(int(i==j)) for j in range(3))) for i,x in enumerate(xx)]
    vals=poly(*dx,mode)
    mid=[A((lo+hi)/2) for lo,hi in box];center=poly(*mid,mode)
    refined=[];scores=[A(0)]*3
    for raw,v,cc in zip(natural,vals,center):
        mean=cc
        for j,(lo,hi) in enumerate(box):mean+=v.g[j]*ball(-(hi-lo)/2,(hi-lo)/2)
        vv=raw.intersection(mean);refined.append(vv)
        if not vv>0:
            for j,(lo,hi) in enumerate(box): scores[j]+=abs(v.g[j])*A((hi-lo)/2)
    # Float conversions select subdivision only, never acceptance.
    split=max(range(3),key=lambda j:float(scores[j].upper()))
    return all(v>0 for v in refined),split,[str(x.lower()) for x in refined]

def roots():
    return [[(F(0),F(11,100)),(F(0),F(1)),(F(0),F(1))],[(F(11,100),F(593,1000)),(F(0),F(1)),(F(0),F(1))]]
def splitbox(box,j):
    lo,hi=box[j];mid=(lo+hi)/2;b0=list(box);b1=list(box);b0[j]=(lo,mid);b1[j]=(mid,hi);return b0,b1

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--seconds',type=int,default=600);ap.add_argument('--output',default='certificate.json');ap.add_argument('--verify');args=ap.parse_args()
    start=time.monotonic()
    if args.verify:
        cert=json.load(open(args.verify));n=0
        for mode,tree in enumerate(cert['trees']):
            stack=[(roots()[mode],tree)]
            while stack:
                box,node=stack.pop()
                if node==1:
                    ok,_,_=evaluate(box,mode);assert ok;n+=1
                else:
                    j,l,r=node;b0,b1=splitbox(box,j);stack.extend([(b0,l),(b1,r)])
        print(json.dumps({'verified_leaves':n,'seconds':time.monotonic()-start}));return
    trees=[None,None];stack=[]
    for mode,box in enumerate(roots()):stack.append((box,mode,trees,mode,0))
    leaves=nodes=0;maxdepth=0
    while stack:
        box,mode,parent,index,depth=stack.pop();nodes+=1
        ok,j,lower=evaluate(box,mode)
        if ok:parent[index]=1;leaves+=1;maxdepth=max(maxdepth,depth)
        else:
            if time.monotonic()-start>args.seconds:
                print(json.dumps({'status':'partial','nodes':nodes,'leaves':leaves,'pending':len(stack)+1,'depth':depth,'box':[[str(a),str(b)] for a,b in box],'lower':lower,'seconds':time.monotonic()-start}));return
            node=[j,None,None];parent[index]=node;b0,b1=splitbox(box,j);stack.extend([(b0,mode,node,1,depth+1),(b1,mode,node,2,depth+1)])
        if nodes%10000==0:print(json.dumps({'nodes':nodes,'leaves':leaves,'pending':len(stack),'seconds':time.monotonic()-start}),flush=True)
    out={'status':'complete','arithmetic':'Arb 128 bit, python-flint 0.9.0','domain':'t in [0,593/1000], a,b in [0,1], split at t=11/100','leaves':leaves,'nodes':nodes,'maxdepth':maxdepth,'trees':trees}
    with open(args.output,'w') as f:json.dump(out,f,separators=(',',':'));f.write('\n')
    print(json.dumps({k:v for k,v in out.items() if k!='trees'}));print('seconds',time.monotonic()-start)
if __name__=='__main__':main()
