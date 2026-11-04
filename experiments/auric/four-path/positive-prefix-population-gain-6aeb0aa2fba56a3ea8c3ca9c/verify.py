"""Independent Fraction-interval / automatic-differentiation verifier.

Reads only this scope's certificate. Uses exact rational arithmetic between
square roots, unlike certify.py's rounded integer operations. For slope
certificates it differentiates the original omega expression automatically.
"""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
from hashlib import sha256
import json

ROOT_SCALE=1<<128


class V:
    def __init__(self, lo, hi=None):
        self.lo=Q(lo)
        self.hi=self.lo if hi is None else Q(hi)
        assert self.lo<=self.hi

    @staticmethod
    def cast(v):
        return v if isinstance(v,V) else V(v)

    def __add__(self,v):
        v=self.cast(v)
        return V(self.lo+v.lo,self.hi+v.hi)
    __radd__=__add__
    def __neg__(self):
        return V(-self.hi,-self.lo)
    def __sub__(self,v):
        return self+-self.cast(v)
    def __rsub__(self,v):
        return self.cast(v)+-self
    def __mul__(self,v):
        v=self.cast(v)
        p=[a*b for a in (self.lo,self.hi) for b in (v.lo,v.hi)]
        return V(min(p),max(p))
    __rmul__=__mul__
    def recip(self):
        assert self.lo>0 or self.hi<0
        return V(1/self.hi,1/self.lo)
    def __truediv__(self,v):
        return self*self.cast(v).recip()
    def __rtruediv__(self,v):
        return self.cast(v)*self.recip()
    def sq(self):
        return V(0 if self.lo<=0<=self.hi else min(self.lo**2,self.hi**2),max(self.lo**2,self.hi**2))
    def sqrt(self):
        assert self.lo>=0
        a=isqrt((self.lo*ROOT_SCALE**2).numerator//(self.lo*ROOT_SCALE**2).denominator)
        b=isqrt((self.hi*ROOT_SCALE**2).numerator//(self.hi*ROOT_SCALE**2).denominator)
        if Q(b*b,ROOT_SCALE**2)!=self.hi:b+=1
        return V(Q(a,ROOT_SCALE),Q(b,ROOT_SCALE))


class Jet:
    def __init__(self,v,d=0):
        self.v=V.cast(v)
        self.d=V.cast(d)
    @staticmethod
    def cast(a):
        return a if isinstance(a,Jet) else Jet(a)
    def __add__(self,a):
        a=self.cast(a)
        return Jet(self.v+a.v,self.d+a.d)
    __radd__=__add__
    def __neg__(self):return Jet(-self.v,-self.d)
    def __sub__(self,a):return self+-self.cast(a)
    def __rsub__(self,a):return self.cast(a)+-self
    def __mul__(self,a):
        a=self.cast(a)
        return Jet(self.v*a.v,self.d*a.v+self.v*a.d)
    __rmul__=__mul__
    def recip(self):return Jet(self.v.recip(),-self.d/self.v.sq())
    def __truediv__(self,a):return self*self.cast(a).recip()
    def __rtruediv__(self,a):return self.cast(a)*self.recip()
    def sq(self):return Jet(self.v.sq(),2*self.v*self.d)
    def sqrt(self):
        r=self.v.sqrt()
        return Jet(r,self.d/(2*r))


def expression(name, box):
    derivative=name=='slope'
    cons=Jet if derivative else V.cast
    xi=cons(V(*box[0])); q=cons(V(*box[1]),1) if derivative else V(*box[1])
    h=cons(V('1/2').sqrt());c=((1+h)/2).sqrt();s2=(1-h)/2;s=s2.sqrt();k=c-h
    x=1-xi/4+q/2; y=(1-x.sq()).sqrt();u=(1-q).sqrt();v=q.sqrt();R=h*x+c
    if name=='gain':
        return k.sq()*u/(8*R)+2*s*v*(1+x)/(3*cons(V(3).sqrt())*x)
    if name=='scalar_derivative':
        Gamma=u*R+h*v*y
        return -(4*h+8*s)-4*Gamma*(h*u-R/u+h*y/v-h*x*v/y)
    mu=(h*k*u/R+h*v/y)/2
    nu=(k*u+(x/2-s2)*v/(h*y))/2
    omega=nu-s2*q*(1+x)/(4*mu*x.sq())
    return omega.d if derivative else omega


def volume(box,active):
    result=Q(1)
    for j in active:result*=box[j][1]-box[j][0]
    return result


def coverage(row):
    domain=[tuple(map(Q,b)) for b in row['domain']]
    boxes=[[tuple(map(Q,b)) for b in leaf['box']] for leaf in row['leaves']]
    active=[j for j,(a,b) in enumerate(domain) if a<b]
    for box in boxes:assert all(a<=c<=d<=b for (a,b),(c,d) in zip(domain,box))
    assert sum((volume(b,active) for b in boxes),Q(0))==volume(domain,active)
    if not active:
        assert len(boxes)==1
    else:
        j=active[0]
        boxes.sort(key=lambda b:b[j][0])
        for index,box in enumerate(boxes):
            for other in boxes[index+1:]:
                if other[j][0]>=box[j][1]:break
                assert any(min(box[i][1],other[i][1])<=max(box[i][0],other[i][0]) for i in active)


def main():
    path=Path(__file__).with_name('certificate.json')
    data=json.loads(path.read_text())
    results=[]
    for row in data['certificates']:
        coverage(row)
        initial=[tuple(map(Q,b)) for b in row['domain']]
        low=None if row['lower'] is None else Q(row['lower'])
        high=None if row['upper'] is None else Q(row['upper'])
        todo=[([tuple(map(Q,b)) for b in leaf['box']],0) for leaf in row['leaves']]
        verified=[]
        while todo:
            box,depth=todo.pop()
            image=expression(row['expression'],box)
            if (low is None or image.lo>=low) and (high is None or image.hi<=high):
                verified.append(image)
                continue
            assert depth<20,(row['name'],box,depth)
            widths=[(b-a)/(d-c) if c<d else Q(-1) for (a,b),(c,d) in zip(box,initial)]
            axis=max(range(2),key=lambda j:widths[j])
            a,b=box[axis];assert a<b
            mid=(a+b)/2
            left=box.copy();right=box.copy();left[axis]=(a,mid);right[axis]=(mid,b)
            todo.extend([(left,depth+1),(right,depth+1)])
        # Small rational result enclosures, rounded out to 10**-12.
        scale=10**12
        lo=min(v.lo for v in verified);hi=max(v.hi for v in verified)
        L=(lo*scale).numerator//(lo*scale).denominator
        H=-(((-hi*scale).numerator)//((-hi*scale).denominator))
        summary={'name':row['name'],'input_leaves':len(row['leaves']),
                 'verified_subleaves':len(verified),'lower':str(Q(L,scale)),
                 'upper':str(Q(H,scale))}
        results.append(summary)
        print(summary,flush=True)
    report={'verification':'PASS','arithmetic':'exact Fraction intervals; 128-bit outward square roots; original-expression automatic differentiation',
            'certificate_sha256':sha256(path.read_bytes()).hexdigest(),
            'verifier_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'results':results}
    Path(__file__).with_name('verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: independent arithmetic, automatic derivatives and complete closed-domain coverage')


if __name__=='__main__':main()
