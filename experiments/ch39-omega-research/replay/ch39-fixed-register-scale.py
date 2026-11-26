"""Exact-rational calibration of the common scale on the fixed 46-height register.

The finite comparisons use the archived dyadic height and target tables, the
four shell phase rows, and the archived integer preconditioner.
"""
from fractions import Fraction as Q
from itertools import product
from math import prod
from pathlib import Path
from runpy import run_path
import json


class Table:
    def __init__(self, root):
        table = json.loads((root/"ch39-five-tick-record.json").read_text())
        pre = json.loads((root/"ch39-native-positive-inverse.json").read_text())
        word = run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
        law = run_path(str(root/"ch39-two-prime-height-law.py"))
        matrix = run_path(str(root/"ch39-native-positive-integers.py"))["enclose_matrix"]
        self.n = len(table["rows"])
        self.unit = word["DEN"]
        self.den = (1 << pre["depth"])*self.unit
        self.v = pre["numerators"]
        self.vnorm = Q(max(sum(map(abs,row)) for row in self.v), 1 << pre["depth"])
        height_den = 1 << table["bits"]
        self.z = [Q(height_den**2,r[1]**2) for r in table["rows"]]
        self.input_step = Q(1,height_den)
        self.y = [Q(a,height_den) for a in table["levels"]]
        self.b = [sum(Q(k,1 << pre["depth"])*x for k,x in zip(row,self.y)) for row in self.v]
        self.center = matrix(table,word,law["rate_numerators"]())
        self.R = self.transform(self.center,True)
        self.Rnorm = self.norm(self.R)
        assert self.Rnorm < Q(1,2**240)
        self.taylor = []
        for power in range(1,5):
            rows = []
            for j in range(42):
                vs = self.observation(j)
                rows.append([self.cell(v*(v-1)**power) for v in vs])
            rows += [[(0,0)]*self.n for _ in range(4)]
            self.taylor.append(self.transform(rows,False))
        self.B = self.taylor[0]

    def observation(self,j,s=1):
        return [1/(1+j*j*s*z) for z in self.z]

    def cell(self,x):
        scaled=x*self.unit
        return scaled//1,-((-scaled)//1)

    def transform(self,matrix,remove_identity):
        cols=list(zip(*matrix))
        result=[]
        for i,row in enumerate(self.v):
            entries=[]
            for j,col in enumerate(cols):
                low=high=-self.den if i==j and remove_identity else 0
                for a,(l,h) in zip(row,col):
                    if a<0:
                        l,h=h,l
                    low+=a*l
                    high+=a*h
                entries.append((low,high))
            result.append(entries)
        return result

    def norm(self,matrix):
        return Q(max(sum(max(-l,h) for l,h in row) for row in matrix),self.den)

    def approximate_product(self,matrix,vector):
        # Product with a dyadic vector, preserving cancellations before enclosure.
        unit=self.unit
        integers=[(v*unit)//1 for v in vector]
        gap=int(any(v*unit!=a for v,a in zip(vector,integers)))
        result=[]
        error=Q(0)
        for row in matrix:
            low=high=0
            rounding=0
            for (l,h),v in zip(row,integers):
                rounding+=gap*max(-l,h)
                if v<0:
                    l,h=h,l
                low+=v*l
                high+=v*h
            mid=(low+high)//(2*self.den)
            result.append(Q(mid,unit))
            error=max(error,Q(max(mid*self.den-low,high-mid*self.den)+rounding,self.den*unit))
        return result,error


def restrict_scale(t):
    windows=[]
    centers=[Q(z*2**80//1,2**80) for z in t.z[:7]]
    for i,count,m,rel in [(0,7,8,Q(".000000087")),(1,6,9,Q(".000099")),(2,6,9,Q(".0026"))]:
        degree=2*(count-1+m)+2
        js=[1+round(Q(40*k,degree)) for k in range(degree+1)]
        assert len(set(js))==degree+1
        coeff=[]
        for j in js:
            v=Q(-1,j*j)
            p=v**m*prod(v-c for k,c in enumerate(centers[:count]) if k!=i)
            coeff.append(p*p*(rel**2-(v/centers[i]-1)**2)/prod(1+h*h*v for h in js if h!=j))
        assert sum(a*t.y[j] for a,j in zip(coeff,js))>t.input_step*sum(map(abs,coeff))
        # Numerator is nonpositive outside the windows, and denominators are positive.
        windows.append([((1-rel)*centers[i]/z,(1+rel)*centers[i]/z*(1+t.input_step)**2) for z in t.z])
    triples=[idx for idx in product(range(t.n),repeat=3)
             if max(w[i][0] for w,i in zip(windows,idx)) < min(w[i][1] for w,i in zip(windows,idx))]
    assert triples==[(0,1,2)]
    radius=Q(1,10**7)
    lo,hi=windows[0][0]
    assert lo>1-radius and hi<1+radius
    # Extract a lower bound for the first mass with an increasing rational function.
    js=range(1,6)
    amount=Q(1,prod(j*j for j in js))
    for j in js:
        x=Q(-1,j*j)
        coeff=x**5/prod(1+h*h*x for h in js if h!=j)
        amount+=t.y[j]*coeff-t.input_step*abs(coeff)
    upper=[]
    for z in t.z[:2]:
        x=(1+radius)*z
        upper.append(x**5/prod(1+j*j*x for j in js))
    assert (amount-upper[1])/(upper[0]-upper[1])>Q(1,5)
    # Bound the scaled divided difference row, for every positive solution.
    for i,sign in [(42,1),(43,-1)]:
        row_error=t.vnorm*(radius**4/(1-radius)+5*t.input_step)
        bounds=[]
        for k in range(t.n):
            low,high=t.B[i][k]
            lower=Q(low if sign==1 else -high,t.den)
            tail=sum(radius**d*Q(max(-l,h),t.den)
                     for d,M in enumerate(t.taylor[1:],1) for l,h in [M[i][k]])
            bounds.append(lower-tail-row_error)
        slope=(bounds[0]-min(bounds[1:]))/5+min(bounds[1:])
        assert bounds[0]>min(bounds[1:]) and slope>0
        # b_i = w_i + (Rw)_i + delta*((VM_s - VM_1)/delta * w)_i.
        assert (t.b[i]+t.Rnorm+2*t.vnorm*t.input_step)/slope < Q(1,10**48)
    return Q(1,10**48)


def local_family(t,radius):
    # Derivatives in x=(s-1)/radius; the matrix at the center is I+R.
    beta=radius*t.norm(t.B)
    remainder=t.Rnorm+t.n*t.vnorm*(radius**2/(1-radius)+2*t.input_step)
    derivative_remainder=2*t.n*t.vnorm*(radius**2/(1-radius)**2+radius*t.input_step)
    contraction=beta+remainder
    assert contraction<Q(1,4)
    first,eps=t.approximate_product(t.B,t.b)
    v=[-radius*a for a in first]
    first_error=radius*eps
    bv,error=t.approximate_product(t.B,v)
    bv_bound=radius*(max(map(abs,bv))+error)
    bnorm=max(map(abs,t.b)); vnorm=max(map(abs,v))
    # Residual of the affine trial b+x*v, then the differentiated system.
    err=(first_error+bv_bound+remainder*(bnorm+vnorm)+t.vnorm*t.input_step)/(1-contraction)
    deriv=(first_error+2*bv_bound+beta*err+derivative_remainder*(bnorm+vnorm+err)+remainder*vnorm)/(1-contraction)
    assert max(err,deriv)<Q(1,10**17)
    assert min(t.b)>Q(".00079")
    return v,err,deriv


def actual_weights(t,delta):
    rows=[[t.cell(x) for x in t.observation(j,1+delta)] for j in range(42)] + t.center[42:]
    residual=t.transform(rows,True)
    rho=t.norm(residual)
    assert rho<Q(1,400)
    iteration=t.b[:]
    err=max(map(abs,t.b))*rho/(1-rho)
    for _ in range(32):
        value,uncertainty=t.approximate_product(residual,iteration)
        iteration=[a-v for a,v in zip(t.b,value)]
        err=err*rho+uncertainty
    assert err<Q(1,2**254)
    return iteration,err


def feasible_endpoints(t,v,err,deriv):
    # Both roots are simple; all remaining weights stay positive between these brackets.
    left=(Q("-.456887937193"),Q("-.456887937192"))
    right=(Q(".035498173876"),Q(".035498173877"))
    assert v[43]>deriv and v[42]<-deriv
    assert t.b[43]+left[0]*v[43]<-err and t.b[43]+left[1]*v[43]>err
    assert t.b[42]+right[0]*v[42]>err and t.b[42]+right[1]*v[42]<-err
    assert all(min(t.b[i]+x*v[i] for x in (left[0],right[1]))>err
               for i in range(t.n) if i not in (42,43))
    return left[0],right[1]


def check_holdout(t,radius,v,err,deriv,xlo,xhi,root):
    u=t.observation(42)
    # Compress the old 46-by-46 inverse onto the held-out row, prior to interval evaluation.
    aa=[(f*t.unit)//1 for f in u]
    pre_depth=t.den//t.unit
    dual=[Q(sum(n*row[j] for n,row in zip(aa,t.v)),t.unit*pre_depth) for j in range(t.n)]
    size=sum(map(abs,dual))
    base=sum(a*y for a,y in zip(dual,t.y))
    record=json.loads((root/"ch39-minimum-tick-finite-inputs.json").read_text())
    data_den=2**record["depth"]
    assert all(Q(level,data_den)==y for level,y in zip(record["levels"],t.y))
    expected=Q(record["levels"][42],data_den)
    # Error both at the target and in the compressed relation u - dual*M_1.
    noise=t.n*t.Rnorm+Q(1,t.unit)+(2*size+2)*t.input_step
    cross=[]
    for i,z in enumerate(t.z):
        ds=[]
        for j in range(42):
            f=1/(1+j*j*z)
            ds.append(f*(f-1))
        cross.append(u[i]*(u[i]-1)-sum(a*x for a,x in zip(dual,ds)))
    gradient_noise=(size+1)*(2*t.input_step+radius/(1-radius))
    firstlo=t.b[0]-abs(v[0])-err
    firsthi=t.b[0]+abs(v[0])+err
    crosslo=min(cross[1:]); crosshi=max(cross[1:])
    assert cross[0]>crosshi
    slow=firstlo*(cross[0]-crosslo)+crosslo-gradient_noise
    shigh=firsthi*(cross[0]-crosshi)+crosshi+gradient_noise
    assert Q(".4749341")<slow<shigh<Q(".4749344")
    deficit_lo=expected-base-noise-xhi*radius*shigh
    deficit_hi=expected-base+noise-xlo*radius*shigh
    assert Q("3.72e-49")<deficit_lo<deficit_hi<Q("6.07e-49")

    r=t.observation(1682)
    value=sum(a*b for a,b in zip(r,t.b))
    slope=sum(a*b for a,b in zip(r,v))
    cap=t.n*err+t.n*(max(map(abs,t.b))+max(map(abs,v))+err)*(radius/(1-radius)+t.input_step)
    slope_noise=t.n*deriv+t.n*(max(map(abs,v))+deriv)*(radius/(1-radius)+t.input_step)
    slope_noise+=t.n*(max(map(abs,t.b))+max(map(abs,v))+err)*radius/(1-radius)**2
    assert Q(".00000988")<slope-slope_noise<slope+slope_noise<Q(".00000989")
    lo,hi=value+xlo*slope-cap,value+xhi*slope+cap
    assert Q(".0321412")<lo<hi<Q(".0321462")
    observed=json.loads((root/"ch39-fixed-register-extra-target.json").read_text())
    actual=Q(observed["floor"],2**observed["bits"])
    assert observed["j"]==1682 and Q(".03597252")<actual
    assert actual-hi>Q(".0038263")
    print("Every positive scale and every nonnegative reweighting on these 46 height cells fails the next arithmetic reading.")
    print("The predicted deficit at j=42 is in (3.72,6.07) times 10^-49. At 1682 it exceeds .0038263.")
    return slope


def main():
    root=Path(__file__).resolve().parent
    t=Table(root)
    radius=restrict_scale(t)
    v,err,deriv=local_family(t,radius)
    lo,hi=feasible_endpoints(t,v,err,deriv)
    check_holdout(t,radius,v,err,deriv,lo,hi,root)
    print("Feasible interval ends at weights for heights 291 and 278, in that order.")


if __name__=="__main__":
    main()

