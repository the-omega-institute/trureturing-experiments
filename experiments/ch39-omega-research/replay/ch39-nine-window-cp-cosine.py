"""以一个相对奇异谱迹在六质量三角度全窗中限制夸克 CP 相位。"""
from pathlib import Path
from itertools import product
from runpy import run_path
from flint import arb

class Dual:
    """Arb 区间上的一元导数，其他输入同时在给定盒内变化。"""
    def __init__(self, v, d=0):
        self.v,self.d = arb(v),arb(d)
    @staticmethod
    def cast(a):
        return a if isinstance(a,Dual) else Dual(a)
    def __add__(self,b):
        b=self.cast(b);return Dual(self.v+b.v,self.d+b.d)
    __radd__=__add__
    def __neg__(self):return Dual(-self.v,-self.d)
    def __sub__(self,b):return self+-self.cast(b)
    def __rsub__(self,b):return self.cast(b)+-self
    def __mul__(self,b):
        b=self.cast(b);return Dual(self.v*b.v,self.d*b.v+self.v*b.d)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=self.cast(b);return Dual(self.v/b.v,(self.d*b.v-self.v*b.d)/b.v**2)
    def __rtruediv__(self,b):return self.cast(b)/self
    def __pow__(self,p):return Dual(self.v**p,p*self.v**(p-1)*self.d)
    def sin(self):return Dual(self.v.sin(),self.d*self.v.cos())
    def cos(self):return Dual(self.v.cos(),-self.d*self.v.sin())


def affine_trace(values, inverse=False):
    orientation = -1 if inverse else 1
    masses=[[(v*v/(group[(i+1)%3]*group[(i+2)%3]))**(arb(2)*sign/3)
             for i,v in enumerate(group)] for group,sign in ((values[:3],orientation),(values[3:6],-orientation))]
    u,d=masses
    a,b,z=[v.sin() for v in values[6:9]]
    ca,cb,cz=[v.cos() for v in values[6:9]]
    x,y,w=a*a,b*b,z*z
    m0=((ca*ca*cz*cz,x*cz*cz,w),
        (x*cb*cb+ca*ca*y*w,ca*ca*cb*cb+x*y*w,y*cz*cz),
        (x*y+ca*ca*cb*cb*w,ca*ca*y+x*cb*cb*w,cb*cb*cz*cz))
    base=sum(u[i]*d[j]*m0[i][j] for i,j in product(range(3),repeat=2))
    slope=2*a*b*z*ca*cb*(u[1]-u[2])*(d[0]-d[1])
    return base,slope


def phase(values,target):
    base,slope = affine_trace(values)
    return (-base+target)/slope

def verify(replay):
    m,u,d=run_path(str(replay/'ch39-two-zero-orbit-condition.py'))['calibration']()
    J=u*d.inv()
    trace=sum(abs(x)**2 for x in J.entries())/abs(J.det())**(arb(2)/3)
    ref,w=m['ref'][:9],m['tol'][:9]
    box=[v+arb(0,r.upper()) for v,r in zip(ref,w)]
    signs = (-1,-1,1,-1,1,1,1,1)
    for i,sign in enumerate(signs):
        result = phase([Dual(v,int(i==j)) for j,v in enumerate(box)],trace)
        assert sign*result.d>0
    extrema = []
    for pm in (-1,1):
        values = [v+pm*sign*r for v,r,sign in zip(ref[:8],w[:8],signs,strict=True)]
        for k in range(32):
            left = ref[8]+w[8]*(arb(2)*k/32-1)
            right = ref[8]+w[8]*(arb(2)*(k+1)/32-1)
            result = phase([Dual(x) for x in values]+[Dual(left.union(right),1)],trace)
            assert -pm*result.d>0
        extrema.append(phase(values+[ref[8]-w[8]],trace))
    low, high = extrema
    assert arb(".30755") < low < high < arb(".51619")
    assert arb("1.02839") < high.acos() < low.acos() < arb("1.25818")
    # 在角和其补角严格为正的标准约定下，J 的绝对值。
    a,b,z = box[6:9]
    absolute_j = a.sin()*a.cos()*b.sin()*b.cos()*z.sin()*z.cos()**2*(1-high**2).sqrt()
    assert absolute_j > arb(".000029174")
    print("单正迹约束的余弦端点",low.str(12),high.str(12),"J 绝对值下界",absolute_j.lower().str(12))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
