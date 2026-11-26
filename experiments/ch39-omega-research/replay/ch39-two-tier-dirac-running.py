"""两档磁壳的反常检验和矢量状重层的尖阈值两圈规范交点。"""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


ROOT = Path(__file__).parent


def matter():
    # (三代总份数，色维数，弱维数，左手超荷，色基本或反基本的三次符号)
    reps = ((6,3,2,F(1,6),1), (3,3,1,-F(2,3),-1), (3,1,1,F(1),0))
    anomaly = [F(0)]*5
    one = [F(0)]*3
    two = [[F(0)]*3 for _ in range(3)]
    adjoint = (0,2,3)
    for count,color,weak,y,sign in reps:
        dim = color*weak
        indices = (F(3,5)*dim*y*y, F(color,2) if weak==2 else 0,
                   F(weak,2) if color==3 else 0)
        casimirs = (F(3,5)*y*y,F(3,4) if weak==2 else 0,F(4,3) if color==3 else 0)
        for i,amount in enumerate((weak*sign,indices[2]*y,indices[1]*y,dim*y,dim*y**3)):
            anomaly[i] += count*amount
        for i in range(3):
            one[i] += count*F(4,3)*indices[i]
            for j in range(3):
                two[i][j] += count*indices[i]*(4*casimirs[j]+(F(20,3)*adjoint[i] if i==j else 0))
    assert anomaly == [9,0,F(3,2),3,F(1,2)]
    assert one == [F(32,5),12,10]
    assert two == [[F(173,25),F(9,5),16],[F(3,5),147,48],[2,18,190]]
    # 每个标准模型多重态的 X=2Y-B，包括三个家族的全部轻 Weyl 场。
    family = [(6,F(1,6),F(1,3)),(3,-F(2,3),-F(1,3)),(3,F(1,3),-F(1,3)),
              (2,-F(1,2),0),(1,F(1),0)]
    x = [2*y-baryon for _,y,baryon in family]
    assert x == [0,-1,1,-1,2]
    assert sum(dim*xi for (dim,_,_),xi in zip(family,x)) == 0
    assert 3*sum(dim*xi**3 for (dim,_,_),xi in zip(family,x)) == 18
    assert (3*x[0]+x[3])%3 == 2
    return one,two


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

def integrate_segments(initial, segments, D, C, record=None, *, min_inverse=8):
    pi = arb.pi()
    correction = [arb(0) for _ in range(3)]
    logs = [arb(x).log() for x in (".98",".016",".010")]
    baseline = initial[:]
    start_time = arb(0)
    for time,coeff,matrix,steps in segments:
        dt = time/steps
        for k in range(steps):
            interval = (arb(k)+arb(0).union(1))*dt
            # 常数宽度包络只用于证明单步内解存在并给导数定界。
            xpad, ypad = arb(".4")*dt.upper(), arb(".12")*dt.upper()
            inverse = [a-s*interval/(2*pi)+v+arb(0,xpad) for a,s,v in zip(baseline,coeff,correction)]
            assert min_inverse>0 and all(v>min_inverse for v in inverse)
            gauge = [1/v for v in inverse]
            y = [(2*(v+arb(0,ypad))).exp() for v in logs]
            if record is not None:
                record(start_time+k*dt,dt,gauge,y)
            trace = 3*(y[0]+y[1])+y[2]
            dx = [-dot(row,gauge)/(8*pi**2)+dot(yr,y)/(32*pi**3) for row,yr in zip(matrix,D)]
            dy = [(arb(3)/2*diff+trace-4*pi*dot(row,gauge))/(16*pi**2)
                  for diff,row in zip((y[0]-y[1],y[1]-y[0],y[2]),C)]
            for values,speeds,pad in [(correction,dx,xpad),(logs,dy,ypad)]:
                for i,speed in enumerate(speeds):
                    assert abs(dt*speed) < pad
                    values[i] += dt*speed
        baseline = [a-s*time/(2*pi) for a,s in zip(baseline,coeff)]
        start_time += time
    return correction,[v+h for v,h in zip(baseline,correction)]



def verify(observer=None):
    sm = run_path(str(ROOT/"ch39-two-loop-shell-separation.py"))
    ctx.prec = 192
    q = sm["q"]
    extra_one, extra_two = matter()
    b, B, D, C = sm["BETA"], sm["B"], sm["D"], sm["C"]
    d = list(map(q,extra_one))
    b_high = [a+k for a,k in zip(b,d)]
    high = [[a+q(k) for a,k in zip(row,extra)] for row,extra in zip(B,extra_two)]
    pi = arb.pi()
    em,weak = arb("127.930"),arb(".23129")
    initial = [(1-weak)*em*3/5,weak*em,1/arb(".1180")]
    slopes = arb_mat([[b[i]-b[2],b_high[i]-b_high[2]] for i in (0,1)])
    sharp_one = slopes.solve(arb_mat([[2*pi*(initial[i]-initial[2])] for i in (0,1)]))
    one_low, one_high = sharp_one[0,0],sharp_one[1,0]
    one_m = arb("91.1876")*one_low.exp()
    one_mu = arb("91.1876")*(one_low+one_high).exp()
    assert arb("1.99047e12") < one_m < arb("1.99049e12")
    assert arb("2.40849e15") < one_mu < arb("2.40851e15")
    print("相同尖阈值处方的一圈轻档与交点",one_m.str(13),one_mu.str(13))
    def integrate(times,steps,record=None):
        segments = [(time,coeff,matrix,steps) for time,coeff,matrix in zip(times,(b,b_high),(B,high))]
        return integrate_segments(initial,segments,D,C,record)

    # 在时间矩形上反演精确的一圈仿射项，两圈积分进入连续映射右端。
    times = [arb("24.25").union(arb("24.62")),arb("5.80").union(arb("6.14"))]
    corr,_ = integrate(times,512)
    rhs = arb_mat([[2*pi*(initial[i]-initial[2]+corr[i]-corr[2])] for i in (0,1)])
    mapped = slopes.solve(rhs)
    assert all(t.contains_interior(mapped[i,0]) for i,t in enumerate(times))
    for _ in range(3):
        times = [times[i].intersection(mapped[i,0]) for i in range(2)]
        corr,ended = integrate(times,1024)
        rhs = arb_mat([[2*pi*(initial[i]-initial[2]+corr[i]-corr[2])] for i in (0,1)])
        mapped = slopes.solve(rhs)
    times = [times[i].intersection(mapped[i,0]) for i in range(2)]
    _,ended = integrate(times,1024 if observer is None else 4096,observer)
    common = ended[0].intersection(ended[1]).intersection(ended[2])
    m = arb("91.1876")*times[0].exp()
    mu = arb("91.1876")*sum(times).exp()
    assert arb("3.63e12") < m < arb("3.81e12")
    assert arb("1.38e15") < mu < arb("1.54e15")
    assert arb("32.78") < common < arb("32.88")
    print("两段正演化时间",*[v.str(9) for v in times])
    print("轻档 GeV、交点 GeV、共同逆耦合的区间",
          *[(v.lower().str(10),v.upper().str(10)) for v in (m,mu,common)])
    return m,mu,common


if __name__ == "__main__":
    verify()
