"""检验与两个原生 Gram 算子的任意非零实线性组合对易的源协方差。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def verify():
    m = run_path(str(Path(__file__).with_name("ch39-two-zero-quark-window.py")))
    adj = m["adj"]
    matrices = [m["prepare"](q,arb(x),arb(y)) for q,x,y in (
        (-4,"-0.680061861483031","0.873223965353125"),
        (2,"0.565247837594779","-0.630343671648397"),
    )]
    factors = [a/abs(a.det())**(arb(1)/3) for a in matrices]
    up, down = [adj(f)*f for f in factors]
    def frob(a):
        return sum(abs(z)**2 for z in a.entries()).sqrt().upper()
    up_norm, down_norm = [frob(h).sqrt() for h in (up,down)]
    comm = abs((up*down-down*up).det())
    boxes = [v+arb(0,w.upper()) for v,w in zip(m["ref"],m["tol"],strict=True)]
    def normalized_squared(y):
        factor = (y[0]*y[1]*y[2])**(arb(2)/3)
        return [v*v/factor for v in y]
    a,b = [normalized_squared(row) for row in (boxes[:3],boxes[3:6])]
    def discriminant(row):
        x,y,z = row
        return (z-x)*(z-y)*(y-x)
    v = acb_mat(m["load"]("shell-clock-ckm-blocks")["standard_matrix"](
        *(x.sin() for x in boxes[6:9]),boxes[9]))
    cp = abs((v[0,1]*v[1,1].conjugate()*v[1,2]*v[0,2].conjugate()).imag)
    observed_comm = 2*cp*discriminant(a)*discriminant(b)
    split = arb(".1")
    # 两个最重谱项的和秩至多二，非零两根之积的绝对值等于 p*q*loss。
    loss = boxes[7].sin()**2+boxes[7].cos()**2*boxes[8].sin()**2
    p, q = split*a[2].upper(), b[2].upper()
    tail = split*a[1].upper()+b[1].upper()
    spread = p+q+2*tail
    gaps_product = spread*(p*q*loss.upper()+2*tail*(p+q)+4*tail**2)
    reduced_comm = comm*gaps_product / observed_comm
    target = sum(a)/sum(b)
    assert target > arb("74.25")
    print("低权重端目标判据的统一上界", reduced_comm.upper().str(12))

    assert frob(down) < 1024
    shift = acb_mat([[1024*int(i==j) for j in range(3)] for i in range(3)])
    for sign in (1,-1):
        rate = frob(up-sign*down)
        for start,stop,low in ((arb(0),split,True),(split,arb(1),False)):
            stack = [(start,stop,0)]
            count = 0
            while stack:
                left,right,depth = stack.pop()
                mid = (left+right)/2
                half = (right-left)/2
                levels,Q = m["eig"](mid*up+sign*(1-mid)*down+shift)
                change = rate*half
                if low:
                    gaps = [levels[j]-levels[i]-2*change for i,j in ((0,1),(0,2),(1,2))]
                    passes = min(gaps)>0 and gaps[0]*gaps[1]*gaps[2] > reduced_comm
                else:
                    passes = True
                    for i in range(3):
                        vec = acb_mat([[Q[j,i]] for j in range(3)])
                        gap = min(abs(levels[j]-levels[i]).lower() for j in range(3) if i!=j)
                        if not gap > 2*change:
                            passes = False
                            break
                        residual = half*frob((up-sign*down)*vec)
                        distance = arb(2).sqrt()*residual/(gap-change)
                        upper = frob(factors[0]*vec)+distance*up_norm
                        lower = sum(abs(z)**2 for z in (factors[1]*vec).entries()).sqrt().lower()-distance*down_norm
                        if not lower>0 or not upper**2 < arb("74.25")*lower**2:
                            passes = False
                            break
                if passes:
                    count += 1
                else:
                    assert depth < 26,(left,right)
                    stack.extend(((left,mid,depth+1),(mid,right,depth+1)))
            print("制备权重分段、终格数",sign,start.str(5),stop.str(5),count,flush=True)


if __name__ == "__main__":
    verify()
