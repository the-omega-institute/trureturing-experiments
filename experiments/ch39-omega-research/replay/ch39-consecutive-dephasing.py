"""Sharp consecutive equal-clock epoch count excluding the second-moment scale."""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
model=run_path(str(root/"ch39-clock-schur-mediators.py"))
B,root_cell,base=[model[k] for k in ("B","root_cell","base")]
round_out=model["clock"]["rounded"]
word=base["word"]
unit=word["DEN"]


def positive(matrix):
    # Fraction-free elimination. The successive pivots are leading principal minors.
    previous=1
    while len(matrix)>1:
        pivot=matrix[0][0]
        assert pivot>0
        out=[]
        for i in range(1,len(matrix)):
            row=[]
            for j in range(1,len(matrix)):
                v,r=divmod(pivot*matrix[i][j]-matrix[i][0]*matrix[0][j],previous)
                assert r==0
                row.append(v)
            out.append(row)
        matrix,previous=out,pivot
    assert matrix[0][0]>0


def coordinates(native):
    return tuple(B(Q(l,unit),Q(h,unit)) for l,h in native)


def dual_check(n,row,p,phases):
    y=[Q(v,2**32) for v in row]
    assert all(v>0 for v in y) and sum(y)<Q("2.410")
    re,im=[],[]
    for t in range(n):
        a,b=B(int(t==0)),B(0)
        for v,price,x in zip(p,y,phases[t]):
            u,w=coordinates(x)
            a=round_out(a-v*u/B(n*price))
            b=round_out(b-v*w/B(n*price))
        re.append(a)
        im.append(round_out(b*root_cell(B(17))))
    grid=2**56
    def approx(v):
        floor=v.lo*grid//1
        assert v.hi*grid<floor+2
        return floor
    real=[approx(v) for v in re]
    imaginary=[approx(v) for v in im]
    assert im[0].lo==im[0].hi==0
    matrix=[]
    for i in range(2*n):
        row=[]
        for j in range(2*n):
            t=i%n-j%n
            if (i<n)==(j<n):v=real[abs(t)]
            else:v=(1 if i>=n else -1)*(1 if t>=0 else -1)*imaginary[abs(t)]
            row.append(v-(4*n if i==j else 0))
        matrix.append(row)
    positive(matrix)


def main():
    ids,z,p,*_=model["calibration"]()
    oracle=base["data"]("ch39-consecutive-dephasing.json")
    ns=run_path(str(root/"ch39-native-clock-dephasing.py"))["exponents"](ids)
    u=list(map(word["power"],ns))
    phases=[[word["ONE"]]*len(ids)]
    for _ in range(29):
        phases.append([word["multiply"](a,b) for a,b in zip(phases[-1],u)])
    amp=list(map(root_cell,p))
    target=sum((v*x for v,x in zip(p,z)),B(0))*sum((v/x for v,x in zip(p,z)),B(0))
    assert Q("2.410")<target.lo<=target.hi<Q("2.411")
    for n in range(1,30):
        if n<15:
            register=list(map(word["power"],oracle["lower"][str(n)]))
            score=B(0)
            for row in phases[:n]:
                real,imag=B(0),B(0)
                for x,y,s in zip(row,register,amp,strict=True):
                    a,b=coordinates(word["multiply"](x,y))
                    real=round_out(real+s*a)
                    imag=round_out(imag+s*b)
                score=round_out(score+real*real+B(17)*imag*imag)
            assert score.lo>n*target.hi
        elif str(n) in oracle["upper"]:
            dual_check(n,oracle["upper"][str(n)],p,phases)
        else:
            score=B(1)
            for i,j in combinations(range(len(p)),2):
                def gap(row):
                    q=word["multiply"](row[i],word["conjugate"](row[j]))[0]
                    return B(Q(unit-q[1],unit),Q(unit-q[0],unit))
                size=root_cell(gap(phases[n])/gap(u))
                score=round_out(score+B(Q(2,n))*amp[i]*amp[j]*size)
            assert score.hi<target.lo
    print("One to fourteen consecutive prime-2 clock epochs admit the second-moment scale.")
    print("All contiguous runs of at least fifteen epochs exclude it, in every starting phase.")


if __name__=="__main__":
    main()
