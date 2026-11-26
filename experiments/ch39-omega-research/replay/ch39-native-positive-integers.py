"""Integer-only existence and positivity certificate for the exact finite joint table."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path
import json


def enclose_matrix(record, words, rates):
    den = 1 << record["bits"]
    unit = words["DEN"]
    rows = []
    for j in range(42):
        row=[]
        for _,h,_,_,_ in record["rows"]:
            a,b = h*h*unit,h*h+(j*den)**2
            row.append((a//b,-((-a)//b)))
        rows.append(row)
    angular=[[] for _ in range(4)]
    microstep=-43235114615559591407865398948523
    for _,h,_,_,_ in record["rows"]:
        m,n=[microstep*((a*h+(den<<23))//(den<<24)) for a in rates]
        for row,exponent in zip(angular,(m,n,n-m,n+m)):
            row.append(words["power"](exponent)[0])
    return rows+angular


def main():
    root = Path(__file__).resolve().parent
    words = run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    law = run_path(str(root/"ch39-two-prime-height-law.py"))
    record = json.loads((root/"ch39-five-tick-record.json").read_text())
    pre = json.loads((root/"ch39-native-positive-inverse.json").read_text())
    matrix = enclose_matrix(record,words,law["rate_numerators"]())
    den = words["DEN"] * 2**pre["depth"]
    target = record["levels"] + [0]*4
    trialden = (1 << record["bits"])*2**pre["depth"]
    largest = 0
    trial = []
    for i,row in enumerate(pre["numerators"]):
        assert len(row)==len(matrix)==len(pre["numerators"])==46
        row_error = 0
        for j in range(46):
            left = right = -int(i==j)*den
            for a,col in zip(row,matrix):
                lo,hi = col[j]
                if a<0:
                    lo,hi=hi,lo
                left+=a*lo
                right+=a*hi
            row_error+=max(abs(left),abs(right))
        largest=max(largest,row_error)
        trial.append(sum(a*y for a,y in zip(row,target)))
    assert largest < den
    eps=Q(largest,den)
    change=eps*max(abs(Q(v,trialden)) for v in trial)/(1-eps)
    assert min(Q(v,trialden) for v in trial)-change > Q(79,100000)
    print("Integer interval products prove that all 46 exact rational weights exceed 79/100000.")
    # Obtain the inherited 96-bit Gram and Schur pivot from the SAME exact levels.
    assert record["levels"][0] == 1 << record["bits"]
    levels=[Q(x,1<<record["bits"]) for x in record["levels"]]
    for j in range(1,4):
        for i in range(j):
            raw=(i*levels[i]+j*levels[j])/(i+j)
            k=record["gram_bins"][j][i]
            corr=raw*raw/(levels[i]*levels[j])
            assert Q(k,2**96)**2 < corr < Q(k+1,2**96)**2
    pivot=run_path(str(root/"ch39-finite-rh-pivot-mass.py"))
    assert pivot["bins"] == record["gram_bins"]
    pivot["main"]()


if __name__=="__main__":
    main()

