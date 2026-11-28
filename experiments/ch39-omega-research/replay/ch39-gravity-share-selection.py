"""比较源布局中四条原有双支点线对引力双矩迁移的灵敏度与精确份额。"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parent


def verify():
    data = run_path(str(ROOT/"ch39-one-share-gravity-repair.py"))
    numbers = (*data["HOLDOUT"],*data["FORWARD"])
    values = {}
    for line in (16,285,301,323):
        with redirect_stdout(StringIO()):
            slope,corrections,error = data["verify"](numbers,line=line)
        values[line] = slope,corrections
        print("源线、沿牛顿等式的真空导数、最大份额变化与预测相对误差",line,slope.str(12),
              max(abs(v).upper() for v in corrections.values()).str(12),error.str(12),flush=True)
    order = sorted(values,key=lambda i:float(abs(values[i][0])),reverse=True)
    assert order==[16,301,285,323]
    assert all(abs(values[i][0])>abs(values[j][0]) for i,j in zip(order,order[1:]))
    for target in numbers:
        assert all(abs(values[i][1][target])<abs(values[j][1][target]) for i,j in zip(order,order[1:]))
    print("源绝对导数由大到小、二十四个目标所需份额绝对值由小到大的共同排序",order)


if __name__=="__main__":
    verify()
