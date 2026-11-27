"""原九层见证的连续变形能否删去一个或两个低截止。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def probe_list():
    more=[(32000000,"y"),(100000000,"y"),(320000000,"y"),(1000000000,"y"),
          (300000,"col"),(100000,"col"),(25000,"col"),(5600,"y"),(10000,"y"),
          (3160000000,"y"),(42000,"col"),(56000000,"y"),(316000,"col"),(10000000000,"y")]
    shifted=[(p,ch,"-.22" if ch=="y" else ".58","0") for p,ch in (
        (32000000,"y"),(100000000,"y"),(320000000,"y"),(1000000000,"y"),
        (3160000000,"y"),(10000000000,"y"),(42000,"col"),(316000,"col"),
        (5600000,"col"),(100000,"col"),(750000,"col"),(17800000000,"y"))]
    return more,shifted


def verify(root):
    core=run_path(str(root/"ch39-fourth-third-pivot-six-low-path-obstruction.py"))
    source=run_path(str(root/"ch39-fourth-third-pivot-sparse-positive-path.py"))
    rows=[i+1 for i in source["ROWS"]]
    anchors=dict(zip(rows,map(arb,source["ANCHORS"]),strict=True))
    free={rows[i] for i in source["ACTIVE"]}
    assert free=={1,16,18}
    more,shifted=probe_list()
    lo,cut,constants,high_error,bounds=core["integral_table"](root,more,shifted)
    leaves=bounds.shape[1]//2
    log_span=(cut/lo).log()
    def grid(value,lower):
        pos=(value/lo).log()*leaves/log_span
        end=pos.lower() if lower else pos.upper()
        return max(0,min(leaves-1,int(end.floor().unique_fmpz())))
    # 六个移除位置允许所有余下低层自由变动，其余三个移除位置允许剩余锚独立变动十倍。
    domains=[]
    for removed in rows:
        support=[]
        for i in rows:
            if i==removed:
                continue
            low,high=0,leaves-1
            if removed in {3,14,15} and i not in free:
                low,high=grid(anchors[i]/10,True),grid(anchors[i]*10,False)
            support.append((i,low,high))
        domains.append(support)
    core["integer_cover"](root,constants,high_error,bounds,8,domains)
    # 七枚时，遍历十九层的所有支撑与物理允许的全部截止。
    core["integer_cover"](root,constants,high_error,bounds,7)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
