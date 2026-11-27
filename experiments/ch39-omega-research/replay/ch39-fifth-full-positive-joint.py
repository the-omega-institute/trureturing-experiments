"""第五壳完整夸克、第三枢轴五式与连续能标正动能的同一第四步。"""
from pathlib import Path
from runpy import run_path

NAME="ch39-fifth-full-positive-joint"
FRAME=((3,0,"i"),(3,0,"r"),(3,1,"i"),(2,0,"r"),(2,1,"r"),(2,0,"i"))
COEFFICIENTS=(303683517812,-42243366948,-45877802961,-296083391333,33782979269888,90145683704739,
              27285170980,-1772989573,7583036013,-23404272160,14764351651298,381538561071)


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    original=lib("fifth-third-pivot-three-high-positive")
    _,_,_,_,boxes,_,_=lib("fifth-third-pivot-em-colour-path")["solve"](
        root,original["CENTERS"],original["AXES"],tighten=True)
    directions=[(int(i),int(j),part) for i,j,part in
                (line.split() for line in (root/(NAME+"-directions.tsv")).read_text().splitlines())]
    rows=[i for i in range(19) if i not in (2,3,14)]
    assert len(set(directions))==len(directions)==78
    def frame(root, lib, blocks):
        return lib("fifth-free-frame-blocks")["solve"](root,lib,blocks,NAME,COEFFICIENTS,FRAME)
    sites=lib("fifth-flavor-cutoff-chart")["solve"](
        root,NAME,(0,0,0,1),directions,frame=frame,rows=rows,
        target_sites=[boxes[i] for i in rows],seed_rays=True,high_targets=[(i,boxes[i]) for i in (2,3)])
    assert all(s.overlaps(b) for s,b in zip(sites,boxes,strict=True))
    # 前十八个截止等式由联合方程精确强制。第 15 层由两条荷六重列的归一化与源迹相同精确给出。
    original["verify"](root)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
