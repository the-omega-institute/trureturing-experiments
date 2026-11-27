"""保持七维平凡字符子空间的连续五壳十五式共同见证。"""
from runpy import run_path
from pathlib import Path

from flint import arb

PARAMETERS = (
    "-0.625", "0.125", "1.625", "1.375",
    "0.25", "-0.75", "-1.5", "-0",
    "-0.0198041831456782", "1.875", "-0.25", "-2.25",
    "-1.875", "0.25", "1.75", "0.875",
    "1.25", "0.875", "1.125", "-0.25",
    "0.5", "0.75", "-1.00118244060565", "-0.25",
    "0.125", "-0.125", "0.625", "-0",
    "0.56202075843254", "-1.5", "3.125", "0.75",
    "1.49490555444976", "0", "-0.25", "0",
    "0", "0.625", "-0", "0",
    "-0.17086541468578", "-0.625", "-2.43416080543375", "-1.125",
    "-0.75", "1.125", "0.625", "0.625",
    "0.375", "-1.75", "-0.375", "-0.25",
    "-0.696255916609773", "-0.625", "-2", "1.25",
    "0.910987701666681", "2.5", "-0.125", "0.625",
    "-0.75", "-0.875", "0.25", "0.875",
    "-0.413191542285008", "-0.875", "1.5", "1.125",
    "1.06819633140711", "2.125", "0.125", "0.5",
    "0.125", "-0.25", "-0.125", "-1",
)
ACTIVE = (8, 22, 28, 32, 40, 42, 52, 56, 64, 68)
LINES = (0,0,0,0,15)
ENERGIES = ("522.398135900","534.208610459","544.705683199","544.755583188","919.160621983")


def experiment(root,fifth,columns,levels,*,return_chart=False):
    base = run_path(str(root/"ch39-cyclic-two-characters-nine.py"))
    chart = base["experiment"](root,fifth,columns,levels,return_chart=True,
                              parameters=PARAMETERS,active=ACTIVE,lines=LINES,energies=ENERGIES,
                              uniform_grid=1,ratio_bound=".525",shared_source_error=True,
                              merge_nontrivial=True,align_source=True)

    forms,blocks,*_ = chart
    # 原第二非平凡字符是合并块的后六个坐标。单独保持该字符的酉会保留此迹。
    source = base["adj"](forms[1])*(fifth*columns)
    original_trace, image_trace = [
        sum(abs(block[i,j])**2 for i in range(6,12) for j in range(2))
        for block in (source,blocks[1])
    ]
    print("第二非平凡字符的源、像 Gram 迹",original_trace.str(12),image_trace.str(12),flush=True)
    assert original_trace<arb(".000001") and image_trace>arb(".35")
    if return_chart:
        return chart
    return sum((form*block for form,block in zip(forms[1:],blocks[1:],strict=True)),forms[0]*blocks[0]),levels


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    joint=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*joint["verify"](root,return_operator="spectral"))
