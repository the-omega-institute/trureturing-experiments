"""重放第八壳预登记的三个支撑指标；可传入冻结提交的 replay 目录。"""
from pathlib import Path
from runpy import run_path
import sys

def verify(root):
    evaluate = run_path(str(root / "ch39-cross-shell-pivot-support.py"))["evaluate"]
    scores = evaluate(root, 2_306_888_100_093)
    for pivot, (_, _, value) in enumerate(scores, 1):
        expected_positive = pivot in (1, 2)
        if value > 0:
            outcome = "命中" if expected_positive else "失配"
        elif value <= 0:
            outcome = "失配" if expected_positive else "命中"
        else:
            outcome = "未决"
        print("枢轴、S3下端、S3上端、预登记判定", pivot, value.lower().str(25),
              value.upper().str(25), outcome, flush=True)
        if value > 0:
            print("三个规范端点成立时的色逆耦合上界", (-10*value).upper().str(25), flush=True)
    return scores


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent)
