"""同六个目标的三种规则；复用已公布的区间，只重算尚缺的九格。"""
from argparse import ArgumentParser
from decimal import Decimal
from pathlib import Path
import re
from runpy import run_path


SCALE = 10**7
TOLERANCE = 15*SCALE//100
TARGETS = (11, 13, 15, 17, 19, 21)
RULES = ("照搬", "单钟", "双钟")


def published(root, name, prefix):
    lines = [row for row in (root/name).read_text().splitlines() if row.startswith(prefix)]
    assert len(lines) == 1
    pairs = re.findall(r"\[(-?\d+\.\d{7}),\s*(-?\d+\.\d{7})\]", lines[0])
    return [tuple(int(Decimal(n)*SCALE) for n in pair) for pair in pairs]


def previous_results(root, target):
    first = "三个相邻零点对的804边时钟首次完整结果.md"
    second = "三个相邻零点对的双时钟首次完整结果.md"
    if target <= 15:
        forecast, actual, _ = published(root, first, f"| {target} |")
        known = {"单钟": (forecast, "事前")}
    else:
        forecast, actual, _ = published(root, second, f"| {target-1}→{target} |")
        known = {"双钟": (forecast, "事前")}
    if target == 11:
        carry, = published(root, "第十壳单时钟的首次预测与完整读数.md", "| 完整第十壳 |")
        known["照搬"] = (carry, "已有源壳")
    if target == 13:
        swaps = "第十三壳单804时钟失配的四项交换实验.md"
        carry, _ = published(root, swaps, "| 十二 | 十二 |")
        dual, _ = published(root, swaps, "| 十三 | 十三 |")
        known["照搬"] = (carry, "已有交换")
        known["双钟"] = (dual, "已有交换")
    return actual, known


def sign(pair):
    low, high = pair
    return 1 if low > 0 else -1 if high < 0 else 0


def evaluation(forecast, actual):
    p, q = forecast
    a, b = actual
    low, high = max(0, p-b, a-q), max(abs(p-b), abs(q-a))
    ps, actual_sign = sign(forecast), sign(actual)
    same = "未定" if 0 in (ps, actual_sign) else "命中" if ps == actual_sign else "失配"
    tolerance = "命中" if high < TOLERANCE else "失配" if low >= TOLERANCE else "未定"
    return (low, high), same, tolerance


def compare(root, target):
    measure = run_path(str(root/"ch39-tenth-804-clock-holdout.py"))["score_bounds"]
    actual, known = previous_results(root, target)
    results = {}
    for rule in RULES:
        if rule not in known:
            if rule == "照搬":
                lo, hi, scale = measure(root, "actual", previous_shell=target-2)
            else:
                lo, hi, scale = measure(root, "forecast", previous_shell=target-1, update_724=(rule == "双钟"))
            assert scale == SCALE
            known[rule] = ((lo, hi), "本次回算")
        forecast, origin = known[rule]
        error, same, tolerance = evaluation(forecast, actual)
        print("对照", target, rule, *actual, *forecast, *error, same, tolerance, origin, flush=True)
        results[rule] = (error, same, tolerance)
    return results


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("targets", type=int, nargs="*", help="省略时计算全部六壳")
    targets = cli.parse_args().targets or TARGETS
    assert all(t in TARGETS for t in targets)
    root = Path(__file__).resolve().parent
    rows = [compare(root, target) for target in targets]
    for rule in RULES:
        data = [row[rule] for row in rows]
        signs = sum(s == "命中" for _,s,_ in data)
        passes = sum(t == "命中" for _,_,t in data)
        joint = sum(s == t == "命中" for _,s,t in data)
        n = len(data)
        mae_low = sum(e[0] for e,_,_ in data)//n
        mae_high = (sum(e[1] for e,_,_ in data)+n-1)//n
        print("规则、总数、符号、容差、同时命中、平均绝对误差两界分子、分母",
              rule,n,signs,passes,joint,mae_low,mae_high,SCALE)
