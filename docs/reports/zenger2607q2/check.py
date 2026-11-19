"""Exact finite-model checks for Intuitionistic Dynamic Logic, Question 2."""

from functools import cache
from itertools import product


BOTTOM = ("bottom",)
P = ("p",)
Q = ("q",)
ANTECEDENT = ("and", Q, ("master", ("imp", Q, ("and", P, ("box", ("box", Q))))))


def reachable(successors, start):
    seen = {start}
    pending = [start]
    while pending:
        for end in successors[pending.pop()] - seen:
            seen.add(end)
            pending.append(end)
    return seen


def evaluator(order, modal, valuation):
    """Return formula truth sets using the two relations in the thesis."""
    worlds = range(len(order))
    universe = frozenset(worlds)
    steps = [set().union(*(modal[v] for v in order[w])) for w in worlds]
    closure = [reachable(steps, w) for w in worlds]

    @cache
    def truth(formula):
        op, *args = formula
        match op:
            case "bottom":
                return frozenset()
            case "p" | "q":
                return valuation[op]
            case "and":
                return truth(args[0]) & truth(args[1])
            case "or":
                return truth(args[0]) | truth(args[1])
            case "imp":
                good = (universe - truth(args[0])) | truth(args[1])
                return frozenset(w for w in worlds if order[w] <= good)
            case "box" | "master":
                relation = steps if op == "box" else closure
                good = truth(args[0])
                return frozenset(w for w in worlds if relation[w] <= good)
            case _:
                raise ValueError(f"Unknown connective: {op}")

    return truth


def repeat_box(count, formula):
    for _ in range(count):
        formula = ("box", formula)
    return formula


def thesis_example():
    # Worlds in Figure 3.1: w, v, f(w), f(v).
    order = [{0, 1}, {1}, {2}, {3}]
    modal = [{2}, {3}, {2}, {3}]
    truth = evaluator(order, modal, {"p": frozenset({2}), "q": frozenset({3})})
    lhs = ("box", ("or", P, Q))
    rhs = ("or", ("box", P), ("box", Q))
    assert 0 in truth(lhs)
    assert 0 not in truth(("box", P)) and 0 not in truth(("box", Q))
    assert 0 not in truth(("imp", lhs, rhs))
    print("Fidelity: Lemma 3.10, Figure 3.1 falsified at w as stated")


def all_rows(n):
    subsets = [
        frozenset(v for v in range(n) if mask & (1 << v))
        for mask in range(1 << n)
    ]
    return subsets, product(subsets, repeat=n)


def check_even_consequences():
    count = 0
    for n in (1, 2):
        worlds = range(n)
        subsets, orders = all_rows(n)
        for order in orders:
            if not all(
                w in order[w]
                and all(order[v] <= order[w] for v in order[w])
                and all(v == w or w not in order[v] for v in order[w])
                for w in worlds
            ):
                continue
            persistent_sets = [
                present for present in subsets
                if all(order[w] <= present for w in present)
            ]
            for modal in product(subsets, repeat=n):
                for p, q in product(persistent_sets, repeat=2):
                    truth = evaluator(order, modal, {"p": p, "q": q})
                    a_worlds = truth(ANTECEDENT)
                    for k in range(7):
                        assert a_worlds <= truth(repeat_box(2 * k, P))
                    count += 1
    print(f"Even-step consequences: {count} dynamic models, 0 <= k <= 6")


def chain(m, mark_even=False):
    worlds = range(m + 2)
    return evaluator(
        [{w} for w in worlds],
        [{min(w + 1, m + 1)} for w in worlds],
        {
            "p": frozenset(w for w in worlds if w != m),
            "q": frozenset(w for w in worlds if mark_even and w % 2 == 0),
        },
    )


def formulas(max_size):
    by_size = [[], [P, BOTTOM]]
    for size in range(2, max_size + 1):
        layer = [(op, a) for op in ("box", "master") for a in by_size[size - 1]]
        for left_size in range(1, size - 1):
            layer.extend(
                (op, a, b)
                for op in ("and", "or", "imp")
                for a in by_size[left_size]
                for b in by_size[size - 1 - left_size]
            )
        by_size.append(layer)
    return [formula for layer in by_size for formula in layer]


@cache
def height(formula):
    op, *args = formula
    match op:
        case "p":
            return 1
        case "bottom":
            return 0
        case "box":
            return 1 + height(args[0])
        case "master":
            return height(args[0])
        case "and" | "or" | "imp":
            return max(height(args[0]), height(args[1]))
        case _:
            raise ValueError(f"Unknown connective: {op}")


def check_family():
    max_defect = 17
    models = [chain(m) for m in range(max_defect + 1)]
    extensions = {m: chain(m, mark_even=True) for m in range(1, max_defect + 1, 2)}
    for m, truth in extensions.items():
        assert 0 in truth(ANTECEDENT)
        for k in range(0, max_defect + 1):
            assert 0 in truth(repeat_box(2 * k, P))
    for m in range(0, max_defect + 1, 2):
        assert 0 not in models[m](repeat_box(m, P))

    tested = formulas(max_size=7)
    transfers = 0
    for formula in tested:
        threshold = height(formula)
        truth_at_root = [0 in truth(formula) for truth in models]
        assert all(value == truth_at_root[threshold] for value in truth_at_root[threshold:])
        m = 2 * threshold + 2
        on_odd_extension = 0 in extensions[m + 1](formula)
        assert on_odd_extension == truth_at_root[m + 1] == truth_at_root[m]
        if on_odd_extension:
            assert 0 not in models[m](repeat_box(m, P))
            transfers += 1
    print(f"Parity models: odd A-witnesses and even countermodels for 0 <= m <= {max_defect}")
    print(f"Stabilization: {len(tested)} syntax trees, size <= 7, defects through {max_defect}")
    print(f"Paired countermodels: {transfers} transfers with the candidate true in the odd model")


def main():
    thesis_example()
    check_even_consequences()
    check_family()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
