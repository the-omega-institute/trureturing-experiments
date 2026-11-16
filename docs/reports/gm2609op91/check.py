#!/usr/bin/env python3
"""Exact operational checks for arXiv:2609.11326v3, Open Problem 9.1.

Finite tests of compiler formulas, not a decision procedure for disruption and
not a proof of either non-semidecidability assertion. All evidence uses integers.
The register-machine fixtures have proved-by-inspection halt/cycle behaviour;
cycle detection below independently checks their execution on self input.
"""
from dataclasses import dataclass
from itertools import product

# A deterministic natural-register machine, with an explicit program counter.
# Instructions are (op, argument); HALT returns the register. A machine step
# includes HALT, so a halting time is always a positive integer.
EVENTS = (
    (("jmp", 0),),
    (("halt", 0),),
    (("inc", 0), ("halt", 0)),
    (("inc", 0), ("inc", 0), ("halt", 0)),
    (("jmp", 1), ("jmp", 0)),
    (("set", 0), ("jz", 3), ("jmp", 0), ("halt", 0)),
    tuple([("inc", 0)] * 7 + [("halt", 0)]),
)
EXPECTED_TIMES = (None, 1, 2, 3, None, 3, 8)


def step(code, pc, reg):
    op, arg = code[pc]
    if op == "halt":
        return None
    if op == "jmp":
        return arg, reg
    if op == "inc":
        return pc + 1, reg + 1
    if op == "set":
        return pc + 1, arg
    if op == "jz":
        return (arg if reg == 0 else pc + 1), reg
    raise AssertionError(op)


def halt_time(e):
    """Exact only for these finite fixtures: halt, or a repeated full state."""
    code, pc, reg, ticks, seen = EVENTS[e], 0, e, 0, set()
    while (pc, reg) not in seen:
        assert ticks < 64, "fixture not resolved by halt or a repeated state"
        seen.add((pc, reg))
        ticks += 1
        nxt = step(code, pc, reg)
        if nxt is None:
            return ticks
        pc, reg = nxt
    return None


def halted_within(e, budget):
    """Bounded simulation; a repeated state safely shortens a pure cycle."""
    code, pc, reg, seen = EVENTS[e], 0, e, set()
    for _ in range(budget):
        if (pc, reg) in seen:
            return False
        seen.add((pc, reg))
        nxt = step(code, pc, reg)
        if nxt is None:
            return True
        pc, reg = nxt
    return False


BEHAVIOURS = ("bottom", "identity", "zero", "one", "successor")


@dataclass(frozen=True)
class Program:
    op: str
    args: tuple = ()


def base(name):
    return Program(name)


def gate(e, target):
    return Program("gate", (e, target))


def evaluate(program, z):
    """Small-step interpreter for bases, inert wrappers and nested gates.

    A repeated complete interpreter state detects divergence in this fixture
    language. This does not solve halting for arbitrary register machines.
    """
    current, continuation, seen = program, (), set()
    for _ in range(256):
        state = (current, continuation, z)
        if state in seen:
            return None
        seen.add(state)
        op, args = current.op, current.args
        if op == "wrap":
            current = args[0]
        elif op == "gate":
            e, target = args
            continuation = (target,) + continuation
            current = Program("event", (e, 0, e))
        elif op == "event":
            e, pc, reg = args
            nxt = step(EVENTS[e], pc, reg)
            if nxt is None:
                assert continuation
                current, continuation = continuation[0], continuation[1:]
            else:
                current = Program("event", (e, *nxt))
        elif op == "bottom":
            current = Program("bottom")
        elif op == "identity":
            assert not continuation
            return z
        elif op == "zero":
            assert not continuation
            return 0
        elif op == "one":
            assert not continuation
            return 1
        elif op == "successor":
            assert not continuation
            return z + 1
        else:
            raise AssertionError(op)
    raise AssertionError("fixture exceeded interpreter bound")


def wrapper_index(t, e):
    # Disjoint finite data codes, with zero reserved for a non-wrapper.
    return 1 + len(EVENTS) * BEHAVIOURS.index(t) + e


def decode_wrapper(x):
    if 1 <= x <= len(EVENTS) * len(BEHAVIOURS):
        q, e = divmod(x - 1, len(EVENTS))
        return BEHAVIOURS[q], e
    return None


def delta(x):
    decoded = decode_wrapper(x)
    return decoded[1] if decoded is not None else 0


def disruptive(target, x):
    return gate(delta(x), base(target))


def compiler_A(e, target, x):
    # The compiler returns syntax immediately. The unbounded wait is run time.
    return gate(e, disruptive(target, x))


def compiler_B(e, target, x):
    # Exactly x simulated steps at compile time, or a detected finite cycle.
    return base("zero") if halted_within(e, x) else disruptive(target, x)


def name_of_output(target, x, times):
    return target if times[delta(x)] is not None else "bottom"


def expected_value(name, z):
    return {"bottom": None, "identity": z, "zero": 0,
            "one": 1, "successor": z + 1}[name]


def main():
    times = tuple(halt_time(e) for e in range(len(EVENTS)))
    assert times == EXPECTED_TIMES
    bounded_checks = 0
    for e, budget in product(range(len(EVENTS)), range(41)):
        assert halted_within(e, budget) == (times[e] is not None and times[e] <= budget)
        bounded_checks += 1
    print("event halt times:", times)
    print("bounded-simulation comparisons:", bounded_checks)

    # Fidelity: source Section 3.2, identity property and an inert wrapper.
    fidelity_checks = 0
    for e, z in product(range(len(EVENTS)), range(9)):
        inert = Program("wrap", (base("identity"), e))
        assert evaluate(inert, z) == z
        out = disruptive("identity", wrapper_index("identity", e))
        assert evaluate(out, z) == (z if times[e] is not None else None)
        fidelity_checks += 1
    print("source identity/bottom instrumentation cases:", fidelity_checks)

    # Test every nontrivial property of five distinct computable behaviours,
    # every possible accepted base, both bottom-membership cases, and all
    # fixture source/event inputs. This checks the formulas independently of
    # any decidability assumption about an arbitrary infinite property.
    property_count = base_cases = output_checks = tail_checks = 0
    inputs = range(48)
    for mask in range(1, (1 << len(BEHAVIOURS)) - 1):
        prop = {b for i, b in enumerate(BEHAVIOURS) if mask & (1 << i)}
        property_count += 1
        target = (next(b for b in BEHAVIOURS if b not in prop)
                  if "bottom" in prop else next(b for b in BEHAVIOURS if b in prop))
        positive = "bottom" not in prop
        for t in prop:
            base_cases += 1
            for e in range(len(EVENTS)):
                x = wrapper_index(t, e)
                assert delta(x) == e
                out_name = name_of_output(target, x, times)
                assert (out_name in prop) == ((times[e] is not None) == positive)
                for z in range(5):
                    assert evaluate(Program("wrap", (base(t), e)), z) == expected_value(t, z)
                    assert evaluate(disruptive(target, x), z) == expected_value(out_name, z)
                    output_checks += 1
        for e, x, z in product(range(len(EVENTS)), inputs, range(5)):
            d_name = name_of_output(target, x, times)
            a_name = d_name if times[e] is not None else "bottom"
            b_name = ("zero" if times[e] is not None and times[e] <= x else d_name)
            assert evaluate(compiler_A(e, target, x), z) == expected_value(a_name, z)
            assert evaluate(compiler_B(e, target, x), z) == expected_value(b_name, z)
            output_checks += 2
        for e, h in enumerate(times):
            if h is None:
                for x in inputs:
                    assert compiler_B(e, target, x) == disruptive(target, x)
            else:
                finite_table = {x: name_of_output(target, x, times) in prop for x in range(h)}
                for x in inputs:
                    expected_member = finite_table[x] if x < h else "zero" in prop
                    actual_name = "zero" if halted_within(e, x) else name_of_output(target, x, times)
                    assert (actual_name in prop) == expected_member
                    if x >= h:
                        assert compiler_B(e, target, x) == base("zero")
                        tail_checks += 1
    print("nontrivial finite behaviour properties:", property_count)
    print("accepted-base cases:", base_cases)
    print("operational output comparisons:", output_checks)
    print("constant-tail comparisons:", tail_checks)
    print("scope: finite compiler/fidelity checks; infinite undecidability requires the proposed argument")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
