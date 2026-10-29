#!/usr/bin/env python3
"""Exact diagnostics for KPW Problem 7.3; no Lean settlement is asserted.

Run with Python 3 and the standard library only. All checks use finite sets,
integer counts and exact graph reachability. The infinite mathematical argument
is stated separately in the preregistration, not inferred from a cutoff.
"""

from collections import Counter, deque
from functools import lru_cache
from itertools import islice, permutations, product


PAYLOADS = "abcde"
VERTICES = range(6)
# (source, label, target). The order also fixes the annotation numbering.
EDGES = (
    (4, "a", 4), (4, "b", 4),
    (5, "a", 5), (5, "b", 5),
    (0, "c", 1), (0, "c", 5),
    (1, "a", 2), (2, "b", 1),
    (1, "d", 3), (4, "d", 3),
) + tuple(sorted({(s, "e", 0) for s in VERTICES}
                 | {(0, "e", t) for t in VERTICES}
                 | {(v, "e", v) for v in range(1, 5)}))
# Symbols 0..24 are annotated edges; 25..29 are plain payload symbols.
SYMBOLS = range(30)


def payload(symbol):
    return EDGES[symbol][1] if symbol < 25 else PAYLOADS[symbol - 25]


@lru_cache(None)
def choices(symbol):
    return (symbol,) if symbol < 25 else tuple(
        j for j, edge in enumerate(EDGES) if edge[1] == payload(symbol))


@lru_cache(None)
def realize(symbols):
    """One compatible underlying edge path, or None, by independent DP."""
    states = {v: () for v in VERTICES}
    for symbol in symbols:
        nxt = {}
        for j in choices(symbol):
            s, _, t = EDGES[j]
            if s in states:
                nxt.setdefault(t, states[s] + (j,))
        states = nxt
    return next(iter(states.values()), None)


def plain_word(word):
    return tuple(25 + PAYLOADS.index(c) for c in word)


def local_legal(symbols):
    """Finite forbidden rules of the 30-symbol SFT, window size at most 4."""
    return not any(all(s >= 25 for s in symbols[i:i + 3])
                   for i in range(len(symbols) - 2)) and all(
        realize(symbols[i:i + n]) is not None
        for n in range(1, 5) for i in range(len(symbols) - n + 1))


def reachable(start, edges):
    seen, todo = {start}, [start]
    while todo:
        s = todo.pop()
        for a, _, b in edges:
            if a == s and b not in seen:
                seen.add(b)
                todo.append(b)
    return seen


def source_fidelity():
    # Section 6: the ternary example's graph is exactly the stated 2-SFT.
    source_edges = {("a", "b"), ("b", "c"), ("c", "a"), ("c", "c")}
    forbidden = {"aa", "ac", "ba", "bb", "cb"}
    assert {a + b for a, b in product("abc", repeat=2)
            if (a, b) not in source_edges} == forbidden
    assert all(any(a == v for a, _ in source_edges) and
               any(b == v for _, b in source_edges) for v in "abc")
    # Thus all four allowed edges have bi-infinite extensions. This is an
    # exhaustive reproduction of the source's finite graph presentation.
    print("FIDELITY ternary SFT: allowed={ab,bc,ca,cc}; forbidden={aa,ac,ba,bb,cb}")

    # Golden mean example: compare the local source test to the all-subfactor
    # Parikh test on every binary word of length at most 12.
    count = 0
    for n in range(13):
        attainable = {sum(w) for w in product((0, 1), repeat=n)
                      if all(w[j:j + 2] != (1, 1) for j in range(n - 1))}
        assert attainable == set(range((n + 1) // 2 + 1))
        for w in product((0, 1), repeat=n):
            local = all(w[j:j + 2] != (1, 1) for j in range(n - 1))
            abelian = all(sum(w[i:j]) <= (j - i + 1) // 2
                          for i in range(n + 1) for j in range(i, n + 1))
            assert local == abelian
            count += 1
    print(f"FIDELITY golden mean: all {count} words of lengths 0..12 pass equality")


def graph_certificates():
    assert len(EDGES) == len(set(EDGES)) == 25
    assert all(reachable(v, EDGES) == set(VERTICES) for v in VERTICES)
    no_reset = tuple(edge for edge in EDGES if edge[1] != "e")
    # The states on any segment after c and before d are exactly {1,2}.
    after_c = {t for _, label, t in no_reset if label == "c"}
    before_d = {s for s, label, _ in no_reset if label == "d"}
    forward = set().union(*(reachable(v, no_reset) for v in after_c))
    backward = {v for v in VERTICES
                if reachable(v, no_reset) & before_d}
    core = forward & backward
    assert core == {1, 2}
    potential = {1: 0, 2: 1}
    weight = {"a": 1, "b": -1, "c": 0, "d": 0}
    for s, c, t in no_reset:
        if s in core and t in core:
            assert weight[c] == potential[t] - potential[s]
    assert {(s, t) for s, c, t in no_reset
            if c == "c" and t in backward} == {(0, 1)}
    assert {(s, t) for s, c, t in no_reset
            if c == "d" and s in forward} == {(1, 3)}
    assert not reachable(3, no_reset) & {0, 1}
    print("GRAPH vertices=6 edges=25 strongly_connected=yes")
    print("GRAPH no-reset c-to-d core={1,2}; potential={1:0,2:1}; a-b telescopes to 0")


def check_lift_memory():
    # Exhaust the finite combined automaton: suffix of length at most 3,
    # together with the end states reachable through the WHOLE prefix.
    # This detects a locally valid but globally unrealizable word of ANY
    # length, rather than checking only windows whose size defines the SFT.
    initial = ((), frozenset(VERTICES))
    seen, todo, transitions = {initial}, deque([initial]), 0
    while todo:
        suffix, ends = todo.popleft()
        for symbol in SYMBOLS:
            window = suffix + (symbol,)
            if not local_legal(window):
                continue
            new_ends = frozenset(EDGES[j][2] for j in choices(symbol)
                                 if EDGES[j][0] in ends)
            assert new_ends, (suffix, ends, symbol)
            nxt = (window[-3:], new_ends)
            transitions += 1
            if nxt not in seen:
                seen.add(nxt)
                todo.append(nxt)
    print(f"LIFT exact finite memory: reachable_states={len(seen)} transitions={transitions}; no local/global mismatch")


def representative(word):
    """Preserve counts and both endpoint letters of a target block hull."""
    if "c" in word and "d" in word:
        i, j = word.index("c"), word.index("d")
        middle = word[i + 1:j]
        assert middle.count("a") == middle.count("b")
        result = word[:i + 1] + "ab" * middle.count("a") + word[j:]
    else:
        result = word
    assert Counter(word) == Counter(result)
    assert word[0] == result[0] and word[-1] == result[-1]
    path = realize(plain_word(result))
    assert path is not None, result
    # Preserve the canonical code at the first and last block; annotations
    # supply all intervening hidden edges. At most two consecutive plains.
    decorated = tuple(25 + PAYLOADS.index(c) if i in (0, len(result) - 1)
                      else path[i] for i, c in enumerate(result))
    assert local_legal(decorated)
    assert realize(decorated) is not None
    # Extend the actual edge path on both sides using annotated e edges.
    s, t = EDGES[path[0]][0], EDGES[path[-1]][2]
    reset_in = EDGES.index((0, "e", s))
    reset_out = EDGES.index((t, "e", 0))
    reset_loop = EDGES.index((0, "e", 0))
    extended = (reset_loop, reset_in) + decorated + (reset_out, reset_loop)
    assert local_legal(extended)
    return decorated


def codebook(r):
    identity = tuple(range(r))
    alternatives = list(islice((p for p in permutations(range(r)) if p != identity), 25))
    assert len(alternatives) == 25
    result = {s: ("#", payload(s)) + tuple(f"t{j}" for j in (
        alternatives[s] if s < 25 else identity)) for s in SYMBOLS}
    assert len(set(result.values())) == 30
    assert all(len(c) == r + 2 and c.count("#") == 1 and c[0] == "#"
               and Counter(c[2:]) == Counter(f"t{j}" for j in range(r))
               for c in result.values())
    return result


def encode(symbols, book):
    return tuple(c for symbol in symbols for c in book[symbol])


def trim_check():
    hulls, cuts = 0, 0
    for r in range(5, 9):
        book, q = codebook(r), r + 2
        for k in range(9):
            base = "eec" + "a" * k + "b" * k + "dee"
            for i in range(len(base)):
                for j in range(i + 1, len(base) + 1):
                    word = base[i:j]
                    decorated = representative(word)
                    canonical = encode(plain_word(word), book)
                    actual = encode(decorated, book)
                    assert Counter(canonical) == Counter(actual)
                    for left in range(q):
                        for right in range(q):
                            stop = len(canonical) - right
                            if left > stop:
                                continue
                            assert Counter(canonical[left:stop]) == Counter(actual[left:stop]), (r, k, word, left, right)
                            cuts += 1
                    hulls += 1
        print(f"CODE r={r}: alphabet={r + 6}; symbols=30; block_length={q}; max_forbidden_length={4*q}")
    print(f"TRIMMING k=0..8 r=5..8: hulls={hulls} partial-block cuts={cuts}; all endpoint-preserving witnesses valid")


@lru_cache(None)
def parikh_path_from(state, counts):
    if not any(counts):
        return True
    for s, c, t in EDGES:
        if s == state and c != "e":
            j = "abcd".index(c)
            if counts[j]:
                rem = list(counts)
                rem[j] -= 1
                if parikh_path_from(t, tuple(rem)):
                    return True
    return False


def obstruction_checks():
    count = 0
    for m in range(13):
        for n in range(13):
            possible = any(parikh_path_from(v, (m, n, 1, 1)) for v in VERTICES)
            assert possible == (m == n), (m, n)
            count += 1
    # Concrete pairwise distinguishing continuations: prefix c a^i,
    # continuation b^i d. Each distinct prefix j is rejected by the invariant.
    for i, j in product(range(13), repeat=2):
        if i != j:
            assert any(parikh_path_from(v, (i, i, 1, 1)) for v in VERTICES)
            assert not any(parikh_path_from(v, (j, i, 1, 1)) for v in VERTICES)
    print(f"OBSTRUCTION exact Parikh DP: m,n=0..12 cases={count}; possible iff m=n")
    print("DISTINGUISHING 13 prefixes: 156 ordered pairs separated by b^i d")
    for n in range(13):
        paths = {tuple(EDGES.index((4, c, 4)) for c in w)
                 for w in product("ab", repeat=n)}
        assert len(paths) == 2 ** n
    print("GROWTH annotated loops at 4: exactly 2^n distinct code words for n=0..12")


def main():
    source_fidelity()
    graph_certificates()
    check_lift_memory()
    trim_check()
    obstruction_checks()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
