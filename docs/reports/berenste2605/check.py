"""Exact Type A Demazure computations for arXiv:2605.09761v2, §3.1."""

from itertools import combinations


def require(condition, description):
    if not condition:
        raise AssertionError(description)


def identity(rank):
    return tuple(range(1, rank + 2))


def word_product(permutation, word):
    """Apply the paper's right Demazure rule for adjacent transpositions."""
    values = list(permutation)
    for i in word:
        if values[i - 1] < values[i]:
            values[i - 1], values[i] = values[i], values[i - 1]
    return tuple(values)


def longest_word(support):
    """A reduced reversal word for each maximal interval of the support."""
    word = []
    for end in sorted(support):
        start = end
        while start - 1 in support:
            start -= 1
        word.extend(range(end, start - 1, -1))
    return tuple(word)


def length(permutation):
    return sum(a > b for a, b in combinations(permutation, 2))


def orthogonal(left, right):
    return all(abs(i - j) > 1 for i in left for j in right)


def meets_claim(left, right):
    only_left, only_right, shared = left - right, right - left, left & right
    return left <= right or right <= left or all(
        orthogonal(a, b) for a, b in combinations((only_left, only_right, shared), 2)
    )


def alternating_length(rank, first, second):
    target = word_product(identity(rank), longest_word(first | second))
    words = (longest_word(first), longest_word(second))
    value = identity(rank)
    for count in range(1, length(target) + 2):
        value = word_product(value, words[(count - 1) % 2])
        if value == target:
            return count
    raise AssertionError("Alternating product did not reach the longest element")


def check_source_example():
    rank = 4
    left, right = frozenset({3, 4}), frozenset({2, 4})
    start = identity(rank)
    x, y = longest_word(left), longest_word(right)
    target = (1, 5, 4, 3, 2)
    require(word_product(start, x + y + x) == target, "Example 3.5: xyx")
    require(
        word_product(start, y + x + y) == word_product(start, (2, 3, 4, 3, 2)),
        "Example 3.5: yxy",
    )
    require(word_product(start, y + x + y + x) == target, "Example 3.5: yxyx")
    mu = alternating_length(rank, left, right), alternating_length(rank, right, left)
    require(mu == (3, 4), "Example 3.5: alternation minima")
    print(f"Example 3.5: mu(J,K), mu(K,J) = {mu}")


def check_family():
    for r in range(2, 11):
        rank = r + 2
        left, right = frozenset(range(1, r + 1)), frozenset({r, r + 2})
        start = identity(rank)
        x, y = longest_word(left), longest_word(right)
        xp, yp = word_product(start, x), word_product(start, y)
        common = tuple(range(r + 1, 0, -1)) + (r + 3, r + 2)
        require(xp == common[:-2] + (r + 2, r + 3), f"A{rank}: longest left")
        yp_expected = list(start)
        for i in sorted(right):
            yp_expected[i - 1], yp_expected[i] = yp_expected[i], yp_expected[i - 1]
        require(yp == tuple(yp_expected), f"A{rank}: longest right")
        require(length(xp) == len(x) and length(yp) == len(y), f"A{rank}: reduced words")
        require(word_product(xp, x) == xp and word_product(yp, y) == yp, f"A{rank}: idempotence")
        require(word_product(xp, y) == word_product(yp, x) == common, f"A{rank}: commutation")
        require(not meets_claim(left, right), f"A{rank}: support condition fails")
        mu = alternating_length(rank, left, right), alternating_length(rank, right, left)
        require(mu == (2, 2), f"A{rank}: alternation minima")
        print(f"A{rank}: r={r}, mu={mu}, common product={common}")


def check_small_ranks():
    for rank, expected in enumerate((0, 0, 0, 8, 64, 364), start=1):
        vertices = range(1, rank + 1)
        supports = [
            frozenset(subset) for size in range(rank + 1) for subset in combinations(vertices, size)
        ]
        words = {s: longest_word(s) for s in supports}
        values = {s: word_product(identity(rank), words[s]) for s in supports}
        count = sum(
            word_product(values[left], words[right]) == word_product(values[right], words[left])
            and not meets_claim(left, right)
            for left in supports
            for right in supports
        )
        require(count == expected, f"A{rank}: ordered-pair count {count}, expected {expected}")
        print(f"A{rank}: {count} violations among {len(supports) ** 2} ordered support pairs")


if __name__ == "__main__":
    check_source_example()
    check_family()
    check_small_ranks()
    print("ALL_CHECKS_OK")
