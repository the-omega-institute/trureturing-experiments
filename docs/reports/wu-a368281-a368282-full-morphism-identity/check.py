#!/usr/bin/env python3
"""Exact finite checks for Wu's OEIS A368281/A368282 conjecture.

The word and recurrence are computed independently.  These experiments are
finite evidence, not a proof of the assertion for every parameter and index.
OEIS inputs retain CC BY-SA 4.0; see README.md for attribution and offsets.
"""

from hashlib import sha256
from itertools import accumulate
from pathlib import Path
import sys

DATA = Path(__file__).resolve().parent / "data"
LIMIT = 10_000
HASHES = {
    "b005206.txt": "965977a39f34fb90ac4a772b8e3ab77282803f7733f3678dd977b7889b1c3d22",
    "b005614.txt": "6b199f7a2b2e0c5b26641305e0f2e035c6dc25bb1c3f44ae79eb9e62fefbd1b2",
    "b285431.txt": "c5218c27af1e849fc48850247745f8d9f79664f28f4034d93f9079f714f66705",
    "b286389.txt": "241b74574cbf0b406fc5c40487b5961a5cf6a80db11887093d616ad9e340baae",
    "b366860.txt": "9358c67389f9f7eb86b04780293073b5dd5c4b166953dcff584a3175422b5992",
    "b366870.txt": "c273d2f847d88cfcbc85d896a8d644582a1868516367d07645bb21b373b6ae75",
    "b366871.txt": "7d8e7a9bc1235cd2703c6c9c5b6165da845c23048565fcb978ad4d954ffdc152",
    "b368281.txt": "2a54713c828dad37778461e1327142c53d2548b88f46fe6a5f5a081574c9ddd9",
    "b368282.txt": "f4030264f140817478aad2fec310e8d646c875fa604c0e737f9eeb89377e87d4",
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read_bfile(name):
    raw = (DATA / name).read_bytes()
    require(sha256(raw).hexdigest() == HASHES[name], f"input hash: {name}")
    result = {}
    for line in raw.decode("utf-8").splitlines():
        fields = line.split("#", 1)[0].split()
        if not fields:
            continue
        require(len(fields) == 2, f"malformed b-file line: {name}")
        index, value = map(int, fields)
        require(index not in result, f"duplicate index: {name}")
        result[index] = value
    return result


def substitute(m, word):
    one_block = [1] * m
    return [digit for letter in word for digit in
            (one_block + [0] if letter == 1 else one_block)]


def fixed_prefix(m, length):
    require(m >= 1, "denominator must be positive")
    word = [1]
    while len(word) < length:
        longer = substitute(m, word)
        require(len(longer) > len(word), f"growth at m={m}")
        require(longer[:len(word)] == word, f"nested prefix at m={m}")
        word = longer
    return word[:length]


def recursive_row(m, length):
    require(m >= 1, "denominator must be positive")
    row = [0]
    for k in range(1, length + 1):
        q = row[k - 1] // m
        require(0 <= q < k, f"recurrence causality at {(m, k)}")
        value = k - row[q]
        require(0 <= value <= k, f"recurrence bound at {(m, k)}")
        row.append(value)
    return row


def check_blocks(m, word, sums):
    """Locate actual output positions, without using recurrence values."""
    position = 0
    cases = [0, 0]
    for j, letter in enumerate(word, start=1):
        for t in range(1, m + 1):
            if position == len(word):
                return cases
            position += 1
            require(word[position - 1] == 1, f"block one at {(m, position)}")
            require(sums[position - 1] == m * (j - 1) + t - 1,
                    f"one-count formula at {(m, position)}")
            require(sums[position - 1] // m == j - 1,
                    f"one quotient at {(m, position)}")
            require(position - sums[position] == sums[j - 1],
                    f"one-case zero bridge at {(m, position)}")
            cases[0] += 1
        if letter == 1:
            if position == len(word):
                return cases
            position += 1
            require(word[position - 1] == 0, f"block zero at {(m, position)}")
            require(sums[position - 1] == m * j,
                    f"terminal one-count at {(m, position)}")
            require(sums[position - 1] // m == j,
                    f"terminal quotient at {(m, position)}")
            require(position - sums[position] == sums[j],
                    f"terminal zero bridge at {(m, position)}")
            cases[1] += 1
    require(position == len(word), "unfinished block traversal")
    return cases


def antidiagonal_cells(count, word_table):
    result = []
    diagonal = 1
    while len(result) < count:
        for m in range(1, diagonal + 1):
            k = diagonal - m + (1 if word_table else 0)
            result.append((m, k))
        diagonal += 1
    return result[:count]


def main():
    require(sys.version_info >= (3, 12), "repository requires Python >= 3.12")
    inputs = {name: read_bfile(name) for name in HASHES}
    words, rows, sums = {}, {}, {}
    one_cases = zero_cases = 0
    for m in range(1, 65):
        words[m] = fixed_prefix(m, LIMIT)
        rows[m] = recursive_row(m, LIMIT)
        sums[m] = [0, *accumulate(words[m])]
        require(rows[m] == sums[m], f"full identity through {LIMIT}, m={m}")
        require(substitute(m, words[m])[:LIMIT] == words[m],
                f"fixed-prefix fidelity, m={m}")
        c1, c0 = check_blocks(m, words[m], sums[m])
        one_cases += c1
        zero_cases += c0
        for k in range(1, LIMIT + 1):
            q = sums[m][k - 1] // m
            require(sums[m][k] == k - sums[m][q],
                    f"prefix-sum recurrence at {(m, k)}")
        if m <= 6:
            print(f"IDENTITY m={m} k=0..{LIMIT} S({LIMIT})={sums[m][-1]} OK")

    # Source Table 1, Shallit arXiv:2310.14252v2; OEIS k_0 is an adjoined zero.
    require(rows[2][:12] == [0, 1, 2, 2, 3, 4, 4, 5, 6, 7, 8, 8],
            "Shallit Table 1 g_n")
    require([0, *words[2][:11]] == [0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0],
            "Shallit Table 1 k_n")
    print("FIDELITY Shallit Table 1 g_0..g_11 and k_0..k_11 OK")

    for name, m, offset in [("b005614.txt", 1, 0), ("b285431.txt", 2, 1)]:
        count = 0
        for index, value in inputs[name].items():
            position = index - offset
            if 0 <= position < LIMIT:
                require(words[m][position] == value, f"word b-file {name}:{index}")
                count += 1
        print(f"BFILE word m={m} {name} comparisons={count} offset={offset} OK")

    for m, name in enumerate(["b005206.txt", "b286389.txt", "b366860.txt",
                              "b366870.txt", "b366871.txt"], start=1):
        require(set(range(LIMIT + 1)) <= inputs[name].keys(), f"coverage {name}")
        for k in range(LIMIT + 1):
            require(rows[m][k] == inputs[name][k], f"recurrence b-file {name}:{k}")
            require(sums[m][k] == inputs[name][k], f"partial-sum b-file {name}:{k}")
        print(f"BFILE recurrence and partial sums m={m} {name} k=0..{LIMIT} OK")

    for name, is_word in [("b368281.txt", True), ("b368282.txt", False)]:
        cells = antidiagonal_cells(len(inputs[name]), is_word)
        for index, (m, k) in enumerate(cells, start=1):
            expected = inputs[name][index]
            actual = words[m][k - 1] if is_word else rows[m][k]
            require(actual == expected, f"table b-file {name}:{index} cell {(m, k)}")
        row6 = [k for m, k in cells if m == 6]
        print(f"BFILE table {name} all={len(cells)} cells OK; m=6 k={min(row6)}..{max(row6)}")

    print("PUBLISHED_DATA_LIMIT table b-files have 87 entries each; no separate m=6")
    print("recurrence b-file located; m=6 k=0..10000 checked by independent constructions")
    print(f"SMALL_CASES m=1..64 k=0..{LIMIT} identities={64 * (LIMIT + 1)} OK")
    print(f"BLOCK_CASES ones={one_cases} terminal_zeros={zero_cases} OK")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
