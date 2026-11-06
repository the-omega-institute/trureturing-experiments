#!/usr/bin/env python3
"""Exact finite-group check for the generalized Cayley regular-set claim.

The accompanying smallgroups.dat is emitted by GAP 4.16.1 + SmallGrp and contains
all 42 groups of order at most 16, their involutory automorphisms, and proper
subgroup masks.  This program performs the settlement test independently using
only integer arithmetic and finite set operations.
"""
from pathlib import Path


def parse_database(path):
    lines = path.read_text().splitlines()
    groups = []
    pos = 0
    while pos < len(lines):
        if not lines[pos].startswith("GROUP "):
            pos += 1
            continue
        _, order, index = lines[pos].split()
        order, index = int(order), int(index)
        pos += 1
        multiplication = []
        for _ in range(order):
            multiplication.append([int(x) - 1 for x in lines[pos].split()])
            pos += 1
        automorphisms = []
        subgroups = []
        while pos < len(lines) and not lines[pos].startswith("END"):
            if lines[pos].startswith("AUT "):
                automorphisms.append([int(x) - 1 for x in lines[pos].split()[1:]])
            elif lines[pos].startswith("SUB "):
                subgroups.append(int(lines[pos].split()[1]))
            pos += 1
        if pos == len(lines):
            raise ValueError("unterminated group record")
        pos += 1
        groups.append((order, index, multiplication, automorphisms, subgroups))
    return groups


def identity_and_inverses(multiplication):
    order = len(multiplication)
    identity = next(i for i, row in enumerate(multiplication) if row[i] == i)
    inverses = [
        next(j for j in range(order) if multiplication[i][j] == identity)
        for i in range(order)
    ]
    for i in range(order):
        if multiplication[i][inverses[i]] != identity:
            raise AssertionError("invalid inverse")
    return identity, inverses


def orbit_partition(mapping, elements):
    unseen = set(elements)
    result = []
    while unseen:
        start = min(unseen)
        orbit = []
        current = start
        while current in unseen:
            unseen.remove(current)
            orbit.append(current)
            current = mapping[current]
        result.append(orbit)
    return result


def is_subgroup(multiplication, identity, subset):
    return identity in subset and all(
        multiplication[x][y] in subset for x in subset for y in subset
    )


def right_cosets(multiplication, subgroup, order):
    seen = set()
    cosets = []
    for x in range(order):
        if x in subgroup:
            continue
        coset = frozenset(multiplication[h][x] for h in subgroup)
        if coset not in seen:
            seen.add(coset)
            cosets.append(coset)
    return cosets


def is_cyclic_subgroup(multiplication, subgroup, identity):
    for generator in subgroup:
        powers = set()
        current = identity
        while current not in powers:
            powers.add(current)
            current = multiplication[current][generator]
        if powers == subgroup:
            return True
    return False


def uniform_subset_exists(orbits, cosets, target):
    if target == 0:
        return True
    orbit_masks = []
    for orbit in orbits:
        mask = 0
        for element in orbit:
            mask |= 1 << element
        orbit_masks.append(mask)
    coset_masks = [sum(1 << element for element in coset) for coset in cosets]
    for choice in range(1 << len(orbit_masks)):
        selected = 0
        for bit, orbit_mask in enumerate(orbit_masks):
            if choice & (1 << bit):
                selected |= orbit_mask
        if all((selected & coset_mask).bit_count() == target for coset_mask in coset_masks):
            return True
    return False


def check_instance(order, multiplication, alpha, subgroup_mask):
    identity, inverses = identity_and_inverses(multiplication)
    subgroup = {i for i in range(order) if subgroup_mask & (1 << i)}
    omega_set = {
        multiplication[alpha[inverses[g]]][g] for g in range(order)
    }
    involution = [alpha[inverses[g]] for g in range(order)]

    admissible_inside = [g for g in subgroup if g not in omega_set]
    inside_orbits = orbit_partition(involution, admissible_inside)
    actual_a = {
        sum(len(inside_orbits[i]) for i in range(len(inside_orbits)) if choice & (1 << i))
        for choice in range(1 << len(inside_orbits))
    }
    fixed_admissible = any(involution[g] == g for g in admissible_inside)
    expected_a = {
        a for a in range(len(admissible_inside) + 1)
        if a % 2 == 0 or fixed_admissible
    }
    if actual_a != expected_a:
        return ("internal", actual_a, expected_a)

    outside_admissible = [
        g for g in range(order) if g not in subgroup and g not in omega_set
    ]
    outside_orbits = orbit_partition(involution, outside_admissible)
    cosets = right_cosets(multiplication, subgroup, order)
    capacities = [len(coset - omega_set) for coset in cosets]
    lmin = min(capacities)
    perfect_code = uniform_subset_exists(outside_orbits, cosets, 1)
    actual_b = set()
    orbit_masks = []
    for orbit in outside_orbits:
        mask = 0
        for element in orbit:
            mask |= 1 << element
        orbit_masks.append(mask)
    coset_masks = [sum(1 << element for element in coset) for coset in cosets]
    for choice in range(1 << len(orbit_masks)):
        selected = 0
        for bit, orbit_mask in enumerate(orbit_masks):
            if choice & (1 << bit):
                selected |= orbit_mask
        values = [(selected & coset_mask).bit_count() for coset_mask in coset_masks]
        if len(set(values)) == 1 and values[0] <= lmin:
            actual_b.add(values[0])
    expected_b = {
        b for b in range(lmin + 1)
        if b % 2 == 0 or perfect_code
    }
    if actual_b != expected_b:
        return ("external", actual_b, expected_b, lmin, perfect_code)
    return None


def main():
    groups = parse_database(Path(__file__).with_name("smallgroups.dat"))
    if len(groups) != 42:
        raise AssertionError(f"expected 42 groups, got {len(groups)}")
    total_automorphisms = sum(len(record[3]) for record in groups)
    if total_automorphisms != 740:
        raise AssertionError(f"unexpected involution count: {total_automorphisms}")

    checked = 0
    skipped_non_subgroup_w = 0
    trivial_h_cases = 0
    fixed_coset_cases = 0
    cyclic_index_two_cases = 0
    order_counts = {}
    failures = []

    for order, index, multiplication, automorphisms, subgroup_masks in groups:
        identity, inverses = identity_and_inverses(multiplication)
        for alpha in automorphisms:
            omega_set = {
                multiplication[alpha[inverses[g]]][g] for g in range(order)
            }
            if not is_subgroup(multiplication, identity, omega_set):
                skipped_non_subgroup_w += 1
                continue
            for subgroup_mask in subgroup_masks:
                subgroup = {i for i in range(order) if subgroup_mask & (1 << i)}
                if subgroup_mask == (1 << order) - 1:
                    continue
                if any(alpha[g] not in subgroup for g in subgroup):
                    continue
                checked += 1
                order_counts[order] = order_counts.get(order, 0) + 1
                if len(subgroup) == 1:
                    trivial_h_cases += 1
                involution = [alpha[inverses[g]] for g in range(order)]
                for x in range(order):
                    if x in subgroup:
                        continue
                    right = {multiplication[h][x] for h in subgroup}
                    image = {involution[g] for g in right}
                    if image == right:
                        fixed_coset_cases += 1
                        break
                failure = check_instance(order, multiplication, alpha, subgroup_mask)
                if failure is not None:
                    failures.append((order, index, failure))

            # Fidelity check for the source's cyclic-index-two family: every
            # alpha-invariant cyclic A of index two contains omega_alpha(G).
            for subgroup_mask in subgroup_masks:
                subgroup = {i for i in range(order) if subgroup_mask & (1 << i)}
                if len(subgroup) * 2 != order or len(subgroup) <= 1:
                    continue
                if any(alpha[g] not in subgroup for g in subgroup):
                    continue
                if not is_cyclic_subgroup(multiplication, subgroup, identity):
                    continue
                omega_set = {
                    multiplication[alpha[inverses[g]]][g] for g in range(order)
                }
                if not omega_set.issubset(subgroup):
                    failures.append((order, index, ("cyclic-index-two", omega_set, subgroup)))
                cyclic_index_two_cases += 1

    if failures:
        raise AssertionError(f"counterexamples: {failures[:3]}")
    print(f"GROUPS_ORDER_LE_16 {len(groups)}")
    print(f"INVOLUTORY_AUTOMORPHISMS {total_automorphisms}")
    print("CHECKED_INSTANCES_BY_ORDER " + " ".join(
        f"{order}:{order_counts.get(order, 0)}" for order in sorted(order_counts)
    ))
    print(f"CHECKED_W_SUBGROUP_ALPHA_INVARIANT_PROPER_H {checked}")
    print(f"SKIPPED_W_NOT_SUBGROUP_AUTOMORPHISMS {skipped_non_subgroup_w}")
    print(f"TRIVIAL_H_CASES {trivial_h_cases}")
    print(f"ALPHA_FIXED_RIGHT_COSET_CASES {fixed_coset_cases}")
    print(f"CYCLIC_INDEX_TWO_FIDELITY_CASES {cyclic_index_two_cases}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
