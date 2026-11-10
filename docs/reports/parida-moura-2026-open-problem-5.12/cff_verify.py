#!/usr/bin/env python3
"""Emit exact bit-mask corroboration for d=3,...,13 as JSON."""
import json
from check import check_masks, construction_masks, middle


def main():
    cases = []
    for d in range(3, 14):
        m = middle(d) + 1
        edges = construction_masks(d)
        checks = check_masks(edges, d + 2, True)
        assert len(edges) == m
        assert checks == 2 * m * (m - 1)
        assert middle(d) < m <= middle(d + 1)
        assert middle(d + 1) < 2 * m
        cases.append(dict(d=d, r=d // 2, m=m, t_e=d + 1, t=d + 2,
                          outside_checks=checks))
    total = sum(c['outside_checks'] for c in cases)
    assert total == 8203100
    print(json.dumps(dict(arxiv='2605.12634v3', cases=cases,
                          total_outside_checks=total,
                          status='finite-corroboration-only'), indent=2))


if __name__ == '__main__':
    main()
