"""Exact finite C1 packing vocabulary and deterministic bounded-probe inputs.

This module supplies no universal C1 theorem. Graph vertices are 0,...,n-1.
Every returned/verified packing may rechoose every set. Cross-set edges are allowed.
"""
from itertools import combinations, combinations_with_replacement
from math import isqrt
from random import Random
from functools import lru_cache

class NodeLimit(Exception):
    """Inconclusive computation, never a nonpacking result."""

def adjacency(n, edges):
    a = [0] * n
    seen = set()
    for u, v in edges:
        if not (0 <= u < v < n) or (u, v) in seen:
            raise ValueError('edges must be distinct unordered nonloop pairs')
        seen.add((u, v))
        a[u] |= 1 << v
        a[v] |= 1 << u
    return a

def independent(a, mask):
    while mask:
        b = mask & -mask
        mask -= b
        if a[b.bit_length() - 1] & mask:
            return False
    return True

def verify_packing(a, sizes, masks):
    if len(sizes) != len(masks):
        return False
    used = 0
    for size, mask in zip(sizes, masks):
        if mask < 0 or mask >= 1 << len(a):
            return False
        if mask.bit_count() != size or used & mask or not independent(a, mask):
            return False
        used |= mask
    return True

def exact_packing(a, sizes, max_nodes=1000000):
    """Return masks in sorted size order, or None after complete finite exclusion.

    NodeLimit means unknown. At every recursive stage all independent subsets
    of the required size are generated once, from increasing vertex indices.
    """
    sizes = tuple(sorted(sizes, reverse=True))
    nodes = 0
    @lru_cache(None)
    def rec(avail, rest):
        nonlocal nodes
        nodes += 1
        if nodes > max_nodes:
            raise NodeLimit()
        if not rest:
            return ()
        if avail.bit_count() < sum(rest):
            return None
        def subsets(candidates, need, chosen):
            if not need:
                yield chosen
                return
            while candidates.bit_count() >= need:
                b = candidates & -candidates
                candidates -= b
                yield from subsets(candidates & ~a[b.bit_length()-1], need-1, chosen | b)
        for mask in subsets(avail, rest[0], 0):
            tail = rec(avail ^ mask, rest[1:])
            if tail is not None:
                return (mask,) + tail
        return None
    return rec((1 << len(a))-1, sizes)

def parameters():
    for ell in (4,5,6):
        for kl in (2,3):
            for prior in combinations_with_replacement(range(kl,kl+3),ell-1):
                k = prior + (kl,)
                for t in range(7):
                    yield ell,k,t,sum(k)+t

def graph_inputs():
    """Yield (seed,index,ell,k,t,n,edges); the count unit is one generated input.

    Sampling may repeat a graph. It is not an enumeration of labeled graphs.
    """
    pr = list(parameters())
    for seed in (261038430,1,2,3):
        rng = Random(seed)
        for z in range(7500):
            ell,k,t,n = pr[rng.randrange(len(pr))]
            bound = (2*ell-1)*(t+1)
            pairs = list(combinations(range(n),2))
            target = min(len(pairs), max(0,bound-1-rng.randrange(max(1,bound//3))))
            if z % 4 == 0:
                edges = rng.sample(pairs,target)
            else:
                order = list(range(n)); rng.shuffle(order)
                edges = []; remaining = target; idx = 0
                while remaining and idx < n-1:
                    c = min(n-idx,max(2,(1+isqrt(1+8*remaining))//2))
                    local = list(combinations(order[idx:idx+c],2))
                    edges.extend(tuple(sorted(e)) for e in local)
                    remaining -= len(local); idx += c
                    if c < 2:
                        break
                es = set(edges)
                while len(es) < target:
                    es.add(pairs[rng.randrange(len(pairs))])
                edges = list(es)
                swaps = (z//4) % max(1,min(12,target))
                for _ in range(swaps):
                    if not edges:
                        break
                    old = edges.pop(rng.randrange(len(edges))); es.remove(old)
                    e = pairs[rng.randrange(len(pairs))]
                    while e in es:
                        e = pairs[rng.randrange(len(pairs))]
                    es.add(e); edges.append(e)
            yield seed,z,ell,k,t,n,edges
