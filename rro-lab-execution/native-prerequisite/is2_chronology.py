"""Unexecuted complete IS2 consumer; callers supply installation and backing.

This file has no command-line execution entry. It uses the fixed Direct image,
actual acquired cache, computed output and original connection throughout.
Limits are demands supplied by the certificate, never resource offers.
"""
from __future__ import annotations

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent)]
import immutable_spine as im
from is2_service_certificate import ServiceDirect
from framing import Reader, arithmetic_input, parse_arithmetic, parse_reply, v


class ChronologyDirect(im.Direct):
    """Same transition and observer; replace the tiny witness's fixed caps."""
    ref = ServiceDirect.ref
    value = ServiceDirect.value
    observe = ServiceDirect.observe
    cache = ServiceDirect.cache
    receipts = ServiceDirect.receipts

    def __init__(self, program, capacity, kernel, connection, native_limit):
        assert program.fingerprint() == 'db4fb22177743a3aafd286c5431b42d09460ab47cfe4ca488d40dd1f28600b4a'
        self.native_limit = native_limit
        super().__init__(program, '', capacity, kernel,
                         connection=connection, service_return=True)
        self.handoffs = []
        self.return_cuts = []
        self.boundary_observer = self.observe

    def tick(self):
        assert self.steps < self.native_limit
        super().tick()

    def drive(self):
        # Direct.run deliberately joins the whole output on every invocation.
        # Its ignored return still has a transient allocation and copying cost.
        self.run(self.native_limit)
        assert self.done and not self.refused

    deliver = ServiceDirect.deliver
    feed = ServiceDirect.feed


def acquisition(n, cache):
    source_bits = [dict(row).get(h, 0) for row in cache for h in range(15)]
    words = ['', 'LL', 'LRR', 'LLR', 'LL', 'LR', 'LLL', 'LRL']
    words += ['R'*j+'LL' for j in range(15*n)
              if (j % 15, source_bits[j]) in cache[j//15]]
    words += ['LL', 'R'*15+'LL', 'R'*45+'LL', 'L'*130,
              'R'*(15*n+10), 'LL']
    return words


def complete(e, raw, retained_B):
    """Future caller-only invocation; neither constructs nor runs an evaluator now.

    e must be the installed fixed composition, attached to a freshly authorized
    unchanged test-source instance. retained_B is a real separately installed
    mutable B owner; the returned owner tuple keeps its material alive.
    No expected numerical result or reference evaluator is consulted.
    """
    n, theta, expected_cache, cylinders = parse_arithmetic(raw)
    port = e.connection.port
    root = port.root
    assert root is port.source.root and not port.calls
    assert retained_B is not None
    e.drive()
    rd = Reader(raw)
    rd.u()
    rd.rational()
    e.feed(raw[:rd.pos])
    for word in acquisition(n, expected_cache):
        e.feed('01'+v(word.translate(str.maketrans('LR', '01'))))
    paired = e.cache()
    assert len(paired) == 15*n
    attained = tuple(tuple((h, paired[15*j+h][1]) for h in range(15)
                           if paired[15*j+h][0]) for j in range(n))
    acquired = arithmetic_input(n, theta, attained, cylinders)
    assert acquired == raw
    output_start = len(e.output)
    e.feed('00'+v(acquired))
    computed = ''.join(e.output[output_start:])
    kind, payload = parse_reply(computed)
    assert kind == 1
    assert len(e.handoffs) == len(e.return_cuts) == 1
    assert e.handoffs[0]['bits'] == acquired
    assert not e.value('eval.busy')
    before, old_ids = e.receipts()
    cache_before = e.cache()
    e.feed('01'+v('00'))
    after, new_ids = e.receipts()
    assert after[:-1] == before and new_ids[:-1] == old_ids
    assert new_ids[-1] not in old_ids
    assert after[-1] == ('LL', port.calls[-1][1])
    assert e.cache() == cache_before
    e.feed('11')
    assert e.value('closed') and port.root is root
    assert [(w, answer) for w, answer, _, _ in port.calls] == after
    return dict(owners=(e, retained_B, port.source, port, e.connection),
                computed_payload=payload, receipts=after, receipt_ids=new_ids,
                same_original_root=True, fresh_post_evaluation_receipt=True,
                closed=True, full_acceptance=False)
