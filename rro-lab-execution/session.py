"""Unchanged physical source provider and ML-R/Q real-receipt execution."""
from __future__ import annotations

from collections import Counter
from lab import Arena, Buffer, Lab
from framing import Receive, Print, Reader, REPLIES, reply, v

TO_BITS = str.maketrans('LR', '01')
TO_WORD = str.maketrans('01', 'LR')


def prepare_source(bits, *, storage=None):
    """UJB5.3's actual P/S preparer, including all inherited routing overhead."""
    lab = Lab(Arena(storage=storage),capture_cuts=False)
    for j in range(len(bits)):
        p = lab.at('R'*j)
        assert lab.grow('P', p) and lab.grow('S', p)
        assert lab.grow('S', p.fields['left'])
    p = lab.at('R'*len(bits))
    assert lab.grow('P', p) and lab.grow('S', p)
    assert lab.grow('S', p.fields['left'])
    assert lab.grow('S', p.fields['left'].fields['left'])
    for j, bit in enumerate(bits):
        q = lab.at('R'*j+'L')
        if bit == 0:
            assert lab.grow('P', q.fields['right'])
            assert lab.grow('S', q.fields['right'])
        else:
            assert lab.grow('S', q.fields['left'])
            assert lab.grow('P', q.fields['right'])
    assert len(lab.events) == 5*len(bits)+4
    assert lab.persistent()['V'] == 8*len(bits)+8
    lab.validate_routes()
    return lab


class OriginalPort:
    """Only call() crosses the observer boundary. Private fields stay private."""
    def __init__(self, source):
        self.source = source
        self.root = source.root
        self.grant = True
        self.calls = []
        self.observations = Counter()

    def call(self, operand, observer_arena):
        assert self.grant
        wire = Buffer(observer_arena, 'provider_request')
        wire.extend(v(operand.translate(TO_BITS)))
        observer_arena.charges['port_wire_bits'] += len(wire.cells)
        observer_arena.charges['port_logical_operand_bits'] += len(operand)
        observer_arena.charges['original_calls'] += 1
        observer_arena.charges['call_grant_checks'] += 1
        r = Reader(wire.bits())
        real_operand = r.v().translate(TO_WORD)
        r.end()
        assert real_operand == operand
        lab, p, pos, e, edges = self.source, self.root, 0, 0, 0
        while True:
            left, right = lab.a.read(p, 'left'), lab.a.read(p, 'right')
            lnull, rnull = lab.a.test('null', left is None), lab.a.test('null', right is None)
            e += 1
            tag = 'alpha' if lnull else ('beta' if rnull else 'branch')
            if pos == len(real_operand):
                result = tag
                break
            if tag != 'branch':
                result = 'absent'
                break
            p = left if real_operand[pos] == 'L' else right
            pos += 1
            edges += 1
            lab.a.charges['navigation_links'] += 1
            lab.a.charges['reference_register_writes'] += 1
            lab.a.charges['provider_operand_bit_scans'] += 1
        assert e <= len(operand)+1 and edges <= len(operand)
        observer_arena.charges['port_reply_bits'] += 2
        observer_arena.charges['port_returns'] += 1
        self.observations['examined_vertices'] += e
        self.observations['child_field_reads'] += 2*e
        self.observations['navigation_links'] += edges
        self.calls.append((operand, result, e, edges))
        # Caller captures this occupied return latch before it can be cleared.
        latch = Buffer(observer_arena, 'return_latch')
        latch.extend(format(REPLIES.index(result), '02b'))
        observer_arena.charges['port_ends'] += 1
        observer_arena.charges['port_acknowledgements'] += 1
        observer_arena.charges['port_root_returns'] += 1
        lab.a.charges['reference_register_writes'] += 1
        wire.clear()
        return latch


TABLE = {
    '': ('branch', 'branch'), 'L': ('beta', 'branch'),
    'R': ('branch', 'beta'), 'LL': ('absent', 'beta'),
    'LR': ('absent', 'alpha'), 'RL': ('beta', 'absent'),
    'RR': ('alpha', 'absent'),
}


def classify(n, word):
    """Public all-word table; returns designated original bit and two replies."""
    k, j = 15*n, 0
    while j < len(word) and j < k and word[j] == 'R':
        j += 1
    if j == k:
        # Fixed marker ((C,alpha),alpha), with C=(beta,alpha).
        marker = {'': 'branch', 'L': 'branch', 'R': 'alpha', 'LL': 'branch',
                  'LR': 'alpha', 'LLL': 'beta', 'LLR': 'alpha'}
        value = marker.get(word[j:], 'absent')
        return None, (value, value)
    if j == len(word):
        return None, ('branch', 'branch')
    assert word[j] == 'L'
    pair = TABLE.get(word[j+1:], ('absent', 'absent'))
    return (j if pair[0] != pair[1] else None), pair


class Session:
    def __init__(self, n, theta, port):
        assert n >= 4 and 0 <= theta <= 1
        self.n, self.theta, self.port = n, theta, port
        self.a = Arena(registers=24)
        # Explicit physical connection object, retained once in the chain binding.
        self.binding = self.a.alloc('original_connection')
        self.chain = self.a.alloc('receipt_chain')
        for k, value in (('binding', self.binding), ('first', None), ('last', None)):
            self.a.write(self.chain, k, value)
        self.cache = Buffer(self.a, 'cache_header')
        self.cache.extend('0'*(30*n))
        self.records = []
        self.pending = None
        self.phase = 'Ready'
        self.cuts = []

    def cache_tuple(self):
        bits = self.cache.bits()
        return tuple(tuple((h, int(bits[15*self.n+15*j+h])) for h in range(15)
                           if bits[15*j+h] == '1') for j in range(self.n))

    def cut(self):
        self.cuts.append({'phase': self.phase, 'receipts': len(self.records),
                          'pending': self.pending is not None, **self.a.census()})

    def start(self, word):
        assert self.phase == 'Ready' and set(word) <= {'L', 'R'}
        self.phase = 'Receive'
        self.receiver = Receive(self.a)
        wire = '01'+v(word.translate(TO_BITS))
        self.receiver.feed(wire)
        tag, payload = self.receiver.seal()
        assert tag == '01' and payload.translate(TO_WORD) == word
        self.pending = {'handle': self.port, 'operand': self.receiver.buffer}
        self.phase = 'WaitReply'
        self.cut()

    def capture(self):
        assert self.phase == 'WaitReply'
        wire = self.receiver.buffer.bits()
        operand = Reader(wire[2:]).v().translate(TO_WORD)
        latch = self.port.call(operand, self.a)
        captured = Buffer(self.a, 'captured_reply')
        captured.extend(latch.bits())
        self.pending.update(reply=captured, word=operand)
        latch.clear()
        self.phase = 'Captured'
        self.cut()

    def append(self):
        assert self.phase == 'Captured'
        word = self.pending['word']
        payload = v(word.translate(TO_BITS))+self.pending['reply'].bits()
        record = self.a.alloc('receipt')
        for key in ('head', 'tail', 'next'):
            self.a.write(record, key, None)
        backing = Buffer.__new__(Buffer)
        backing.a, backing.header, backing.cells = self.a, record, []
        backing.extend(payload)
        self.pending.update(record=record, payload=backing)
        self.phase = 'ReceiptUnlinked'
        self.cut()
        tail = self.a.read(self.chain, 'last')
        if self.a.test('null', tail is None):
            self.a.write(self.chain, 'first', record)
        else:
            self.a.write(tail, 'next', record)
        self.phase = 'ReceiptLinked'
        self.cut()
        self.a.write(self.chain, 'last', record)
        self.records.append((word, REPLIES[int(self.pending['reply'].bits(), 2)], backing))
        self.phase = 'ReceiptTailAdvanced'
        self.cut()

    def update(self):
        assert self.phase == 'ReceiptTailAdvanced'
        word, value, _ = self.records[-1]
        b, pair = classify(self.n, word)
        if b is not None:
            bit = pair.index(value)
            # ML-T specifies paid linear record selection, including repeated tests.
            p = self.a.read(self.cache.header, 'head')
            for _ in range(b):
                p = self.a.read(p, 'next', True)
                self.a.charges['cache_selection_steps'] += 1
            mask = self.a.read(p, 'bit')
            q = p
            for _ in range(15*self.n):
                q = self.a.read(q, 'next', True)
                self.a.charges['cache_selection_steps'] += 1
            old = self.a.read(q, 'bit')
            if self.a.test('bit', mask == 1):
                assert old == bit
                self.a.charges['consistent_rereads'] += 1
            self.a.write(p, 'bit', 1)
            self.a.write(q, 'bit', bit)
        else:
            assert value == pair[0]
        self.phase = 'BuildReply'

    def finish(self):
        assert self.phase == 'BuildReply'
        flow = Print(self.a, reply(2, self.pending['reply'].bits()))
        self.phase = 'Print'
        self.cut()
        output = flow.finish()
        self.receiver.buffer.clear()
        self.pending['reply'].clear()
        self.pending = None
        self.phase = 'Ready'
        self.cut()
        self.check_receipts()
        return output

    def read(self, word):
        self.start(word)
        self.capture()
        self.append()
        self.update()
        return self.finish()

    def check_receipts(self):
        p = self.a.read(self.chain, 'first')
        found = []
        while p is not None:
            bits, q = '', self.a.read(p, 'head')
            while q is not None:
                bits += str(self.a.read(q, 'bit'))
                q = self.a.read(q, 'next', True)
            r = Reader(bits)
            found.append((r.v().translate(TO_WORD), REPLIES[int(r.take(2), 2)]))
            r.end()
            p = self.a.read(p, 'next', True)
        assert found == [(w, val) for w, val, _ in self.records]
        assert found == [(w, val) for w, val, _, _ in self.port.calls]
        payload_size = sum(2*len(w)+3 for w, _, _ in self.records)
        count = len(self.records)
        actual = [self.chain]+[b.header for _, _, b in self.records]+[
            cell for _, _, b in self.records for cell in b.cells]
        assert len(actual) == 1+count+payload_size
        assert sum(sum(k != 'bit' for k in o.fields) for o in actual) == 3+3*count+payload_size
        assert sum(o.writes for o in actual) == 3+5*count+4*payload_size
        assert len(self.cache.cells) == 30*self.n
        return {'receipts': count, 'payload_bits': payload_size,
                'objects': len(actual), 'reference_slots': 3+3*count+payload_size,
                'initialized_and_append_fields':3+5*count+4*payload_size}
