"""Explicit append-copy binary TM and its primary NWMG thirteen-block lowering.

Eight-bit blocks retain the paper's head/value triple plus monotone cursor,
boundary and handoff metadata. No transition may erase a one. The generator
accepts an arbitrary finite *erasing* binary TM, without arithmetic callbacks.
"""
from __future__ import annotations

WIDTH = 8
HALT = ':halt'
# Type bits 0/1/2 distinguish header/data/end from unused zero blocks.
# Header: active=3, retired=4, new left-guard head=5.
# Data: old head=3, value zero/one=4/5, consumed=6, new head=7.
# End: new right-guard head=3.


def initial_word(word):
    assert set(word) <= {'0', '1'}
    blocks = ['10010000']
    for j, value in enumerate(word or '0'):
        b = [0]*WIDTH
        b[1], b[3], b[4+int(value)] = 1, int(j == 0), 1
        blocks.append(''.join(map(str, b)))
    blocks.append('00100000')
    return ''.join(blocks)


class Assembly:
    """A finite block graph expanded to ordinary non-erasing binary rules."""
    def __init__(self):
        self.ops = {}

    def add(self, label, *op):
        assert label not in self.ops
        self.ops[label] = op

    def compile(self, start):
        rules = {}
        def move(label, direction, destination):
            assert label not in rules
            rules[label] = ((0, direction, destination), (1, direction, destination))
        for label, op in self.ops.items():
            kind = op[0]
            if kind == 'move':
                _, direction, destination = op
                for j in range(WIDTH):
                    move(label if j == 0 else (label, j), direction,
                         destination if j == WIDTH-1 else (label, j+1))
            elif kind in ('test', 'set'):
                bit = op[1]
                assert 0 <= bit < WIDTH
                for j in range(bit):
                    move(label if j == 0 else (label, 'go', j), 'R',
                         (label, 'read') if j+1 == bit else (label, 'go', j+1))
                read = label if bit == 0 else (label, 'read')
                if kind == 'test':
                    rules[read] = tuple((b, 'R', (label, 'back', b, bit+1)) for b in (0, 1))
                    for b in (0, 1):
                        for j in range(bit+1, 0, -1):
                            move((label, 'back', b, j), 'L',
                                 op[2+b] if j == 1 else (label, 'back', b, j-1))
                else:
                    rules[read] = ((1, 'R', (label, 'back', bit+1)),)*2
                    for j in range(bit+1, 0, -1):
                        move((label, 'back', j), 'L', op[2] if j == 1 else (label, 'back', j-1))
            else:
                raise AssertionError(op)
        order = [start]+[q for q in rules if q != start]
        positions = {q: i for i, q in enumerate(order)}
        positions[HALT] = len(order)
        return tuple(tuple((w, d, positions[q]) for w, d, q in rules[k]) for k in order), positions


class Copier:
    def __init__(self, transitions):
        self.source = tuple(tuple(row) for row in transitions)
        if not self.source:
            self.rows, self.boundaries = (), {0: 0}
            self.description = {'block_width': WIDTH, 'source_states': 0,
                                'block_control_states': 0, 'binary_states': 0}
            return
        for row in self.source:
            assert len(row) == 2
            for w, d, q in row:
                assert w in (0, 1) and d in ('L', 'R') and 0 <= q <= len(self.source)
        a = Assembly()
        def entry(q):
            return HALT if q == len(self.source) else ('entry', q)
        for qi, row in enumerate(self.source):
            find = ('head', qi)
            a.add(entry(qi), 'move', 'R', find)
            a.add(find, 'test', 3, ('head-next', qi), ('value', qi))
            a.add(('head-next', qi), 'move', 'R', find)
            a.add(('value', qi), 'test', 5, (qi, 0, 'neighbor'), (qi, 1, 'neighbor'))
            for read, (write, direction, dest) in enumerate(row):
                def lab(s):
                    return (qi, read, s)
                def add(s, *op):
                    a.add(lab(s), *op)
                def test(s, bit, z, o):
                    add(s, 'test', bit, lab(z), lab(o))
                def move(s, d, nxt):
                    add(s, 'move', d, lab(nxt))
                def sets(s, bits, nxt):
                    assert bits
                    for j, bit in enumerate(bits):
                        add(s if j == 0 else (s, j), 'set', bit,
                            lab(nxt if j == len(bits)-1 else (s, j+1)))
                move('neighbor', direction, 'neighbor-kind')
                test('neighbor-kind', 1, 'guard', 'target')
                sets('target', [7], 'old-header')
                sets('guard', [5 if direction == 'L' else 3], 'old-header')
                test('old-header', 0, 'old-header-left', 'left-head')
                move('old-header-left', 'L', 'old-header')
                test('left-head', 5, ('end', 0), ('end', 1))
                for head in (0, 1):
                    test(('end', head), 2, ('end-next', head), ('header-next', head))
                    move(('end-next', head), 'R', ('end', head))
                    move(('header-next', head), 'R', ('new-header', head))
                    sets(('new-header', head), [0], ('left-next', head))
                    move(('left-next', head), 'R', ('left-data', head))
                    sets(('left-data', head), [1, 4]+([3] if head else []), ('left-finish', head))
                    move(('left-finish', head), 'R', 'seek-old')
                test('seek-old', 0, 'seek-left', 'seek-active')
                test('seek-active', 3, 'seek-left', 'seek-retired')
                test('seek-retired', 4, 'scan-next', 'seek-left')
                move('seek-left', 'L', 'seek-old')
                move('scan-next', 'R', 'scan')
                test('scan', 2, 'consumed', 'right-head')
                test('consumed', 6, 'override', 'scan-next')
                test('override', 3, 'copy-value', ('copy-head', write))
                test('copy-value', 5, ('copy-head', 0), ('copy-head', 1))
                for value in (0, 1):
                    test(('copy-head', value), 7, ('consume', value, 0), ('consume', value, 1))
                    for head in (0, 1):
                        tag = ('copy', value, head)
                        sets(('consume', value, head), [6], ('frontier', tag))
                test('right-head', 3, ('frontier', ('right', 0, 0)), ('frontier', ('right', 0, 1)))
                tags = [('copy', v, h) for v in (0, 1) for h in (0, 1)]
                tags += [('right', 0, h) for h in (0, 1)]
                for tag in tags:
                    _, value, head = tag
                    test(('frontier', tag), 0, ('frontier1', tag), ('frontier-next', tag))
                    test(('frontier1', tag), 1, ('frontier2', tag), ('frontier-next', tag))
                    test(('frontier2', tag), 2, ('append', tag), ('frontier-next', tag))
                    move(('frontier-next', tag), 'R', ('frontier', tag))
                    sets(('append', tag), [1, 4+value]+([3] if head else []), ('append-next', tag))
                    move(('append-next', tag), 'R', 'seek-old' if tag[0] == 'copy' else 'end-new')
                sets('end-new', [2], 'activate-left')
                move('activate-left', 'L', 'activate')
                test('activate', 0, 'activate-left', 'active-set')
                sets('active-set', [3], 'retire-left')
                move('retire-left', 'L', 'retire')
                test('retire', 0, 'retire-left', 'retired-set')
                sets('retired-set', [4], 'return-next')
                move('return-next', 'R', 'return')
                add('return', 'test', 0, lab('return-next'), entry(dest))
        self.rows, positions = a.compile(entry(0))
        self.boundaries = {positions[entry(q)]: q for q in range(len(self.source)+1)}
        self.description = {'block_width': WIDTH, 'source_states': len(self.source),
                            'block_control_states': len(a.ops), 'binary_states': len(self.rows)}


def decode_records(read, header):
    """Decode from actual tape reads at the returned active-header coordinate.

    Traversing old headers counts copies and determines the original coordinate
    offset. A caller must supply storage/charges for its read/cursor realization;
    this pure decoder alone is not a full session resource certificate.
    """
    assert header >= 0 and header % WIDTH == 0
    def block(p):
        return tuple(read(p+j) for j in range(WIDTH))
    copies, pos = 0, 0
    while pos < header:
        b = block(pos)
        assert sum(b[:3]) == 1
        if b[0]:
            assert b[3:5] == (1, 1)
            copies += 1
        pos += WIDTH
    b = block(header)
    assert b[:5] == (1, 0, 0, 1, 0)
    pos, coord, marked, heads = header+WIDTH, -copies, [], []
    while True:
        b = block(pos)
        if b[2]:
            assert b[:3] == (0, 0, 1)
            break
        assert b[:3] == (0, 1, 0) and b[4]+b[5] == 1
        assert b[6:] == (0, 0), 'Active source data have not yet been consumed.'
        if b[3]:
            heads.append(coord)
        if b[5]:
            marked.append(coord)
        pos += WIDTH
        coord += 1
    assert len(heads) == 1
    return {'marked': marked, 'head': heads[0], 'state': None,
            'left': -copies, 'right': coord-1, 'copies': copies,
            'header': header, 'end': pos}
