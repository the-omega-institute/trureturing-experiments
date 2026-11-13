"""ML-C/L/B with AA20 physical unary markers and literal ancestor tapes.

Object references are Python object identities, never an integer-to-object map.
The primitive ledger models the conditional ML-I installation, not CPython cost.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field


@dataclass(eq=False)
class Obj:
    kind: str
    fields: dict = field(default_factory=dict)
    writes: int = 0  # Audit metadata, outside the simulated object's fields.


class Arena:
    """Widen actual reference records; retain all allocated material and grants."""
    def __init__(self, registers=16, storage=None):
        self.storage=storage
        self.storage_owner=storage.owner() if storage is not None else None
        self.objects = []
        self.identities = set()
        self.reference_fields = 0
        self.bit_fields = 0
        self.free_pool = None
        self.charges = Counter()
        self.width = 1
        self.registers = registers
        self.peak_bits = 0
        self.peak_objects = 0
        self.peak_reference_slots = 0

    def slots(self):
        return self.registers + self.reference_fields

    def census(self):
        refs = self.slots()
        bits = self.bit_fields
        # One retained current-grant bit per binding; code is reported separately.
        return {'objects': len(self.objects), 'reference_slots': refs, 'bit_fields': bits,
                'grant_bits': len(self.objects), 'reference_width': self.width,
                'bits': refs*self.width + bits + len(self.objects)}

    def sample(self, overlap=0):
        c = self.census()
        self.peak_bits = max(self.peak_bits, c['bits']+overlap)
        self.peak_objects = max(self.peak_objects, c['objects'])
        self.peak_reference_slots = max(self.peak_reference_slots, c['reference_slots'])

    def alloc(self, kind):
        needed = (len(self.objects)+2-1).bit_length()
        if needed > self.width:
            n = self.slots()
            self.charges['reference_widenings'] += 1
            self.charges['reference_copy_bits'] += n*(self.width+needed)
            self.charges['reference_rewrite_fields'] += n
            self.sample(overlap=n*needed)
            self.width = needed
            if self.storage is not None:self.storage.state.width[self.storage_owner]=needed
        o = self.storage.obj(self.storage_owner,kind) if self.storage is not None else Obj(kind)
        self.objects.append(o)
        self.identities.add(o)
        self.charges['allocations'] += 1
        self.charges['grant_deliveries'] += 1
        self.charges['reference_transport_bits'] += self.width
        self.sample()
        return o

    def write(self, o, key, value, located=False):
        assert o in self.identities
        if key not in o.fields:
            if key == 'bit':
                self.bit_fields += 1
            else:
                self.reference_fields += 1
        o.fields[key] = value
        o.writes += 1
        self.charges['field_writes'] += 1
        if key == 'bit':
            assert value in (0, 1, None)
            self.charges['bit_writes'] += 1
        else:
            assert value is None or isinstance(value, Obj)
            self.charges['reference_write_bits'] += self.width
        if located:
            self.charges['located_fields'] += 1
        self.sample()

    def read(self, o, key, navigation=False):
        assert isinstance(o, Obj) and key in o.fields
        self.charges['field_reads'] += 1
        value = o.fields[key]
        self.charges['bit_reads' if key == 'bit' else 'reference_read_bits'] += 1 if key == 'bit' else self.width
        if navigation:
            self.charges['navigation_links'] += 1
        return value

    def test(self, kind, value):
        self.charges[kind+'_tests'] += 1
        return bool(value)

    def node(self, parent=None):
        n = self.alloc('source')
        for k, v in (('parent', parent), ('left', None), ('right', None), ('header', None)):
            self.write(n, k, v, True)
        return n

    def release_bit_cell(self, cell):
        if self.free_pool is None:
            self.free_pool = self.alloc('free_pool_header')
            self.write(self.free_pool,'head',None)
            self.write(self.free_pool,'tail',None)
        old = self.read(self.free_pool,'head')
        self.write(cell,'next',old)
        self.write(self.free_pool,'head',cell)
        self.charges['free_list_writes'] += 2
        if self.test('null',old is None):
            self.write(self.free_pool,'tail',cell)
            self.charges['free_list_writes'] += 1


class Buffer:
    def __init__(self, arena, kind='buffer'):
        self.a = arena
        self.header = arena.alloc(kind)
        arena.write(self.header, 'head', None)
        arena.write(self.header, 'tail', None)
        self.cells = []

    def append(self, bit):
        a = self.a
        cell = a.alloc('bitcell')
        a.write(cell, 'bit', bit)
        a.write(cell, 'next', None)
        tail = a.read(self.header, 'tail')
        if a.test('null', tail is None):
            a.write(self.header, 'head', cell)
        else:
            a.write(tail, 'next', cell)
        a.write(self.header, 'tail', cell)
        self.cells.append(cell)

    def extend(self, bits):
        for b in bits:
            self.append(int(b))

    def bits(self):
        a, out = self.a, []
        p = a.read(self.header, 'head')
        while not a.test('null', p is None):
            out.append(str(a.read(p, 'bit')))
            p = a.read(p, 'next', navigation=True)
        return ''.join(out)

    def clear(self):
        # Retain and physically link reusable cells in the shared paid pool.
        cell = self.a.read(self.header,'head')
        while not self.a.test('null',cell is None):
            nxt = self.a.read(cell,'next',True)
            self.a.write(cell, 'bit', None)
            self.a.write(cell, 'next', None)
            self.a.charges['cleanup_fields'] += 2
            self.a.release_bit_cell(cell)
            cell = nxt
        self.a.write(self.header, 'head', None)
        self.a.write(self.header, 'tail', None)
        self.a.charges['cleanup_fields'] += 2


class Lab:
    def __init__(self, arena=None, capture_cuts=True):
        self.a = arena or Arena()
        self.root = self.a.node()
        # Ownership is audit metadata; actual references remain object identities.
        self.owned = {self.root}
        self.events = []
        self.event_counts = Counter()
        self.maintenance_visits = 0
        self.pending = None
        self.cuts = []
        self.capture_cuts = capture_cuts
        self.observer = None

    def tag(self, n):
        a = self.a
        left, right = a.read(n, 'left'), a.read(n, 'right')
        lnull, rnull = a.test('null', left is None), a.test('null', right is None)
        a.charges['arity_tests'] += 1
        if lnull and rnull:
            return 'alpha'
        if not lnull and rnull:
            return 'beta'
        assert not lnull and not rnull
        return 'branch'

    @staticmethod
    def peek_tag(n):
        """Independent inspector: does not insert audit scans in the driver ledger."""
        left, right = n.fields['left'], n.fields['right']
        if left is None and right is None:
            return 'alpha'
        if left is not None and right is None:
            return 'beta'
        assert left is not None and right is not None
        return 'branch'

    def at(self, word):
        p = self.root
        for d in word:
            if self.tag(p) != 'branch':
                return None
            p = self.a.read(p, 'left' if d == 'L' else 'right', navigation=True)
        return p

    def address(self, target):
        word, p = [], target
        while p is not self.root:
            parent = p.fields['parent']
            assert parent is not None
            word.append('L' if parent.fields['left'] is p else 'R')
            p = parent
        return ''.join(reversed(word))

    def cut(self, phase):
        if self.observer is not None:
            self.observer(phase)
        if self.capture_cuts:
            self.cuts.append({'phase': phase, 'pending': self.pending is not None,
                              'events': len(self.events), **self.a.census()})

    def grow(self, kind, target):
        """Reserve, physical edit, then every proper-ancestor append, in order."""
        a = self.a
        if self.tag(target) != ('alpha' if kind == 'P' else 'beta'):
            a.charges['guard_refusals'] += 1
            return False
        assert self.pending is None
        addr = self.address(target)
        self.pending = {'kind': kind, 'target': target, 'cursor': target, 'fresh': []}
        self.cut('Reserve')
        if kind == 'P':
            marker = a.node(target)
            self.owned.add(marker)
            self.pending['fresh'].append(marker)
            self.cut('ReservedP')
            a.write(target, 'left', marker, True)
        else:
            old = a.read(target, 'left')
            marker, right = a.node(old), a.node(target)
            header = a.alloc('route_header')
            a.write(header, 'head', None, True)
            a.write(header, 'tail', None, True)
            self.owned.update((marker, right, header))
            self.pending['fresh'].extend((marker, right, header))
            self.cut('ReservedS')
            a.write(old, 'left', marker, True)
            a.write(target, 'right', right, True)
            a.write(target, 'header', header, True)
        self.events.append((kind, addr))
        self.event_counts[kind] += 1
        a.charges[kind] += 1
        self.cut('EditComplete')
        child = target
        while True:
            parent = a.read(child, 'parent')
            if a.test('null', parent is None):
                break
            a.charges['maintenance_ancestor_visits'] += 1
            self.maintenance_visits += 1
            left, right = a.read(parent, 'left'), a.read(parent, 'right')
            side = 0 if a.test('equality', child is left) else 1
            assert child is (left if side == 0 else right)
            header = a.read(parent, 'header')
            assert header is not None
            cell = a.alloc('route_cell')
            self.owned.add(cell)
            self.pending.update(cursor=parent, append=cell)
            a.write(cell, 'bit', side, True)
            a.write(cell, 'next', None, True)
            self.cut('AppendInitialized')
            tail = a.read(header, 'tail')
            if a.test('null', tail is None):
                a.write(header, 'head', cell, True)
            else:
                a.write(tail, 'next', cell, True)
            self.cut('AppendLinked')
            a.write(header, 'tail', cell, True)
            self.cut('AppendTailAdvanced')
            a.charges['route_appends'] += 1
            a.charges['navigation_links'] += 1
            child = parent
        self.pending = None
        self.cut('Maintained')
        return True

    def logical(self, n=None):
        n = self.root if n is None else n
        tag = self.peek_tag(n)
        return tag if tag != 'branch' else (self.logical(n.fields['left']), self.logical(n.fields['right']))

    def routing(self):
        out = {}
        for o in self.owned:
            if o.kind != 'source' or o.fields['header'] is None:
                continue
            h = o.fields['header']
            p, bits = h.fields['head'], ''
            while p is not None:
                bits += str(p.fields['bit'])
                p = p.fields['next']
            out[self.address(o)] = bits
        return out

    def validate_routes(self):
        expected = {p: '' for kind, p in self.events if kind == 'S'}
        for _, word in self.events:
            for j in range(len(word)):
                assert word[:j] in expected
                expected[word[:j]] += str(int(word[j] == 'R'))
        assert self.routing() == expected

    def persistent(self):
        objs = [o for o in self.owned if o.kind in ('source', 'route_header', 'route_cell')]
        return {'V': sum(o.kind == 'source' for o in objs),
                'H': sum(o.kind == 'route_header' for o in objs),
                'D': sum(o.kind == 'route_cell' for o in objs),
                'objects': len(objs),
                'reference_slots': sum(sum(k != 'bit' for k in o.fields) for o in objs)}


def encode_comb(minus, plus):
    def arm(bits):
        t = 'beta'
        for b in reversed(bits):
            t = (t, 'beta' if b else 'alpha')
        return t
    return arm(minus), arm(plus)


class Machine:
    def __init__(self, program, bits='', arena=None, audit=True):
        arena = arena or Arena()
        self.audit = audit
        self.lab = Lab(arena, capture_cuts=audit)
        self.program = tuple(program)
        self.pc, self.arm, self.commits = 0, '+', 0
        self.head_navigation = 0
        self.busy = None
        self.trace = []
        self.input = Buffer(self.lab.a, 'input_header')
        self.input.extend(bits)
        backed = self.input.bits()
        l = self.lab
        assert l.grow('P', l.root)
        assert l.grow('S', l.root)
        assert l.grow('P', l.root.fields['right'])
        p = l.root.fields['right']
        for j in range(max(1, len(backed))):
            assert l.grow('S', p)
            if j < len(backed) and backed[j] == '1':
                assert l.grow('P', p.fields['right'])
            p = p.fields['left']
        self.head = l.root.fields['right']
        self.micro_boundaries = []
        l.observer = self.observe_busy
        self.boundary()

    def observe_busy(self, phase):
        assert self.busy is not None
        assert self.pc == self.busy['old_pc'] and self.head is self.busy['old_head']
        self.busy['phase'] = phase
        assert {'proposed_head','new_arm','next_pc'} <= self.busy.keys()
        if not self.audit:
            return
        self.micro_boundaries.append({'phase':phase,'commits':self.commits,
                                      'old_pc':self.pc,'old_head':self.lab.address(self.head),
                                      'pending_edit':self.lab.pending is not None})

    def decode(self, charged=False):
        root = self.lab.root
        tag = self.lab.tag if charged else self.lab.peek_tag
        def link(node, key):
            return self.lab.a.read(node,key,True) if charged else node.fields[key]
        values, lengths = {}, []
        for arm in ('-', '+'):
            p = link(root,'left' if arm == '-' else 'right')
            j = 0
            while tag(p) == 'branch':
                value = tag(link(p,'right'))
                assert value in ('alpha', 'beta')
                values[-j-1 if arm == '-' else j] = int(value == 'beta')
                j += 1
                p = link(p,'left')
            assert tag(p) == 'beta'
            lengths.append(j)
        word = self.lab.address(self.head)
        assert word[0] == ('L' if self.arm == '-' else 'R') and set(word[1:]) <= {'L'}
        h = -len(word) if self.arm == '-' else len(word)-1
        return values, h, self.pc, tuple(lengths)

    def boundary(self):
        l, a = self.lab, self.lab.a
        values, h, pc, (nm, np) = self.decode()
        assert np >= 1
        assert l.logical() == encode_comb([values[-j-1] for j in range(nm)], [values[j] for j in range(np)])
        l.validate_routes()
        c, b = nm+np, sum(values.values())
        p, s = (l.event_counts[k] for k in ('P','S'))
        counts = l.persistent()
        d = 1+np*(np+1)//2+nm*(nm+1)//2+sum(2+(j if j >= 0 else -j-1) for j, v in values.items() if v)
        assert (p, s, counts['V'], counts['H'], counts['D']) == (2+b, 1+c, 5+2*c+b, 1+c, d)
        assert sum(o.writes for o in l.owned) == 3+4*p+11*s+counts['V']+4*d
        assert counts['reference_slots'] == 4*counts['V']+2*counts['H']+d
        assert l.maintenance_visits == d
        self.trace.append({'head': h, 'pc': pc, 'marked': sorted(j for j, v in values.items() if v),
                           'lengths': [nm, np], 'P': p, 'S': s, **counts})

    def step(self):
        assert self.busy is None, 'No instruction is enabled before current Commit'
        if self.pc == len(self.program):
            return False
        l, a = self.lab, self.lab.a
        ins = self.program[self.pc]
        old_head, old_arm, old_pc = self.head, self.arm, self.pc
        new_head, new_arm, new_pc = old_head, old_arm, old_pc+1
        self.busy = {'phase': 'Fetch', 'old_head': old_head, 'old_arm': old_arm,
                     'old_pc': old_pc, 'instruction': ins,
                     'proposed_head':new_head,'new_arm':new_arm,'next_pc':new_pc}
        if ins == 'M' or isinstance(ins, tuple):
            bit = a.read(self.head, 'right')
            tag = l.tag(bit)
            self.busy['phase'] = 'Inspect'
            if ins == 'M' and tag == 'alpha':
                assert l.grow('P', bit)
            elif isinstance(ins, tuple):
                assert ins[0] == 'J' and 0 <= ins[1] < len(self.program)
                new_pc = ins[1] if tag == 'beta' else old_pc+1
        else:
            assert ins in ('L', 'R')
            outward = (old_arm == '+' and ins == 'R') or (old_arm == '-' and ins == 'L')
            if outward:
                new_head = a.read(old_head, 'left', True)
                a.charges['head_navigation_links'] += 1
                self.head_navigation += 1
            else:
                parent = a.read(old_head, 'parent', True)
                a.charges['head_navigation_links'] += 1
                self.head_navigation += 1
                if a.test('equality', parent is l.root):
                    new_arm = '-' if old_arm == '+' else '+'
                    new_head = a.read(l.root, 'left' if new_arm == '-' else 'right', True)
                    a.charges['head_navigation_links'] += 1
                    self.head_navigation += 1
                else:
                    new_head = parent
            self.busy.update(phase='Plan',proposed_head=new_head,new_arm=new_arm,next_pc=new_pc)
            if l.tag(new_head) == 'beta':
                assert l.grow('S', new_head)
            else:
                assert l.tag(new_head) == 'branch'
        self.busy.update(phase='ReturnHead', proposed_head=new_head, next_pc=new_pc, new_arm=new_arm)
        assert l.pending is None
        self.head, self.arm, self.pc = new_head, new_arm, new_pc
        a.charges['reference_register_writes'] += 1
        self.commits += 1
        self.busy = None
        if self.audit:
            self.boundary()
        assert self.head_navigation <= 2*self.commits
        return True

    def run(self, bound):
        for _ in range(bound):
            if not self.step():
                break
        return self.trace[-1]

    def output_bits(self, left, right):
        assert self.pc == len(self.program) and left <= right
        values, _, _, _ = self.decode(charged=True)
        return ''.join(str(values.get(j, 0)) for j in range(left, right+1))


def reference_run(program, bits, bound):
    """Independent Z-indexed abstract B interpreter; no comb or link operations."""
    tape = {j for j, b in enumerate(bits) if b == '1'}
    h, pc, trace = 0, 0, [(0, 0, sorted(tape))]
    locations = []
    for _ in range(bound):
        if pc == len(program):
            break
        locations.append(pc)
        ins = program[pc]
        if ins == 'M':
            tape.add(h)
        elif ins == 'L':
            h -= 1
        elif ins == 'R':
            h += 1
        elif isinstance(ins, tuple):
            pc = ins[1] if h in tape else pc+1
            trace.append((h, pc, sorted(tape)))
            continue
        else:
            raise ValueError(ins)
        pc += 1
        trace.append((h, pc, sorted(tape)))
    return trace, locations
