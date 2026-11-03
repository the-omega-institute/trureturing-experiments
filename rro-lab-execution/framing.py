"""Bit-exact ML.4 formats and paid Receive/Print/End/Ack/Return phases."""
from __future__ import annotations

from fractions import Fraction
from lab import Buffer

REPLIES = ('alpha', 'beta', 'branch', 'absent')


class NeedMore(Exception):
    pass


class BadFrame(Exception):
    pass


def v(word):
    assert set(word) <= {'0', '1'}
    return '1'*len(word)+'0'+word


def u(n):
    assert n >= 0
    return v(bin(n)[2:])


def words(items):
    return u(len(items))+''.join(v(x) for x in items)


def rational(x):
    x = Fraction(x)
    return str(int(x < 0))+u(abs(x.numerator))+u(x.denominator)


class Reader:
    def __init__(self, bits):
        if set(bits)-{'0', '1'}:
            raise BadFrame('nonbinary input')
        self.bits, self.pos = bits, 0

    def take(self, n):
        if len(self.bits)-self.pos < n:
            raise NeedMore
        s = self.bits[self.pos:self.pos+n]
        self.pos += n
        return s

    def v(self):
        size = 0
        while self.take(1) == '1':
            size += 1
        return self.take(size)

    def u(self):
        word = self.v()
        if not word or (len(word) > 1 and word[0] == '0'):
            raise BadFrame('noncanonical unsigned integer')
        return int(word, 2)

    def rational(self):
        sign, n, d = self.take(1), self.u(), self.u()
        if not d or (not n and sign != '0'):
            raise BadFrame('invalid rational')
        return Fraction(-n if sign == '1' else n, d)

    def words(self):
        return [self.v() for _ in range(self.u())]

    def end(self):
        if self.pos != len(self.bits):
            raise BadFrame('trailing bits')


def mask_word(cache):
    entries = [dict(row) for row in cache]
    return ''.join(str(int(h in row)) for row in entries for h in range(15))+''.join(
        str(row[h]) for row in entries for h in range(15) if h in row)


def parse_mask(reader, n):
    mask = reader.take(15*n)
    cache = []
    for j in range(n):
        row = []
        for h in range(15):
            if mask[15*j+h] == '1':
                row.append((h, int(reader.take(1))))
        cache.append(tuple(row))
    return tuple(cache)


def cylinder_word(cache):
    # ML-F leaves the conflict-flag position unspecified; use its leading bit.
    return '1' if cache is None else '0'+mask_word(cache)


def arithmetic_input(n, theta, cache, cylinders):
    assert len(cache) == n
    return u(n)+rational(theta)+mask_word(cache)+words([cylinder_word(c) for c in cylinders])


def parse_arithmetic(bits):
    r = Reader(bits)
    n, theta = r.u(), r.rational()
    if n < 4 or not 0 <= theta <= 1:
        raise BadFrame('N/prior outside the declared domain')
    cache = parse_mask(r, n)
    cylinders = []
    for item in r.words():
        s = Reader(item)
        flag = s.take(1)
        cylinder = None if flag == '1' else parse_mask(s, n)
        s.end()
        cylinders.append(cylinder)
    r.end()
    return n, theta, cache, cylinders


def arithmetic_output(values, ratios):
    return words([rational(x) for x in values])+('0' if ratios is None else '1'+words([rational(x) for x in ratios]))


def parse_arithmetic_output(bits):
    r = Reader(bits)
    def parse_q(word):
        q = Reader(word)
        x = q.rational()
        q.end()
        return x
    values = [parse_q(x) for x in r.words()]
    ratios = None if r.take(1) == '0' else [parse_q(x) for x in r.words()]
    r.end()
    return values, ratios


def reply(kind, payload=''):
    assert 0 <= kind <= 3
    return format(kind, '02b')+v(payload)


def parse_reply(bits):
    r = Reader(bits)
    kind, payload = int(r.take(2), 2), r.v()
    r.end()
    if (kind in (0, 3) and payload) or (kind == 2 and len(payload) != 2):
        raise BadFrame('wrong service payload')
    return kind, payload


def parse_command(bits):
    r = Reader(bits)
    tag = r.take(2)
    if tag in ('00', '01'):
        result = (tag, r.v())
    elif tag == '10':
        variant = r.take(1)
        result = (tag, variant, r.v() if variant == '0' else r.words())
    else:
        result = (tag,)
    r.end()
    return result


class Receive:
    def __init__(self, arena):
        self.a = arena
        self.buffer = Buffer(arena, 'request_header')
        self.latch = Buffer(arena, 'input_latch')
        self.sealed = False

    def feed(self, bits):
        assert not self.sealed
        for bit in bits:
            self.a.charges['input_bits'] += 1
            self.latch.append(int(bit))
            self.buffer.append(int(self.latch.bits()))
            self.latch.clear()

    def seal(self):
        word = self.buffer.bits()
        try:
            result = parse_command(word)
        except NeedMore:
            return None
        self.sealed = True
        return result


class Print:
    def __init__(self, arena, bits):
        self.a = arena
        self.buffer = Buffer(arena, 'answer_header')
        self.buffer.extend(bits)
        self.cursor = arena.read(self.buffer.header, 'head')
        self.emitted = ''
        self.phase = 'Print'

    def bit(self):
        assert self.phase == 'Print'
        if self.a.test('null', self.cursor is None):
            self.phase = 'End'
            return False
        self.emitted += str(self.a.read(self.cursor, 'bit'))
        self.cursor = self.a.read(self.cursor, 'next', True)
        self.a.charges['output_bits'] += 1
        return True

    def end(self):
        assert self.phase == 'End'
        self.a.charges['output_ends'] += 1
        self.phase = 'Ack'

    def ack(self):
        assert self.phase == 'Ack'
        self.a.charges['acknowledgements'] += 1
        self.phase = 'Return'

    def finish(self):
        while self.phase == 'Print':
            self.bit()
        self.end()
        self.ack()
        self.buffer.clear()
        self.phase = 'Ready'
        return self.emitted
