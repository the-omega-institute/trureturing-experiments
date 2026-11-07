"""Primary NWMG pp.6-9 equations (1)-(3), unchanged thirteen-instruction blocks."""
from __future__ import annotations


def compile_non_erasing(transitions):
    """States 0..k-1 are nonhalting; k is halt. Each row has 0/1 rules."""
    program = []
    k = len(transitions)
    for i, row in enumerate(transitions):
        assert len(row) == 2
        program.extend(('R', ('J', 13*i+8)))
        for read, (write, direction, dest) in enumerate(row):
            assert write in (0, 1) and write >= read
            assert direction in ('L', 'R') and 0 <= dest <= k
            if read == 0:
                if direction == 'R':
                    block = ('R', 'M', 'M', 'M', 'M') if write == 0 else ('M', 'R', 'M', 'M', 'M')
                else:
                    block = ('L', 'L', 'L', 'M', 'M') if write == 0 else ('M', 'L', 'L', 'L', 'M')
            else:
                block = ('R', 'M', 'M', 'M') if direction == 'R' else ('L', 'L', 'L', 'M')
            program.extend(block+(('J', 13*dest),))
    program.append('M')
    assert len(program) == 13*k+1
    return tuple(program)


def encode(bits):
    return ''.join('1'+str(int(bit)) for bit in bits)
