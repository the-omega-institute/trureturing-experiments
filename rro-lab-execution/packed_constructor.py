"""Owner-tagged packed host representation for the unchanged ML-I/LAB path.

The packed IDs are serialization handles only.  Source-visible identity is the
(host-owner, slot) pair and is decoded through the retained reference graph.
This module does not supply an evaluator, source answer, or hot-refill path.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import ceil, log2

Shape = tuple[int, int]
HostId = tuple[int, int]
SHAPES: tuple[Shape, ...] = ((4, 8), (4, 21), (1, 0), (0, 1), (0, 32768), (2, 0))


@dataclass(slots=True)
class PackedSlot:
    owner: int
    shape: Shape
    schema: int
    refs: list[HostId | None]
    bits: bytearray
    grant: int = 1

    def put_bit(self, index: int, value: int) -> None:
        assert value in (0, 1) and 0 <= index < self.shape[1]
        byte, offset = divmod(index, 8)
        mask = 1 << offset
        self.bits[byte] = (self.bits[byte] & ~mask) | (value << offset)

    def get_bit(self, index: int) -> int:
        assert 0 <= index < self.shape[1]
        byte, offset = divmod(index, 8)
        return (self.bits[byte] >> offset) & 1


@dataclass(slots=True)
class ControllerHeader:
    pc: HostId | None = None
    old_pc: HostId | None = None
    rx: HostId | None = None
    ry: HostId | None = None
    destination: HostId | None = None
    pending: HostId | None = None
    cursor: HostId | None = None
    phase: int = 0
    selector: int = 0
    waiting: int = 0
    ended: int = 0
    acknowledged: int = 0
    occupied: int = 0


@dataclass(slots=True)
class ReceiptRow:
    occurrence: int
    original_index: int
    operand: tuple[int, ...]
    returned: HostId | None
    mask: HostId | None
    value: HostId | None


@dataclass(slots=True)
class SessionPacking:
    owners: dict[int, list[PackedSlot]] = field(default_factory=dict)
    receipts: list[ReceiptRow] = field(default_factory=list)
    cache_positions: int = 0


class PackedState:
    """Finite owner-separated arenas with literal field operations."""

    def __init__(self, owner_count: int = 4) -> None:
        assert owner_count > 0
        self.owners: dict[int, list[PackedSlot]] = {i: [] for i in range(owner_count)}
        self.width: dict[int, int] = {i: 1 for i in range(owner_count)}
        self.free: dict[int, list[HostId]] = {i: [] for i in range(owner_count)}
        self.header = ControllerHeader()
        self.events: list[tuple] = []

    def allocate(self, owner: int, shape: Shape, schema: int = 0) -> HostId:
        assert owner in self.owners and shape in SHAPES
        slot = PackedSlot(owner, shape, schema, [None] * shape[0], bytearray(ceil(shape[1] / 8)))
        self.owners[owner].append(slot)
        handle = (owner, len(self.owners[owner]) - 1)
        self.events.append(("allocate", handle, shape, schema))
        return handle

    def slot(self, handle: HostId) -> PackedSlot:
        owner, index = handle
        assert owner in self.owners and 0 <= index < len(self.owners[owner])
        return self.owners[owner][index]

    def grant(self, handle: HostId, value: int) -> None:
        assert value in (0, 1)
        self.slot(handle).grant = value
        self.events.append(("grant", handle, value))

    def read_ref(self, handle: HostId, field: int) -> HostId | None:
        slot = self.slot(handle)
        assert 0 <= field < slot.shape[0]
        value = slot.refs[field]
        self.events.append(("read-ref", handle, field, value))
        return value

    def write_ref(self, handle: HostId, field: int, value: HostId | None) -> None:
        slot = self.slot(handle)
        assert 0 <= field < slot.shape[0]
        if value is not None:
            self.slot(value)
        slot.refs[field] = value
        self.events.append(("write-ref", handle, field, value))

    def read_bit(self, handle: HostId, field: int) -> int:
        value = self.slot(handle).get_bit(field)
        self.events.append(("read-bit", handle, field, value))
        return value

    def write_bit(self, handle: HostId, field: int, value: int) -> None:
        self.slot(handle).put_bit(field, value)
        self.events.append(("write-bit", handle, field, value))

    def shape(self, handle: HostId) -> Shape:
        return self.slot(handle).shape

    def null(self, handle: HostId | None) -> bool:
        self.events.append(("null", handle))
        return handle is None

    def release_to_free(self, handle: HostId) -> None:
        self.slot(handle).grant = 1
        self.free[handle[0]].append(handle)
        self.events.append(("cleanup", handle))

    def snapshot(self) -> dict[HostId, tuple]:
        out: dict[HostId, tuple] = {}
        for owner, slots in self.owners.items():
            for index, slot in enumerate(slots):
                out[(owner, index)] = (slot.owner, slot.shape, slot.schema,
                                       tuple(slot.refs), bytes(slot.bits), slot.grant)
        return out

    def modeled_census(self) -> dict[int, dict[str, int]]:
        result: dict[int, dict[str, int]] = {}
        for owner, slots in self.owners.items():
            refs = sum(slot.shape[0] for slot in slots)
            bits = sum(slot.shape[1] for slot in slots)
            result[owner] = {
                "objects": len(slots), "reference_slots": refs, "bit_fields": bits,
                "grant_bits": sum(slot.grant for slot in slots), "reference_width": self.width[owner],
                "bits": refs * self.width[owner] + bits + len(slots),
            }
        return result

    def packed_bits(self) -> int:
        count = sum(len(slots) for slots in self.owners.values())
        index_bits = max(1, ceil(log2(count + 1)))
        owner_bits = max(1, ceil(log2(len(self.owners) + 1)))
        total = 0
        for owner, slots in self.owners.items():
            for slot in slots:
                refs, bits = slot.shape
                total += 1 + 4 + 4 + 1 + slot.schema.bit_length()
                total += refs * (owner_bits + index_bits) + bits
        return total

    def relocate(self, source: list[HostId], target: list[HostId]) -> None:
        """Copy fields through the actual forward map, preserving aliases."""
        assert len(source) == len(target)
        forward = dict(zip(source, target))
        for old, new in zip(source, target):
            old_slot, new_slot = self.slot(old), self.slot(new)
            assert old_slot.shape == new_slot.shape
            for i in range(old_slot.shape[0]):
                value = self.read_ref(old, i)
                self.write_ref(new, i, None if value is None else forward[value])
            for i in range(old_slot.shape[1]):
                self.write_bit(new, i, self.read_bit(old, i))

    def decode_aliases(self) -> dict[HostId, tuple[HostId | None, ...]]:
        """Decoded graph view; equal handles remain equal and null stays null."""
        return {(owner, i): tuple(slot.refs)
                for owner, slots in self.owners.items()
                for i, slot in enumerate(slots)}
