"""Owner-tagged packed host representation for the unchanged ML-I/LAB path.

The packed IDs are serialization handles only.  Source-visible identity is the
(host-owner, slot) pair and is decoded through the retained reference graph.
This module does not supply an evaluator, source answer, or hot-refill path.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import ceil, log2, isqrt
from collections import Counter
from collections.abc import MutableSequence, MutableMapping

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

    def __init__(self, owner_count: int = 4, capture: bool = True) -> None:
        assert owner_count > 0
        self.owners: dict[int, list[PackedSlot]] = {i: [] for i in range(owner_count)}
        self.width: dict[int, int] = {i: 1 for i in range(owner_count)}
        self.free: dict[int, list[HostId]] = {i: [] for i in range(owner_count)}
        self.capture = capture
        self.header = ControllerHeader()
        self.events: list[tuple] = []

    def _record(self, event):
        if self.capture:self.events.append(event)

    def allocate(self, owner: int, shape: Shape, schema: int = 0) -> HostId:
        assert owner in self.owners and shape in SHAPES
        return self.allocate_record(owner, shape, schema)

    def allocate_record(self, owner: int, shape: Shape, schema: int = 0) -> HostId:
        assert owner in self.owners and min(shape) >= 0
        slot = PackedSlot(owner, shape, schema, PackedRefs(shape[0]), bytearray(ceil(shape[1] / 8)))
        self.owners[owner].append(slot)
        handle = (owner, len(self.owners[owner]) - 1)
        self._record(("allocate", handle, shape, schema))
        return handle

    def slot(self, handle: HostId) -> PackedSlot:
        owner, index = handle
        assert owner in self.owners and 0 <= index < len(self.owners[owner])
        return self.owners[owner][index]

    def grant(self, handle: HostId, value: int) -> None:
        assert value in (0, 1)
        self.slot(handle).grant = value
        self._record(("grant", handle, value))

    def read_ref(self, handle: HostId, field: int) -> HostId | None:
        slot = self.slot(handle)
        assert 0 <= field < slot.shape[0]
        value = slot.refs[field]
        self._record(("read-ref", handle, field, value))
        return value

    def write_ref(self, handle: HostId, field: int, value: HostId | None) -> None:
        slot = self.slot(handle)
        assert 0 <= field < slot.shape[0]
        if value is not None:
            self.slot(value)
        slot.refs[field] = value
        self._record(("write-ref", handle, field, value))

    def read_bit(self, handle: HostId, field: int) -> int:
        value = self.slot(handle).get_bit(field)
        self._record(("read-bit", handle, field, value))
        return value

    def write_bit(self, handle: HostId, field: int, value: int) -> None:
        self.slot(handle).put_bit(field, value)
        self._record(("write-bit", handle, field, value))

    def shape(self, handle: HostId) -> Shape:
        return self.slot(handle).shape

    def null(self, handle: HostId | None) -> bool:
        self._record(("null", handle))
        return handle is None

    def release_to_free(self, handle: HostId) -> None:
        self.slot(handle).grant = 1
        self.free[handle[0]].append(handle)
        self._record(("cleanup", handle))

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
                "grant_bits": len(slots), "reference_width": self.width[owner],
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


class PackedRefs(MutableSequence):
    """Null or Cantor-coded owner/slot pairs in a widening byte array.

    Widening is host serialization, independent of the charged owner width.
    No fixed machine-integer ceiling is imposed on the finite source domain.
    """
    __slots__=('data','width','count')
    def __init__(self,count=0):
        self.data=bytearray(count);self.width=1;self.count=count
    @staticmethod
    def encode(value):
        if value is None:return 0
        owner,index=value;s=owner+index
        return s*(s+1)//2+index+1
    @staticmethod
    def decode(n):
        if not n:return None
        n-=1;s=(isqrt(8*n+1)-1)//2;index=n-s*(s+1)//2
        return (s-index,index)
    def __len__(self):return self.count
    def _index(self,i):
        if i<0:i+=self.count
        if not 0<=i<self.count:raise IndexError(i)
        return i
    def __getitem__(self,i):
        if isinstance(i,slice):return [self[j] for j in range(*i.indices(self.count))]
        i=self._index(i);at=i*self.width
        return self.decode(int.from_bytes(self.data[at:at+self.width],'little'))
    def __setitem__(self,i,value):
        if isinstance(i,slice):
            indices=range(*i.indices(self.count));values=list(value)
            assert len(indices)==len(values),'reference shape cannot change through a slice'
            for j,v in zip(indices,values):self[j]=v
            return
        i=self._index(i);n=self.encode(value);width=max(1,(n.bit_length()+7)//8)
        if width>self.width:
            old=self.data;previous=self.width
            self.data=bytearray(self.count*width);self.width=width
            for j in range(self.count):self.data[j*width:j*width+previous]=old[j*previous:(j+1)*previous]
        at=i*self.width;self.data[at:at+self.width]=n.to_bytes(self.width,'little')
    def __delitem__(self,i):raise TypeError('retained reference fields cannot be deleted')
    def insert(self,i,value):
        assert i==self.count
        self.data.extend(b'\0'*self.width);self.count+=1;self[i]=value


class RefView(MutableSequence):
    __slots__=('storage','handle','count')
    def __init__(self,storage,handle,count):self.storage=storage;self.handle=handle;self.count=count
    def __len__(self):return self.count
    def __getitem__(self,i):
        if isinstance(i,slice):return [self[j] for j in range(*i.indices(self.count))]
        if i<0:i+=self.count
        if not 0<=i<self.count:raise IndexError(i)
        return self.storage.read(self.handle,i)
    def __setitem__(self,i,value):
        if isinstance(i,slice):
            indices=range(*i.indices(self.count));values=list(value)
            assert len(indices)==len(values)
            for j,v in zip(indices,values):self[j]=v
            return
        if i<0:i+=self.count
        if not 0<=i<self.count:raise IndexError(i)
        self.storage.write(self.handle,i,value)
    def __delitem__(self,i):raise TypeError('fixed source shape')
    def insert(self,i,v):raise TypeError('fixed source shape')


# These subclasses preserve the consumers' existing type/identity tests.
# All semantic fields below reside in PackedState; a proxy is canonical for
# its retained owner/slot. There is no reference-engine fallback or shadow value.
from native_install import Block, Source
from lab import Obj

class PackedBlock(Block):
    __slots__=('_storage','_handle','_refs')
    def __init__(self,storage,owner,refs,data,width,grant=0):
        self._storage=storage
        self._handle=storage.state.allocate_record(owner,(len(refs)+1,width),100)
        storage.proxies[self._handle]=self
        self._refs=RefView(storage,self._handle,len(refs))
        slot=storage.state.slot(self._handle);slot.bits[:]=data;slot.grant=grant
        self._refs[:]=refs
    @property
    def refs(self):return self._refs
    @property
    def data(self):return self._storage.state.slot(self._handle).bits
    @property
    def width(self):return self._storage.state.slot(self._handle).shape[1]
    @property
    def next(self):return self._storage.read(self._handle,len(self._refs))
    @next.setter
    def next(self,value):self._storage.write(self._handle,len(self._refs),value)
    @property
    def grant(self):return self._storage.state.slot(self._handle).grant
    @grant.setter
    def grant(self,value):self._storage.state.grant(self._handle,value)
    def bit(self,i):
        self._storage.counts['bit_reads']+=1
        return super().bit(i)
    def put(self,i,b):
        self._storage.counts['bit_writes']+=1
        super().put(i,b)

class PackedSource(Source):
    __slots__=('_storage','_handle')
    def __init__(self,storage,owner,block):
        self._storage=storage;self._handle=storage.state.allocate_record(owner,(3,0),101)
        storage.proxies[self._handle]=self;self.block=block
    def field(index):
        def get(self):return self._storage.read(self._handle,index)
        def put(self,value):self._storage.write(self._handle,index,value)
        return property(get,put)
    block=field(0);next=field(1);forward=field(2)
    @property
    def grant(self):return self._storage.state.slot(self._handle).grant
    @grant.setter
    def grant(self,value):self._storage.state.grant(self._handle,value)

class PackedFields(MutableMapping):
    __slots__=('obj','keys')
    def __init__(self,obj):self.obj=obj;self.keys={}
    def __len__(self):return len(self.keys)
    def __iter__(self):return iter(self.keys)
    def __getitem__(self,key):
        index=self.keys[key];obj=self.obj;storage=obj._storage;h=obj._handle
        if key=='bit':
            storage.counts['bit_reads']+=1
            value=storage.state.slot(h).bits[0]
            return None if value==2 else value
        return storage.read(h,index)
    def __setitem__(self,key,value):
        obj=self.obj;storage=obj._storage;h=obj._handle;slot=storage.state.slot(h)
        if key not in self.keys:
            if key=='bit':
                self.keys[key]=-1;slot.bits.append(0);slot.shape=(slot.shape[0],1)
            else:
                self.keys[key]=len(slot.refs);slot.refs.append(None);slot.shape=(len(slot.refs),slot.shape[1])
        if key=='bit':
            assert value in (0,1,None)
            slot.bits[0]=2 if value is None else value;storage.counts['bit_writes']+=1
        else:storage.write(h,self.keys[key],value)
    def __delitem__(self,key):raise TypeError('retained fields cannot be deleted')

class PackedObj(Obj):
    __slots__=('_storage','_handle')
    def __init__(self,storage,owner,kind):
        self._storage=storage;self._handle=storage.state.allocate_record(owner,(0,0),102)
        storage.proxies[self._handle]=self
        self.kind=kind;self.fields=PackedFields(self);self.writes=0

class PackedStorage:
    """Consumed by the real installer/NativeEngine and source/connection Arena.

    The native/installer control objects and C initializer remain their original
    supplied finite graph. The bank holds the exact Block/Source/Obj fields;
    it changes neither the source ledger nor owner widths and census functions.
    """
    def __init__(self):
        self.state=PackedState(1,capture=False);self.proxies={};self.counts=Counter()
        self.installation_owner=0;self.next_owner=1
    def owner(self):
        owner=self.next_owner;self.next_owner+=1
        self.state.owners[owner]=[];self.state.width[owner]=1;self.state.free[owner]=[]
        return owner
    def handle(self,value):
        if value is None:return None
        assert value._storage is self,'cross-storage reference is not a source identity'
        return value._handle
    def read(self,h,i):
        self.counts['reference_reads']+=1
        value=self.state.read_ref(h,i)
        return None if value is None else self.proxies[value]
    def write(self,h,i,value):
        self.counts['reference_writes']+=1
        self.state.write_ref(h,i,self.handle(value))
    def block(self,refs,data,width):return PackedBlock(self,self.installation_owner,refs,data,width)
    def obj(self,owner,kind):return PackedObj(self,owner,kind)
    def import_sources(self,sources):
        # Same public preparation graph, represented with a distinct slot for
        # each old Block and Source. This is preparation, never paid relocation.
        blocks={}
        for s in sources:
            if s.block not in blocks:blocks[s.block]=self.block([None]*len(s.block.refs),s.block.data,s.block.width)
        mapping={s:PackedSource(self,self.installation_owner,blocks[s.block]) for s in sources}
        for old,new in mapping.items():
            new.grant=old.grant;new.block.grant=old.block.grant
            new.next=None if old.next is None else mapping[old.next]
            assert old.forward is None,'import only the uninstalled public source'
            new.block.refs[:]=[None if x is None else mapping[x] for x in old.block.refs]
            assert old.block.next is None
        return [mapping[s] for s in sources],mapping
    def metrics(self):
        return {**self.counts,'slots':len(self.proxies),
            'reference_payload_bytes':sum(len(s.refs.data) for a in self.state.owners.values() for s in a),
            'bit_payload_bytes':sum(len(s.bits) for a in self.state.owners.values() for s in a),
            'owner_count':len(self.state.owners),
            'host_scope':'payload only; canonical proxies, schemas, runtime, controllers and observer maps are additional'}
