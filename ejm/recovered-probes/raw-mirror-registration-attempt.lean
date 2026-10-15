/- Transient source-bound registration candidate; not compiled or delivered.
The original telescope retains every positive dimension n and row count N.
For fixed (n,N), State = MirrorRow × RawRow and Output =
(Partial (Fin n) (Fin 4)) × Int. The proposed guard is checkRawMirror;
left readout is (partialEvent mirror.hash, mirror.lower), and right is
(original row.assignment, original row.lower).
The existing GuardedEqualityRegistrationTemplates.guardedEqRealization uses
one nondependent Output Y for its whole State. It therefore cannot literally
carry the varying original Fin n output across the complete telescope.
DependentFamily can represent this output dependency, but original source
selection, exact source/Law equivalence, and complete guard/left/right role
variation, sensitivity and observational dependence remain to be constructed.
No impossibility is claimed for DependentFamily. No registration status is
claimed, and no guard-only Boolean surrogate is substituted for the source.
-/
import D5.S3.Quantum.Entanglement.EJM.EJMLocalCacheData
import Reg.Support.DependentFamily
import Reg.Support.GuardedEqualityRegistrations

open D5.S0.Certificates.FiniteSymmetryBoundReplay
open D5.S3.Quantum.Entanglement.EJM.EJMPartialBoundData
open D5.S3.Quantum.Entanglement.EJM.EJMLocalCacheData

private abbrev RawMirrorState := MirrorRow × RawRow
private def originalPair {n N : Nat} (hn : 0 < n) (hN : 0 < N)
    (state : RawMirrorState) : Partial (Fin n) (Fin 4) × Int :=
  ((row hn hN state.2).assignment, (row hn hN state.2).lower)
private def mirrorPair (n : Nat) (state : RawMirrorState) :
    Partial (Fin n) (Fin 4) × Int :=
  (partialEvent state.1.hash, state.1.lower)
private def originalGuard (state : RawMirrorState) : Bool :=
  checkRawMirror state.1 state.2
