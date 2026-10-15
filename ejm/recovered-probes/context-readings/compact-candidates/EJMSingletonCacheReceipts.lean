/- GID: D5/S3/Quantum/Entanglement/EJM/EJMSingletonCacheReceipts
   generality: I
   mirror-B: D5/B/S3/Quantum/Entanglement/EJM/EJMSingletonCacheReceipts
   mirror-E: none(waiver:actual-exact-certificate-cache)
   anchors: []
   utility: kind=numeric-reduction; basis=consumer=D5/S3/Quantum/Entanglement/EJM/EJMExactCertificate.result; premises=D5/S3/Quantum/Entanglement/EJM/EJMSingletonCacheReceipts.accepted
   digest: Full finite-domain kernel receipts bind every local mirror and unknown-coordinate transition to the original rows. -/

/- Part of the same complete actual EJM
certificate, not a standalone conditional or numerical delivery. All literal
mirrors, actions and transition receipts require kernel acceptance before use.
The original score, scale, all actual assignments and four value branches are
preserved. Producer checks below are not mathematical premises. -/

import D5.S3.Quantum.Entanglement.EJM.EJMSingletonCacheData

set_option autoImplicit false
set_option relaxedAutoImplicit false
set_option maxRecDepth 100000
-- Async elaboration retained a measured 32 GiB heap for these finite proofs.
-- Finish each unchanged kernel receipt before beginning the next one.
set_option Elab.async false

namespace D5.S3.Quantum.Entanglement.EJM.EJMSingletonCacheReceipts
open D5.S3.Quantum.Entanglement.EJM.EJMPartialBoundData
open D5.S3.Quantum.Entanglement.EJM.EJMPartialBoundTransport
open D5.S3.Quantum.Entanglement.EJM.EJMSingletonCacheData

set_option maxHeartbeats 100000000
set_option trace.profiler true
set_option linter.style.longLine false

private theorem mirror_0 : ∀ k : Fin 128, checkMirrorFields ⟨0 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_128 : ∀ k : Fin 128, checkMirrorFields ⟨128 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_256 : ∀ k : Fin 128, checkMirrorFields ⟨256 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_384 : ∀ k : Fin 128, checkMirrorFields ⟨384 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_512 : ∀ k : Fin 128, checkMirrorFields ⟨512 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_640 : ∀ k : Fin 128, checkMirrorFields ⟨640 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_768 : ∀ k : Fin 128, checkMirrorFields ⟨768 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_896 : ∀ k : Fin 128, checkMirrorFields ⟨896 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_1024 : ∀ k : Fin 128, checkMirrorFields ⟨1024 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_1152 : ∀ k : Fin 128, checkMirrorFields ⟨1152 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_1280 : ∀ k : Fin 128, checkMirrorFields ⟨1280 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_1408 : ∀ k : Fin 128, checkMirrorFields ⟨1408 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_1536 : ∀ k : Fin 128, checkMirrorFields ⟨1536 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_1664 : ∀ k : Fin 128, checkMirrorFields ⟨1664 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_1792 : ∀ k : Fin 128, checkMirrorFields ⟨1792 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_1920 : ∀ k : Fin 128, checkMirrorFields ⟨1920 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_2048 : ∀ k : Fin 128, checkMirrorFields ⟨2048 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_2176 : ∀ k : Fin 128, checkMirrorFields ⟨2176 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_2304 : ∀ k : Fin 128, checkMirrorFields ⟨2304 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_2432 : ∀ k : Fin 128, checkMirrorFields ⟨2432 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_2560 : ∀ k : Fin 128, checkMirrorFields ⟨2560 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_2688 : ∀ k : Fin 128, checkMirrorFields ⟨2688 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_2816 : ∀ k : Fin 128, checkMirrorFields ⟨2816 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_2944 : ∀ k : Fin 128, checkMirrorFields ⟨2944 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_3072 : ∀ k : Fin 128, checkMirrorFields ⟨3072 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_3200 : ∀ k : Fin 128, checkMirrorFields ⟨3200 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_3328 : ∀ k : Fin 128, checkMirrorFields ⟨3328 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_3456 : ∀ k : Fin 128, checkMirrorFields ⟨3456 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_3584 : ∀ k : Fin 128, checkMirrorFields ⟨3584 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_3712 : ∀ k : Fin 128, checkMirrorFields ⟨3712 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_3840 : ∀ k : Fin 128, checkMirrorFields ⟨3840 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_3968 : ∀ k : Fin 128, checkMirrorFields ⟨3968 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_4096 : ∀ k : Fin 128, checkMirrorFields ⟨4096 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_4224 : ∀ k : Fin 128, checkMirrorFields ⟨4224 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_4352 : ∀ k : Fin 128, checkMirrorFields ⟨4352 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_4480 : ∀ k : Fin 128, checkMirrorFields ⟨4480 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_4608 : ∀ k : Fin 128, checkMirrorFields ⟨4608 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_4736 : ∀ k : Fin 128, checkMirrorFields ⟨4736 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_4864 : ∀ k : Fin 128, checkMirrorFields ⟨4864 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_4992 : ∀ k : Fin 128, checkMirrorFields ⟨4992 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_5120 : ∀ k : Fin 128, checkMirrorFields ⟨5120 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_5248 : ∀ k : Fin 128, checkMirrorFields ⟨5248 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_5376 : ∀ k : Fin 128, checkMirrorFields ⟨5376 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_5504 : ∀ k : Fin 128, checkMirrorFields ⟨5504 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_5632 : ∀ k : Fin 128, checkMirrorFields ⟨5632 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_5760 : ∀ k : Fin 128, checkMirrorFields ⟨5760 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_5888 : ∀ k : Fin 128, checkMirrorFields ⟨5888 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_6016 : ∀ k : Fin 128, checkMirrorFields ⟨6016 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_6144 : ∀ k : Fin 128, checkMirrorFields ⟨6144 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_6272 : ∀ k : Fin 128, checkMirrorFields ⟨6272 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_6400 : ∀ k : Fin 128, checkMirrorFields ⟨6400 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_6528 : ∀ k : Fin 128, checkMirrorFields ⟨6528 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_6656 : ∀ k : Fin 128, checkMirrorFields ⟨6656 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_6784 : ∀ k : Fin 128, checkMirrorFields ⟨6784 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_6912 : ∀ k : Fin 128, checkMirrorFields ⟨6912 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_7040 : ∀ k : Fin 128, checkMirrorFields ⟨7040 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_7168 : ∀ k : Fin 128, checkMirrorFields ⟨7168 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_7296 : ∀ k : Fin 128, checkMirrorFields ⟨7296 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_7424 : ∀ k : Fin 128, checkMirrorFields ⟨7424 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_7552 : ∀ k : Fin 128, checkMirrorFields ⟨7552 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_7680 : ∀ k : Fin 128, checkMirrorFields ⟨7680 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_7808 : ∀ k : Fin 128, checkMirrorFields ⟨7808 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_7936 : ∀ k : Fin 128, checkMirrorFields ⟨7936 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_8064 : ∀ k : Fin 128, checkMirrorFields ⟨8064 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_8192 : ∀ k : Fin 128, checkMirrorFields ⟨8192 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_8320 : ∀ k : Fin 128, checkMirrorFields ⟨8320 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_8448 : ∀ k : Fin 128, checkMirrorFields ⟨8448 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_8576 : ∀ k : Fin 113, checkMirrorFields ⟨8576 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_0 : ∀ k : Fin 128, checkTransitions ⟨0 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_128 : ∀ k : Fin 128, checkTransitions ⟨128 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_256 : ∀ k : Fin 128, checkTransitions ⟨256 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_384 : ∀ k : Fin 128, checkTransitions ⟨384 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_512 : ∀ k : Fin 128, checkTransitions ⟨512 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_640 : ∀ k : Fin 128, checkTransitions ⟨640 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_768 : ∀ k : Fin 128, checkTransitions ⟨768 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_896 : ∀ k : Fin 128, checkTransitions ⟨896 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_1024 : ∀ k : Fin 128, checkTransitions ⟨1024 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_1152 : ∀ k : Fin 128, checkTransitions ⟨1152 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_1280 : ∀ k : Fin 128, checkTransitions ⟨1280 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_1408 : ∀ k : Fin 128, checkTransitions ⟨1408 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_1536 : ∀ k : Fin 128, checkTransitions ⟨1536 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_1664 : ∀ k : Fin 128, checkTransitions ⟨1664 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_1792 : ∀ k : Fin 128, checkTransitions ⟨1792 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_1920 : ∀ k : Fin 128, checkTransitions ⟨1920 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_2048 : ∀ k : Fin 128, checkTransitions ⟨2048 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_2176 : ∀ k : Fin 128, checkTransitions ⟨2176 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_2304 : ∀ k : Fin 128, checkTransitions ⟨2304 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_2432 : ∀ k : Fin 128, checkTransitions ⟨2432 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_2560 : ∀ k : Fin 128, checkTransitions ⟨2560 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_2688 : ∀ k : Fin 128, checkTransitions ⟨2688 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_2816 : ∀ k : Fin 128, checkTransitions ⟨2816 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_2944 : ∀ k : Fin 128, checkTransitions ⟨2944 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_3072 : ∀ k : Fin 128, checkTransitions ⟨3072 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_3200 : ∀ k : Fin 128, checkTransitions ⟨3200 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_3328 : ∀ k : Fin 128, checkTransitions ⟨3328 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_3456 : ∀ k : Fin 128, checkTransitions ⟨3456 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_3584 : ∀ k : Fin 128, checkTransitions ⟨3584 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_3712 : ∀ k : Fin 128, checkTransitions ⟨3712 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_3840 : ∀ k : Fin 128, checkTransitions ⟨3840 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_3968 : ∀ k : Fin 128, checkTransitions ⟨3968 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_4096 : ∀ k : Fin 128, checkTransitions ⟨4096 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_4224 : ∀ k : Fin 128, checkTransitions ⟨4224 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_4352 : ∀ k : Fin 128, checkTransitions ⟨4352 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_4480 : ∀ k : Fin 128, checkTransitions ⟨4480 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_4608 : ∀ k : Fin 128, checkTransitions ⟨4608 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_4736 : ∀ k : Fin 128, checkTransitions ⟨4736 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_4864 : ∀ k : Fin 128, checkTransitions ⟨4864 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_4992 : ∀ k : Fin 128, checkTransitions ⟨4992 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_5120 : ∀ k : Fin 128, checkTransitions ⟨5120 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_5248 : ∀ k : Fin 128, checkTransitions ⟨5248 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_5376 : ∀ k : Fin 128, checkTransitions ⟨5376 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_5504 : ∀ k : Fin 128, checkTransitions ⟨5504 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_5632 : ∀ k : Fin 128, checkTransitions ⟨5632 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_5760 : ∀ k : Fin 128, checkTransitions ⟨5760 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_5888 : ∀ k : Fin 128, checkTransitions ⟨5888 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_6016 : ∀ k : Fin 128, checkTransitions ⟨6016 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_6144 : ∀ k : Fin 128, checkTransitions ⟨6144 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_6272 : ∀ k : Fin 128, checkTransitions ⟨6272 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_6400 : ∀ k : Fin 128, checkTransitions ⟨6400 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_6528 : ∀ k : Fin 128, checkTransitions ⟨6528 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_6656 : ∀ k : Fin 128, checkTransitions ⟨6656 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_6784 : ∀ k : Fin 128, checkTransitions ⟨6784 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_6912 : ∀ k : Fin 128, checkTransitions ⟨6912 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_7040 : ∀ k : Fin 128, checkTransitions ⟨7040 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_7168 : ∀ k : Fin 128, checkTransitions ⟨7168 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_7296 : ∀ k : Fin 128, checkTransitions ⟨7296 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_7424 : ∀ k : Fin 128, checkTransitions ⟨7424 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_7552 : ∀ k : Fin 128, checkTransitions ⟨7552 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_7680 : ∀ k : Fin 128, checkTransitions ⟨7680 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_7808 : ∀ k : Fin 128, checkTransitions ⟨7808 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_7936 : ∀ k : Fin 128, checkTransitions ⟨7936 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_8064 : ∀ k : Fin 128, checkTransitions ⟨8064 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_8192 : ∀ k : Fin 128, checkTransitions ⟨8192 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_8320 : ∀ k : Fin 128, checkTransitions ⟨8320 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_8448 : ∀ k : Fin 128, checkTransitions ⟨8448 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_8576 : ∀ k : Fin 113, checkTransitions ⟨8576 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel

private def chunkLength (b : Fin 68) : ℕ :=
  min 128 (8689 - 128 * b.val)

private def chunkIndex (b : Fin 68) (k : Fin (chunkLength b)) : Fin 8689 :=
  ⟨128 * b.val + k.val, by
    have hb := b.isLt
    have hk := k.isLt
    dsimp [chunkLength] at hk
    omega⟩

private theorem mirror_blocks : ∀ b : Fin 68, ∀ k : Fin (chunkLength b),
    checkMirrorFields (chunkIndex b k) = true :=
  Fin.cases mirror_0 (Fin.cases mirror_128 (Fin.cases mirror_256 (Fin.cases mirror_384 (Fin.cases mirror_512 (Fin.cases mirror_640 (Fin.cases mirror_768 (Fin.cases mirror_896 (Fin.cases mirror_1024 (Fin.cases mirror_1152 (Fin.cases mirror_1280 (Fin.cases mirror_1408 (Fin.cases mirror_1536 (Fin.cases mirror_1664 (Fin.cases mirror_1792 (Fin.cases mirror_1920 (
  Fin.cases mirror_2048 (Fin.cases mirror_2176 (Fin.cases mirror_2304 (Fin.cases mirror_2432 (Fin.cases mirror_2560 (Fin.cases mirror_2688 (Fin.cases mirror_2816 (Fin.cases mirror_2944 (Fin.cases mirror_3072 (Fin.cases mirror_3200 (Fin.cases mirror_3328 (Fin.cases mirror_3456 (Fin.cases mirror_3584 (Fin.cases mirror_3712 (Fin.cases mirror_3840 (Fin.cases mirror_3968 (
  Fin.cases mirror_4096 (Fin.cases mirror_4224 (Fin.cases mirror_4352 (Fin.cases mirror_4480 (Fin.cases mirror_4608 (Fin.cases mirror_4736 (Fin.cases mirror_4864 (Fin.cases mirror_4992 (Fin.cases mirror_5120 (Fin.cases mirror_5248 (Fin.cases mirror_5376 (Fin.cases mirror_5504 (Fin.cases mirror_5632 (Fin.cases mirror_5760 (Fin.cases mirror_5888 (Fin.cases mirror_6016 (
  Fin.cases mirror_6144 (Fin.cases mirror_6272 (Fin.cases mirror_6400 (Fin.cases mirror_6528 (Fin.cases mirror_6656 (Fin.cases mirror_6784 (Fin.cases mirror_6912 (Fin.cases mirror_7040 (Fin.cases mirror_7168 (Fin.cases mirror_7296 (Fin.cases mirror_7424 (Fin.cases mirror_7552 (Fin.cases mirror_7680 (Fin.cases mirror_7808 (Fin.cases mirror_7936 (Fin.cases mirror_8064 (
  Fin.cases mirror_8192 (Fin.cases mirror_8320 (Fin.cases mirror_8448 (Fin.cases mirror_8576 (
  (fun b : Fin 0 => Fin.elim0 b)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

private theorem transitions_blocks : ∀ b : Fin 68, ∀ k : Fin (chunkLength b),
    checkTransitions (chunkIndex b k) = true :=
  Fin.cases transitions_0 (Fin.cases transitions_128 (Fin.cases transitions_256 (Fin.cases transitions_384 (Fin.cases transitions_512 (Fin.cases transitions_640 (Fin.cases transitions_768 (Fin.cases transitions_896 (Fin.cases transitions_1024 (Fin.cases transitions_1152 (Fin.cases transitions_1280 (Fin.cases transitions_1408 (Fin.cases transitions_1536 (Fin.cases transitions_1664 (Fin.cases transitions_1792 (Fin.cases transitions_1920 (
  Fin.cases transitions_2048 (Fin.cases transitions_2176 (Fin.cases transitions_2304 (Fin.cases transitions_2432 (Fin.cases transitions_2560 (Fin.cases transitions_2688 (Fin.cases transitions_2816 (Fin.cases transitions_2944 (Fin.cases transitions_3072 (Fin.cases transitions_3200 (Fin.cases transitions_3328 (Fin.cases transitions_3456 (Fin.cases transitions_3584 (Fin.cases transitions_3712 (Fin.cases transitions_3840 (Fin.cases transitions_3968 (
  Fin.cases transitions_4096 (Fin.cases transitions_4224 (Fin.cases transitions_4352 (Fin.cases transitions_4480 (Fin.cases transitions_4608 (Fin.cases transitions_4736 (Fin.cases transitions_4864 (Fin.cases transitions_4992 (Fin.cases transitions_5120 (Fin.cases transitions_5248 (Fin.cases transitions_5376 (Fin.cases transitions_5504 (Fin.cases transitions_5632 (Fin.cases transitions_5760 (Fin.cases transitions_5888 (Fin.cases transitions_6016 (
  Fin.cases transitions_6144 (Fin.cases transitions_6272 (Fin.cases transitions_6400 (Fin.cases transitions_6528 (Fin.cases transitions_6656 (Fin.cases transitions_6784 (Fin.cases transitions_6912 (Fin.cases transitions_7040 (Fin.cases transitions_7168 (Fin.cases transitions_7296 (Fin.cases transitions_7424 (Fin.cases transitions_7552 (Fin.cases transitions_7680 (Fin.cases transitions_7808 (Fin.cases transitions_7936 (Fin.cases transitions_8064 (
  Fin.cases transitions_8192 (Fin.cases transitions_8320 (Fin.cases transitions_8448 (Fin.cases transitions_8576 (
  (fun b : Fin 0 => Fin.elim0 b)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

private theorem mirror_fields_accepted (j : Fin 8689) : checkMirrorFields j = true := by
  let b : Fin 68 := ⟨j.val / 128, by have hj := j.isLt; omega⟩
  let k : Fin (chunkLength b) := ⟨j.val % 128, by
    have hj := j.isLt
    dsimp [chunkLength, b]
    omega⟩
  have he : chunkIndex b k = j := by
    apply Fin.ext
    dsimp [chunkIndex, b, k]
    omega
  simpa only [he] using mirror_blocks b k

theorem mirror_accepted (j : Fin 8689) : checkMirror j = true := by
  apply decide_eq_true
  exact EJMLocalCacheData.checked_raw_mirror_eq (by decide) (by decide)
    (mirror j) (singletonTree.read 9 j.val) (mirror_fields_accepted j)


theorem transitions_accepted (j : Fin 8689) : checkTransitions j = true := by
  let b : Fin 68 := ⟨j.val / 128, by have hj := j.isLt; omega⟩
  let k : Fin (chunkLength b) := ⟨j.val % 128, by
    have hj := j.isLt
    dsimp [chunkLength, b]
    omega⟩
  have he : chunkIndex b k = j := by
    apply Fin.ext
    dsimp [chunkIndex, b, k]
    omega
  simpa only [he] using transitions_blocks b k

theorem assignment_eq (j : Fin 8689) :
    assignment j = (singletonRows j).assignment :=
  (of_decide_eq_true (mirror_accepted j)).1

theorem lower_eq (j : Fin 8689) : lower j = (singletonRows j).lower :=
  (of_decide_eq_true (mirror_accepted j)).2

/-- This consequence uses the original row assignments on both sides. -/
theorem transition_accepted (j : Fin 8689) (i : Fin 9) (a : Fin 4)
    (unknown : (singletonRows j).assignment i = none) :
    D5.S0.Certificates.FiniteLocalRowCache.checkTransition singletonPartialAct
      singletonRows transition j i a = true := by
  have hu : assignment j i = none := by rw [assignment_eq, unknown]
  have h := of_decide_eq_true (transitions_accepted j) i hu a
  rw [assignment_eq, assignment_eq] at h
  exact decide_eq_true h

theorem accepted : (∀ j, checkMirror j = true) ∧ (∀ j, checkTransitions j = true) :=
  ⟨mirror_accepted, transitions_accepted⟩

#print axioms assignment_eq
#print axioms lower_eq
#print axioms transition_accepted

end D5.S3.Quantum.Entanglement.EJM.EJMSingletonCacheReceipts
