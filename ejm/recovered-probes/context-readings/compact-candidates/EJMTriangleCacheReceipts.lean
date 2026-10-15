/- GID: D5/S3/Quantum/Entanglement/EJM/EJMTriangleCacheReceipts
   generality: I
   mirror-B: D5/B/S3/Quantum/Entanglement/EJM/EJMTriangleCacheReceipts
   mirror-E: none(waiver:actual-exact-certificate-cache)
   anchors: []
   utility: kind=numeric-reduction; basis=consumer=D5/S3/Quantum/Entanglement/EJM/EJMExactCertificate.result; premises=D5/S3/Quantum/Entanglement/EJM/EJMTriangleCacheReceipts.accepted
   digest: Full finite-domain kernel receipts bind every local mirror and unknown-coordinate transition to the original rows. -/

/- Part of the same complete actual EJM certificate,
not a standalone conditional or numerical delivery. All literal
mirrors, actions and transition receipts require kernel acceptance before use.
The original score, scale, all actual assignments and four value branches are
preserved. Producer checks below are not mathematical premises. -/

import D5.S3.Quantum.Entanglement.EJM.EJMTriangleCacheData

set_option autoImplicit false
set_option relaxedAutoImplicit false
set_option maxRecDepth 100000
-- Bounded kernel reductions preserve every original row and transition.
-- Complete each receipt before beginning the next finite block.
set_option Elab.async false

namespace D5.S3.Quantum.Entanglement.EJM.EJMTriangleCacheReceipts
open D5.S3.Quantum.Entanglement.EJM.EJMPartialBoundData
open D5.S3.Quantum.Entanglement.EJM.EJMPartialBoundTransport
open D5.S3.Quantum.Entanglement.EJM.EJMTriangleCacheData

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
private theorem mirror_8576 : ∀ k : Fin 128, checkMirrorFields ⟨8576 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_8704 : ∀ k : Fin 128, checkMirrorFields ⟨8704 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_8832 : ∀ k : Fin 128, checkMirrorFields ⟨8832 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_8960 : ∀ k : Fin 128, checkMirrorFields ⟨8960 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_9088 : ∀ k : Fin 128, checkMirrorFields ⟨9088 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_9216 : ∀ k : Fin 128, checkMirrorFields ⟨9216 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_9344 : ∀ k : Fin 128, checkMirrorFields ⟨9344 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_9472 : ∀ k : Fin 128, checkMirrorFields ⟨9472 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_9600 : ∀ k : Fin 128, checkMirrorFields ⟨9600 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_9728 : ∀ k : Fin 128, checkMirrorFields ⟨9728 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_9856 : ∀ k : Fin 128, checkMirrorFields ⟨9856 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_9984 : ∀ k : Fin 128, checkMirrorFields ⟨9984 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_10112 : ∀ k : Fin 128, checkMirrorFields ⟨10112 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_10240 : ∀ k : Fin 128, checkMirrorFields ⟨10240 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_10368 : ∀ k : Fin 128, checkMirrorFields ⟨10368 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_10496 : ∀ k : Fin 128, checkMirrorFields ⟨10496 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_10624 : ∀ k : Fin 128, checkMirrorFields ⟨10624 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_10752 : ∀ k : Fin 128, checkMirrorFields ⟨10752 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_10880 : ∀ k : Fin 128, checkMirrorFields ⟨10880 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_11008 : ∀ k : Fin 128, checkMirrorFields ⟨11008 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_11136 : ∀ k : Fin 128, checkMirrorFields ⟨11136 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_11264 : ∀ k : Fin 128, checkMirrorFields ⟨11264 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_11392 : ∀ k : Fin 128, checkMirrorFields ⟨11392 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_11520 : ∀ k : Fin 128, checkMirrorFields ⟨11520 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_11648 : ∀ k : Fin 128, checkMirrorFields ⟨11648 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_11776 : ∀ k : Fin 128, checkMirrorFields ⟨11776 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_11904 : ∀ k : Fin 128, checkMirrorFields ⟨11904 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_12032 : ∀ k : Fin 128, checkMirrorFields ⟨12032 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_12160 : ∀ k : Fin 128, checkMirrorFields ⟨12160 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_12288 : ∀ k : Fin 128, checkMirrorFields ⟨12288 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_12416 : ∀ k : Fin 128, checkMirrorFields ⟨12416 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_12544 : ∀ k : Fin 128, checkMirrorFields ⟨12544 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_12672 : ∀ k : Fin 128, checkMirrorFields ⟨12672 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_12800 : ∀ k : Fin 128, checkMirrorFields ⟨12800 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_12928 : ∀ k : Fin 128, checkMirrorFields ⟨12928 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_13056 : ∀ k : Fin 128, checkMirrorFields ⟨13056 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_13184 : ∀ k : Fin 128, checkMirrorFields ⟨13184 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_13312 : ∀ k : Fin 128, checkMirrorFields ⟨13312 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_13440 : ∀ k : Fin 128, checkMirrorFields ⟨13440 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_13568 : ∀ k : Fin 128, checkMirrorFields ⟨13568 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_13696 : ∀ k : Fin 128, checkMirrorFields ⟨13696 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_13824 : ∀ k : Fin 128, checkMirrorFields ⟨13824 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_13952 : ∀ k : Fin 128, checkMirrorFields ⟨13952 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_14080 : ∀ k : Fin 128, checkMirrorFields ⟨14080 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_14208 : ∀ k : Fin 128, checkMirrorFields ⟨14208 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_14336 : ∀ k : Fin 128, checkMirrorFields ⟨14336 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_14464 : ∀ k : Fin 128, checkMirrorFields ⟨14464 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_14592 : ∀ k : Fin 128, checkMirrorFields ⟨14592 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_14720 : ∀ k : Fin 128, checkMirrorFields ⟨14720 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_14848 : ∀ k : Fin 128, checkMirrorFields ⟨14848 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_14976 : ∀ k : Fin 128, checkMirrorFields ⟨14976 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_15104 : ∀ k : Fin 128, checkMirrorFields ⟨15104 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_15232 : ∀ k : Fin 128, checkMirrorFields ⟨15232 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_15360 : ∀ k : Fin 128, checkMirrorFields ⟨15360 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_15488 : ∀ k : Fin 128, checkMirrorFields ⟨15488 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_15616 : ∀ k : Fin 128, checkMirrorFields ⟨15616 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_15744 : ∀ k : Fin 128, checkMirrorFields ⟨15744 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_15872 : ∀ k : Fin 128, checkMirrorFields ⟨15872 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_16000 : ∀ k : Fin 128, checkMirrorFields ⟨16000 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_16128 : ∀ k : Fin 128, checkMirrorFields ⟨16128 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_16256 : ∀ k : Fin 128, checkMirrorFields ⟨16256 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_16384 : ∀ k : Fin 128, checkMirrorFields ⟨16384 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_16512 : ∀ k : Fin 128, checkMirrorFields ⟨16512 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_16640 : ∀ k : Fin 128, checkMirrorFields ⟨16640 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_16768 : ∀ k : Fin 128, checkMirrorFields ⟨16768 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_16896 : ∀ k : Fin 128, checkMirrorFields ⟨16896 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_17024 : ∀ k : Fin 128, checkMirrorFields ⟨17024 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_17152 : ∀ k : Fin 128, checkMirrorFields ⟨17152 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_17280 : ∀ k : Fin 128, checkMirrorFields ⟨17280 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_17408 : ∀ k : Fin 128, checkMirrorFields ⟨17408 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_17536 : ∀ k : Fin 128, checkMirrorFields ⟨17536 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_17664 : ∀ k : Fin 128, checkMirrorFields ⟨17664 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_17792 : ∀ k : Fin 128, checkMirrorFields ⟨17792 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_17920 : ∀ k : Fin 128, checkMirrorFields ⟨17920 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_18048 : ∀ k : Fin 128, checkMirrorFields ⟨18048 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_18176 : ∀ k : Fin 128, checkMirrorFields ⟨18176 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_18304 : ∀ k : Fin 128, checkMirrorFields ⟨18304 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_18432 : ∀ k : Fin 128, checkMirrorFields ⟨18432 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_18560 : ∀ k : Fin 128, checkMirrorFields ⟨18560 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_18688 : ∀ k : Fin 128, checkMirrorFields ⟨18688 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_18816 : ∀ k : Fin 128, checkMirrorFields ⟨18816 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_18944 : ∀ k : Fin 128, checkMirrorFields ⟨18944 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_19072 : ∀ k : Fin 128, checkMirrorFields ⟨19072 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_19200 : ∀ k : Fin 128, checkMirrorFields ⟨19200 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_19328 : ∀ k : Fin 128, checkMirrorFields ⟨19328 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_19456 : ∀ k : Fin 128, checkMirrorFields ⟨19456 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_19584 : ∀ k : Fin 128, checkMirrorFields ⟨19584 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_19712 : ∀ k : Fin 128, checkMirrorFields ⟨19712 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_19840 : ∀ k : Fin 128, checkMirrorFields ⟨19840 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_19968 : ∀ k : Fin 128, checkMirrorFields ⟨19968 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_20096 : ∀ k : Fin 128, checkMirrorFields ⟨20096 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_20224 : ∀ k : Fin 128, checkMirrorFields ⟨20224 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_20352 : ∀ k : Fin 128, checkMirrorFields ⟨20352 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_20480 : ∀ k : Fin 128, checkMirrorFields ⟨20480 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_20608 : ∀ k : Fin 128, checkMirrorFields ⟨20608 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_20736 : ∀ k : Fin 128, checkMirrorFields ⟨20736 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_20864 : ∀ k : Fin 128, checkMirrorFields ⟨20864 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_20992 : ∀ k : Fin 128, checkMirrorFields ⟨20992 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_21120 : ∀ k : Fin 128, checkMirrorFields ⟨21120 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_21248 : ∀ k : Fin 128, checkMirrorFields ⟨21248 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_21376 : ∀ k : Fin 128, checkMirrorFields ⟨21376 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_21504 : ∀ k : Fin 128, checkMirrorFields ⟨21504 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_21632 : ∀ k : Fin 128, checkMirrorFields ⟨21632 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_21760 : ∀ k : Fin 128, checkMirrorFields ⟨21760 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_21888 : ∀ k : Fin 128, checkMirrorFields ⟨21888 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_22016 : ∀ k : Fin 128, checkMirrorFields ⟨22016 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_22144 : ∀ k : Fin 128, checkMirrorFields ⟨22144 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_22272 : ∀ k : Fin 128, checkMirrorFields ⟨22272 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_22400 : ∀ k : Fin 128, checkMirrorFields ⟨22400 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_22528 : ∀ k : Fin 128, checkMirrorFields ⟨22528 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_22656 : ∀ k : Fin 128, checkMirrorFields ⟨22656 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_22784 : ∀ k : Fin 128, checkMirrorFields ⟨22784 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_22912 : ∀ k : Fin 128, checkMirrorFields ⟨22912 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_23040 : ∀ k : Fin 128, checkMirrorFields ⟨23040 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_23168 : ∀ k : Fin 128, checkMirrorFields ⟨23168 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_23296 : ∀ k : Fin 128, checkMirrorFields ⟨23296 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_23424 : ∀ k : Fin 128, checkMirrorFields ⟨23424 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_23552 : ∀ k : Fin 128, checkMirrorFields ⟨23552 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_23680 : ∀ k : Fin 128, checkMirrorFields ⟨23680 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_23808 : ∀ k : Fin 128, checkMirrorFields ⟨23808 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_23936 : ∀ k : Fin 128, checkMirrorFields ⟨23936 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_24064 : ∀ k : Fin 128, checkMirrorFields ⟨24064 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_24192 : ∀ k : Fin 128, checkMirrorFields ⟨24192 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_24320 : ∀ k : Fin 128, checkMirrorFields ⟨24320 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_24448 : ∀ k : Fin 128, checkMirrorFields ⟨24448 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_24576 : ∀ k : Fin 128, checkMirrorFields ⟨24576 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_24704 : ∀ k : Fin 128, checkMirrorFields ⟨24704 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_24832 : ∀ k : Fin 128, checkMirrorFields ⟨24832 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_24960 : ∀ k : Fin 128, checkMirrorFields ⟨24960 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_25088 : ∀ k : Fin 128, checkMirrorFields ⟨25088 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_25216 : ∀ k : Fin 128, checkMirrorFields ⟨25216 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_25344 : ∀ k : Fin 128, checkMirrorFields ⟨25344 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_25472 : ∀ k : Fin 128, checkMirrorFields ⟨25472 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_25600 : ∀ k : Fin 128, checkMirrorFields ⟨25600 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_25728 : ∀ k : Fin 128, checkMirrorFields ⟨25728 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_25856 : ∀ k : Fin 128, checkMirrorFields ⟨25856 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_25984 : ∀ k : Fin 128, checkMirrorFields ⟨25984 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_26112 : ∀ k : Fin 128, checkMirrorFields ⟨26112 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_26240 : ∀ k : Fin 128, checkMirrorFields ⟨26240 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_26368 : ∀ k : Fin 128, checkMirrorFields ⟨26368 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_26496 : ∀ k : Fin 128, checkMirrorFields ⟨26496 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_26624 : ∀ k : Fin 128, checkMirrorFields ⟨26624 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_26752 : ∀ k : Fin 128, checkMirrorFields ⟨26752 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_26880 : ∀ k : Fin 128, checkMirrorFields ⟨26880 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_27008 : ∀ k : Fin 128, checkMirrorFields ⟨27008 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_27136 : ∀ k : Fin 128, checkMirrorFields ⟨27136 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_27264 : ∀ k : Fin 128, checkMirrorFields ⟨27264 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_27392 : ∀ k : Fin 128, checkMirrorFields ⟨27392 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_27520 : ∀ k : Fin 128, checkMirrorFields ⟨27520 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_27648 : ∀ k : Fin 128, checkMirrorFields ⟨27648 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_27776 : ∀ k : Fin 128, checkMirrorFields ⟨27776 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_27904 : ∀ k : Fin 128, checkMirrorFields ⟨27904 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_28032 : ∀ k : Fin 128, checkMirrorFields ⟨28032 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_28160 : ∀ k : Fin 128, checkMirrorFields ⟨28160 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_28288 : ∀ k : Fin 128, checkMirrorFields ⟨28288 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_28416 : ∀ k : Fin 128, checkMirrorFields ⟨28416 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_28544 : ∀ k : Fin 128, checkMirrorFields ⟨28544 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_28672 : ∀ k : Fin 128, checkMirrorFields ⟨28672 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_28800 : ∀ k : Fin 128, checkMirrorFields ⟨28800 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_28928 : ∀ k : Fin 128, checkMirrorFields ⟨28928 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_29056 : ∀ k : Fin 128, checkMirrorFields ⟨29056 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_29184 : ∀ k : Fin 128, checkMirrorFields ⟨29184 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_29312 : ∀ k : Fin 128, checkMirrorFields ⟨29312 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_29440 : ∀ k : Fin 128, checkMirrorFields ⟨29440 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_29568 : ∀ k : Fin 128, checkMirrorFields ⟨29568 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_29696 : ∀ k : Fin 128, checkMirrorFields ⟨29696 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_29824 : ∀ k : Fin 128, checkMirrorFields ⟨29824 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_29952 : ∀ k : Fin 128, checkMirrorFields ⟨29952 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_30080 : ∀ k : Fin 128, checkMirrorFields ⟨30080 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_30208 : ∀ k : Fin 128, checkMirrorFields ⟨30208 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_30336 : ∀ k : Fin 128, checkMirrorFields ⟨30336 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_30464 : ∀ k : Fin 128, checkMirrorFields ⟨30464 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_30592 : ∀ k : Fin 128, checkMirrorFields ⟨30592 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_30720 : ∀ k : Fin 128, checkMirrorFields ⟨30720 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_30848 : ∀ k : Fin 128, checkMirrorFields ⟨30848 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_30976 : ∀ k : Fin 128, checkMirrorFields ⟨30976 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_31104 : ∀ k : Fin 128, checkMirrorFields ⟨31104 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_31232 : ∀ k : Fin 128, checkMirrorFields ⟨31232 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_31360 : ∀ k : Fin 128, checkMirrorFields ⟨31360 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_31488 : ∀ k : Fin 128, checkMirrorFields ⟨31488 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_31616 : ∀ k : Fin 128, checkMirrorFields ⟨31616 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_31744 : ∀ k : Fin 128, checkMirrorFields ⟨31744 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_31872 : ∀ k : Fin 128, checkMirrorFields ⟨31872 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_32000 : ∀ k : Fin 128, checkMirrorFields ⟨32000 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_32128 : ∀ k : Fin 128, checkMirrorFields ⟨32128 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_32256 : ∀ k : Fin 128, checkMirrorFields ⟨32256 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_32384 : ∀ k : Fin 128, checkMirrorFields ⟨32384 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_32512 : ∀ k : Fin 128, checkMirrorFields ⟨32512 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_32640 : ∀ k : Fin 128, checkMirrorFields ⟨32640 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_32768 : ∀ k : Fin 128, checkMirrorFields ⟨32768 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_32896 : ∀ k : Fin 128, checkMirrorFields ⟨32896 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_33024 : ∀ k : Fin 128, checkMirrorFields ⟨33024 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_33152 : ∀ k : Fin 128, checkMirrorFields ⟨33152 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_33280 : ∀ k : Fin 128, checkMirrorFields ⟨33280 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_33408 : ∀ k : Fin 128, checkMirrorFields ⟨33408 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_33536 : ∀ k : Fin 128, checkMirrorFields ⟨33536 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_33664 : ∀ k : Fin 128, checkMirrorFields ⟨33664 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_33792 : ∀ k : Fin 128, checkMirrorFields ⟨33792 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_33920 : ∀ k : Fin 128, checkMirrorFields ⟨33920 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_34048 : ∀ k : Fin 128, checkMirrorFields ⟨34048 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_34176 : ∀ k : Fin 128, checkMirrorFields ⟨34176 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_34304 : ∀ k : Fin 128, checkMirrorFields ⟨34304 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_34432 : ∀ k : Fin 128, checkMirrorFields ⟨34432 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_34560 : ∀ k : Fin 128, checkMirrorFields ⟨34560 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_34688 : ∀ k : Fin 128, checkMirrorFields ⟨34688 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_34816 : ∀ k : Fin 128, checkMirrorFields ⟨34816 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_34944 : ∀ k : Fin 128, checkMirrorFields ⟨34944 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_35072 : ∀ k : Fin 128, checkMirrorFields ⟨35072 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_35200 : ∀ k : Fin 128, checkMirrorFields ⟨35200 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_35328 : ∀ k : Fin 128, checkMirrorFields ⟨35328 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_35456 : ∀ k : Fin 128, checkMirrorFields ⟨35456 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_35584 : ∀ k : Fin 128, checkMirrorFields ⟨35584 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_35712 : ∀ k : Fin 128, checkMirrorFields ⟨35712 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_35840 : ∀ k : Fin 128, checkMirrorFields ⟨35840 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_35968 : ∀ k : Fin 128, checkMirrorFields ⟨35968 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_36096 : ∀ k : Fin 128, checkMirrorFields ⟨36096 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_36224 : ∀ k : Fin 128, checkMirrorFields ⟨36224 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_36352 : ∀ k : Fin 128, checkMirrorFields ⟨36352 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_36480 : ∀ k : Fin 128, checkMirrorFields ⟨36480 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_36608 : ∀ k : Fin 128, checkMirrorFields ⟨36608 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_36736 : ∀ k : Fin 128, checkMirrorFields ⟨36736 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_36864 : ∀ k : Fin 128, checkMirrorFields ⟨36864 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_36992 : ∀ k : Fin 128, checkMirrorFields ⟨36992 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_37120 : ∀ k : Fin 128, checkMirrorFields ⟨37120 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_37248 : ∀ k : Fin 128, checkMirrorFields ⟨37248 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_37376 : ∀ k : Fin 128, checkMirrorFields ⟨37376 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_37504 : ∀ k : Fin 128, checkMirrorFields ⟨37504 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_37632 : ∀ k : Fin 128, checkMirrorFields ⟨37632 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_37760 : ∀ k : Fin 128, checkMirrorFields ⟨37760 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_37888 : ∀ k : Fin 128, checkMirrorFields ⟨37888 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_38016 : ∀ k : Fin 128, checkMirrorFields ⟨38016 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_38144 : ∀ k : Fin 128, checkMirrorFields ⟨38144 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_38272 : ∀ k : Fin 128, checkMirrorFields ⟨38272 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_38400 : ∀ k : Fin 128, checkMirrorFields ⟨38400 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_38528 : ∀ k : Fin 128, checkMirrorFields ⟨38528 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_38656 : ∀ k : Fin 128, checkMirrorFields ⟨38656 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_38784 : ∀ k : Fin 128, checkMirrorFields ⟨38784 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_38912 : ∀ k : Fin 128, checkMirrorFields ⟨38912 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_39040 : ∀ k : Fin 128, checkMirrorFields ⟨39040 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_39168 : ∀ k : Fin 128, checkMirrorFields ⟨39168 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_39296 : ∀ k : Fin 128, checkMirrorFields ⟨39296 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_39424 : ∀ k : Fin 128, checkMirrorFields ⟨39424 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_39552 : ∀ k : Fin 128, checkMirrorFields ⟨39552 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_39680 : ∀ k : Fin 128, checkMirrorFields ⟨39680 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_39808 : ∀ k : Fin 128, checkMirrorFields ⟨39808 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_39936 : ∀ k : Fin 128, checkMirrorFields ⟨39936 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_40064 : ∀ k : Fin 128, checkMirrorFields ⟨40064 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_40192 : ∀ k : Fin 128, checkMirrorFields ⟨40192 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_40320 : ∀ k : Fin 128, checkMirrorFields ⟨40320 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_40448 : ∀ k : Fin 128, checkMirrorFields ⟨40448 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_40576 : ∀ k : Fin 128, checkMirrorFields ⟨40576 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_40704 : ∀ k : Fin 128, checkMirrorFields ⟨40704 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_40832 : ∀ k : Fin 128, checkMirrorFields ⟨40832 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_40960 : ∀ k : Fin 128, checkMirrorFields ⟨40960 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_41088 : ∀ k : Fin 128, checkMirrorFields ⟨41088 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem mirror_41216 : ∀ k : Fin 78, checkMirrorFields ⟨41216 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
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
private theorem transitions_8576 : ∀ k : Fin 128, checkTransitions ⟨8576 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_8704 : ∀ k : Fin 128, checkTransitions ⟨8704 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_8832 : ∀ k : Fin 128, checkTransitions ⟨8832 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_8960 : ∀ k : Fin 128, checkTransitions ⟨8960 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_9088 : ∀ k : Fin 128, checkTransitions ⟨9088 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_9216 : ∀ k : Fin 128, checkTransitions ⟨9216 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_9344 : ∀ k : Fin 128, checkTransitions ⟨9344 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_9472 : ∀ k : Fin 128, checkTransitions ⟨9472 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_9600 : ∀ k : Fin 128, checkTransitions ⟨9600 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_9728 : ∀ k : Fin 128, checkTransitions ⟨9728 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_9856 : ∀ k : Fin 128, checkTransitions ⟨9856 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_9984 : ∀ k : Fin 128, checkTransitions ⟨9984 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_10112 : ∀ k : Fin 128, checkTransitions ⟨10112 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_10240 : ∀ k : Fin 128, checkTransitions ⟨10240 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_10368 : ∀ k : Fin 128, checkTransitions ⟨10368 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_10496 : ∀ k : Fin 128, checkTransitions ⟨10496 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_10624 : ∀ k : Fin 128, checkTransitions ⟨10624 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_10752 : ∀ k : Fin 128, checkTransitions ⟨10752 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_10880 : ∀ k : Fin 128, checkTransitions ⟨10880 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_11008 : ∀ k : Fin 128, checkTransitions ⟨11008 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_11136 : ∀ k : Fin 128, checkTransitions ⟨11136 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_11264 : ∀ k : Fin 128, checkTransitions ⟨11264 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_11392 : ∀ k : Fin 128, checkTransitions ⟨11392 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_11520 : ∀ k : Fin 128, checkTransitions ⟨11520 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_11648 : ∀ k : Fin 128, checkTransitions ⟨11648 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_11776 : ∀ k : Fin 128, checkTransitions ⟨11776 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_11904 : ∀ k : Fin 128, checkTransitions ⟨11904 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_12032 : ∀ k : Fin 128, checkTransitions ⟨12032 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_12160 : ∀ k : Fin 128, checkTransitions ⟨12160 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_12288 : ∀ k : Fin 128, checkTransitions ⟨12288 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_12416 : ∀ k : Fin 128, checkTransitions ⟨12416 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_12544 : ∀ k : Fin 128, checkTransitions ⟨12544 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_12672 : ∀ k : Fin 128, checkTransitions ⟨12672 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_12800 : ∀ k : Fin 128, checkTransitions ⟨12800 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_12928 : ∀ k : Fin 128, checkTransitions ⟨12928 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_13056 : ∀ k : Fin 128, checkTransitions ⟨13056 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_13184 : ∀ k : Fin 128, checkTransitions ⟨13184 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_13312 : ∀ k : Fin 128, checkTransitions ⟨13312 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_13440 : ∀ k : Fin 128, checkTransitions ⟨13440 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_13568 : ∀ k : Fin 128, checkTransitions ⟨13568 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_13696 : ∀ k : Fin 128, checkTransitions ⟨13696 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_13824 : ∀ k : Fin 128, checkTransitions ⟨13824 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_13952 : ∀ k : Fin 128, checkTransitions ⟨13952 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_14080 : ∀ k : Fin 128, checkTransitions ⟨14080 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_14208 : ∀ k : Fin 128, checkTransitions ⟨14208 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_14336 : ∀ k : Fin 128, checkTransitions ⟨14336 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_14464 : ∀ k : Fin 128, checkTransitions ⟨14464 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_14592 : ∀ k : Fin 128, checkTransitions ⟨14592 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_14720 : ∀ k : Fin 128, checkTransitions ⟨14720 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_14848 : ∀ k : Fin 128, checkTransitions ⟨14848 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_14976 : ∀ k : Fin 128, checkTransitions ⟨14976 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_15104 : ∀ k : Fin 128, checkTransitions ⟨15104 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_15232 : ∀ k : Fin 128, checkTransitions ⟨15232 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_15360 : ∀ k : Fin 128, checkTransitions ⟨15360 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_15488 : ∀ k : Fin 128, checkTransitions ⟨15488 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_15616 : ∀ k : Fin 128, checkTransitions ⟨15616 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_15744 : ∀ k : Fin 128, checkTransitions ⟨15744 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_15872 : ∀ k : Fin 128, checkTransitions ⟨15872 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_16000 : ∀ k : Fin 128, checkTransitions ⟨16000 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_16128 : ∀ k : Fin 128, checkTransitions ⟨16128 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_16256 : ∀ k : Fin 128, checkTransitions ⟨16256 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_16384 : ∀ k : Fin 128, checkTransitions ⟨16384 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_16512 : ∀ k : Fin 128, checkTransitions ⟨16512 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_16640 : ∀ k : Fin 128, checkTransitions ⟨16640 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_16768 : ∀ k : Fin 128, checkTransitions ⟨16768 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_16896 : ∀ k : Fin 128, checkTransitions ⟨16896 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_17024 : ∀ k : Fin 128, checkTransitions ⟨17024 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_17152 : ∀ k : Fin 128, checkTransitions ⟨17152 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_17280 : ∀ k : Fin 128, checkTransitions ⟨17280 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_17408 : ∀ k : Fin 128, checkTransitions ⟨17408 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_17536 : ∀ k : Fin 128, checkTransitions ⟨17536 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_17664 : ∀ k : Fin 128, checkTransitions ⟨17664 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_17792 : ∀ k : Fin 128, checkTransitions ⟨17792 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_17920 : ∀ k : Fin 128, checkTransitions ⟨17920 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_18048 : ∀ k : Fin 128, checkTransitions ⟨18048 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_18176 : ∀ k : Fin 128, checkTransitions ⟨18176 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_18304 : ∀ k : Fin 128, checkTransitions ⟨18304 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_18432 : ∀ k : Fin 128, checkTransitions ⟨18432 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_18560 : ∀ k : Fin 128, checkTransitions ⟨18560 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_18688 : ∀ k : Fin 128, checkTransitions ⟨18688 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_18816 : ∀ k : Fin 128, checkTransitions ⟨18816 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_18944 : ∀ k : Fin 128, checkTransitions ⟨18944 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_19072 : ∀ k : Fin 128, checkTransitions ⟨19072 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_19200 : ∀ k : Fin 128, checkTransitions ⟨19200 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_19328 : ∀ k : Fin 128, checkTransitions ⟨19328 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_19456 : ∀ k : Fin 128, checkTransitions ⟨19456 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_19584 : ∀ k : Fin 128, checkTransitions ⟨19584 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_19712 : ∀ k : Fin 128, checkTransitions ⟨19712 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_19840 : ∀ k : Fin 128, checkTransitions ⟨19840 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_19968 : ∀ k : Fin 128, checkTransitions ⟨19968 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_20096 : ∀ k : Fin 128, checkTransitions ⟨20096 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_20224 : ∀ k : Fin 128, checkTransitions ⟨20224 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_20352 : ∀ k : Fin 128, checkTransitions ⟨20352 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_20480 : ∀ k : Fin 128, checkTransitions ⟨20480 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_20608 : ∀ k : Fin 128, checkTransitions ⟨20608 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_20736 : ∀ k : Fin 128, checkTransitions ⟨20736 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_20864 : ∀ k : Fin 128, checkTransitions ⟨20864 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_20992 : ∀ k : Fin 128, checkTransitions ⟨20992 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_21120 : ∀ k : Fin 128, checkTransitions ⟨21120 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_21248 : ∀ k : Fin 128, checkTransitions ⟨21248 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_21376 : ∀ k : Fin 128, checkTransitions ⟨21376 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_21504 : ∀ k : Fin 128, checkTransitions ⟨21504 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_21632 : ∀ k : Fin 128, checkTransitions ⟨21632 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_21760 : ∀ k : Fin 128, checkTransitions ⟨21760 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_21888 : ∀ k : Fin 128, checkTransitions ⟨21888 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_22016 : ∀ k : Fin 128, checkTransitions ⟨22016 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_22144 : ∀ k : Fin 128, checkTransitions ⟨22144 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_22272 : ∀ k : Fin 128, checkTransitions ⟨22272 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_22400 : ∀ k : Fin 128, checkTransitions ⟨22400 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_22528 : ∀ k : Fin 128, checkTransitions ⟨22528 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_22656 : ∀ k : Fin 128, checkTransitions ⟨22656 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_22784 : ∀ k : Fin 128, checkTransitions ⟨22784 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_22912 : ∀ k : Fin 128, checkTransitions ⟨22912 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_23040 : ∀ k : Fin 128, checkTransitions ⟨23040 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_23168 : ∀ k : Fin 128, checkTransitions ⟨23168 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_23296 : ∀ k : Fin 128, checkTransitions ⟨23296 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_23424 : ∀ k : Fin 128, checkTransitions ⟨23424 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_23552 : ∀ k : Fin 128, checkTransitions ⟨23552 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_23680 : ∀ k : Fin 128, checkTransitions ⟨23680 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_23808 : ∀ k : Fin 128, checkTransitions ⟨23808 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_23936 : ∀ k : Fin 128, checkTransitions ⟨23936 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_24064 : ∀ k : Fin 128, checkTransitions ⟨24064 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_24192 : ∀ k : Fin 128, checkTransitions ⟨24192 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_24320 : ∀ k : Fin 128, checkTransitions ⟨24320 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_24448 : ∀ k : Fin 128, checkTransitions ⟨24448 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_24576 : ∀ k : Fin 128, checkTransitions ⟨24576 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_24704 : ∀ k : Fin 128, checkTransitions ⟨24704 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_24832 : ∀ k : Fin 128, checkTransitions ⟨24832 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_24960 : ∀ k : Fin 128, checkTransitions ⟨24960 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_25088 : ∀ k : Fin 128, checkTransitions ⟨25088 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_25216 : ∀ k : Fin 128, checkTransitions ⟨25216 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_25344 : ∀ k : Fin 128, checkTransitions ⟨25344 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_25472 : ∀ k : Fin 128, checkTransitions ⟨25472 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_25600 : ∀ k : Fin 128, checkTransitions ⟨25600 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_25728 : ∀ k : Fin 128, checkTransitions ⟨25728 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_25856 : ∀ k : Fin 128, checkTransitions ⟨25856 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_25984 : ∀ k : Fin 128, checkTransitions ⟨25984 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_26112 : ∀ k : Fin 128, checkTransitions ⟨26112 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_26240 : ∀ k : Fin 128, checkTransitions ⟨26240 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_26368 : ∀ k : Fin 128, checkTransitions ⟨26368 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_26496 : ∀ k : Fin 128, checkTransitions ⟨26496 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_26624 : ∀ k : Fin 128, checkTransitions ⟨26624 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_26752 : ∀ k : Fin 128, checkTransitions ⟨26752 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_26880 : ∀ k : Fin 128, checkTransitions ⟨26880 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_27008 : ∀ k : Fin 128, checkTransitions ⟨27008 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_27136 : ∀ k : Fin 128, checkTransitions ⟨27136 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_27264 : ∀ k : Fin 128, checkTransitions ⟨27264 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_27392 : ∀ k : Fin 128, checkTransitions ⟨27392 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_27520 : ∀ k : Fin 128, checkTransitions ⟨27520 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_27648 : ∀ k : Fin 128, checkTransitions ⟨27648 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_27776 : ∀ k : Fin 128, checkTransitions ⟨27776 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_27904 : ∀ k : Fin 128, checkTransitions ⟨27904 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_28032 : ∀ k : Fin 128, checkTransitions ⟨28032 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_28160 : ∀ k : Fin 128, checkTransitions ⟨28160 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_28288 : ∀ k : Fin 128, checkTransitions ⟨28288 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_28416 : ∀ k : Fin 128, checkTransitions ⟨28416 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_28544 : ∀ k : Fin 128, checkTransitions ⟨28544 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_28672 : ∀ k : Fin 128, checkTransitions ⟨28672 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_28800 : ∀ k : Fin 128, checkTransitions ⟨28800 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_28928 : ∀ k : Fin 128, checkTransitions ⟨28928 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_29056 : ∀ k : Fin 128, checkTransitions ⟨29056 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_29184 : ∀ k : Fin 128, checkTransitions ⟨29184 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_29312 : ∀ k : Fin 128, checkTransitions ⟨29312 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_29440 : ∀ k : Fin 128, checkTransitions ⟨29440 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_29568 : ∀ k : Fin 128, checkTransitions ⟨29568 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_29696 : ∀ k : Fin 128, checkTransitions ⟨29696 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_29824 : ∀ k : Fin 128, checkTransitions ⟨29824 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_29952 : ∀ k : Fin 128, checkTransitions ⟨29952 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_30080 : ∀ k : Fin 128, checkTransitions ⟨30080 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_30208 : ∀ k : Fin 128, checkTransitions ⟨30208 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_30336 : ∀ k : Fin 128, checkTransitions ⟨30336 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_30464 : ∀ k : Fin 128, checkTransitions ⟨30464 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_30592 : ∀ k : Fin 128, checkTransitions ⟨30592 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_30720 : ∀ k : Fin 128, checkTransitions ⟨30720 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_30848 : ∀ k : Fin 128, checkTransitions ⟨30848 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_30976 : ∀ k : Fin 128, checkTransitions ⟨30976 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_31104 : ∀ k : Fin 128, checkTransitions ⟨31104 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_31232 : ∀ k : Fin 128, checkTransitions ⟨31232 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_31360 : ∀ k : Fin 128, checkTransitions ⟨31360 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_31488 : ∀ k : Fin 128, checkTransitions ⟨31488 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_31616 : ∀ k : Fin 128, checkTransitions ⟨31616 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_31744 : ∀ k : Fin 128, checkTransitions ⟨31744 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_31872 : ∀ k : Fin 128, checkTransitions ⟨31872 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_32000 : ∀ k : Fin 128, checkTransitions ⟨32000 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_32128 : ∀ k : Fin 128, checkTransitions ⟨32128 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_32256 : ∀ k : Fin 128, checkTransitions ⟨32256 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_32384 : ∀ k : Fin 128, checkTransitions ⟨32384 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_32512 : ∀ k : Fin 128, checkTransitions ⟨32512 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_32640 : ∀ k : Fin 128, checkTransitions ⟨32640 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_32768 : ∀ k : Fin 128, checkTransitions ⟨32768 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_32896 : ∀ k : Fin 128, checkTransitions ⟨32896 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_33024 : ∀ k : Fin 128, checkTransitions ⟨33024 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_33152 : ∀ k : Fin 128, checkTransitions ⟨33152 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_33280 : ∀ k : Fin 128, checkTransitions ⟨33280 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_33408 : ∀ k : Fin 128, checkTransitions ⟨33408 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_33536 : ∀ k : Fin 128, checkTransitions ⟨33536 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_33664 : ∀ k : Fin 128, checkTransitions ⟨33664 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_33792 : ∀ k : Fin 128, checkTransitions ⟨33792 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_33920 : ∀ k : Fin 128, checkTransitions ⟨33920 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_34048 : ∀ k : Fin 128, checkTransitions ⟨34048 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_34176 : ∀ k : Fin 128, checkTransitions ⟨34176 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_34304 : ∀ k : Fin 128, checkTransitions ⟨34304 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_34432 : ∀ k : Fin 128, checkTransitions ⟨34432 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_34560 : ∀ k : Fin 128, checkTransitions ⟨34560 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_34688 : ∀ k : Fin 128, checkTransitions ⟨34688 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_34816 : ∀ k : Fin 128, checkTransitions ⟨34816 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_34944 : ∀ k : Fin 128, checkTransitions ⟨34944 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_35072 : ∀ k : Fin 128, checkTransitions ⟨35072 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_35200 : ∀ k : Fin 128, checkTransitions ⟨35200 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_35328 : ∀ k : Fin 128, checkTransitions ⟨35328 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_35456 : ∀ k : Fin 128, checkTransitions ⟨35456 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_35584 : ∀ k : Fin 128, checkTransitions ⟨35584 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_35712 : ∀ k : Fin 128, checkTransitions ⟨35712 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_35840 : ∀ k : Fin 128, checkTransitions ⟨35840 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_35968 : ∀ k : Fin 128, checkTransitions ⟨35968 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_36096 : ∀ k : Fin 128, checkTransitions ⟨36096 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_36224 : ∀ k : Fin 128, checkTransitions ⟨36224 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_36352 : ∀ k : Fin 128, checkTransitions ⟨36352 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_36480 : ∀ k : Fin 128, checkTransitions ⟨36480 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_36608 : ∀ k : Fin 128, checkTransitions ⟨36608 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_36736 : ∀ k : Fin 128, checkTransitions ⟨36736 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_36864 : ∀ k : Fin 128, checkTransitions ⟨36864 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_36992 : ∀ k : Fin 128, checkTransitions ⟨36992 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_37120 : ∀ k : Fin 128, checkTransitions ⟨37120 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_37248 : ∀ k : Fin 128, checkTransitions ⟨37248 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_37376 : ∀ k : Fin 128, checkTransitions ⟨37376 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_37504 : ∀ k : Fin 128, checkTransitions ⟨37504 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_37632 : ∀ k : Fin 128, checkTransitions ⟨37632 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_37760 : ∀ k : Fin 128, checkTransitions ⟨37760 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_37888 : ∀ k : Fin 128, checkTransitions ⟨37888 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_38016 : ∀ k : Fin 128, checkTransitions ⟨38016 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_38144 : ∀ k : Fin 128, checkTransitions ⟨38144 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_38272 : ∀ k : Fin 128, checkTransitions ⟨38272 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_38400 : ∀ k : Fin 128, checkTransitions ⟨38400 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_38528 : ∀ k : Fin 128, checkTransitions ⟨38528 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_38656 : ∀ k : Fin 128, checkTransitions ⟨38656 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_38784 : ∀ k : Fin 128, checkTransitions ⟨38784 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_38912 : ∀ k : Fin 128, checkTransitions ⟨38912 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_39040 : ∀ k : Fin 128, checkTransitions ⟨39040 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_39168 : ∀ k : Fin 128, checkTransitions ⟨39168 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_39296 : ∀ k : Fin 128, checkTransitions ⟨39296 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_39424 : ∀ k : Fin 128, checkTransitions ⟨39424 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_39552 : ∀ k : Fin 128, checkTransitions ⟨39552 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_39680 : ∀ k : Fin 128, checkTransitions ⟨39680 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_39808 : ∀ k : Fin 128, checkTransitions ⟨39808 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_39936 : ∀ k : Fin 128, checkTransitions ⟨39936 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_40064 : ∀ k : Fin 128, checkTransitions ⟨40064 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_40192 : ∀ k : Fin 128, checkTransitions ⟨40192 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_40320 : ∀ k : Fin 128, checkTransitions ⟨40320 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_40448 : ∀ k : Fin 128, checkTransitions ⟨40448 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_40576 : ∀ k : Fin 128, checkTransitions ⟨40576 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_40704 : ∀ k : Fin 128, checkTransitions ⟨40704 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_40832 : ∀ k : Fin 128, checkTransitions ⟨40832 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_40960 : ∀ k : Fin 128, checkTransitions ⟨40960 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_41088 : ∀ k : Fin 128, checkTransitions ⟨41088 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel
private theorem transitions_41216 : ∀ k : Fin 78, checkTransitions ⟨41216 + k.val, by have h := k.isLt; omega⟩ = true := by decide +kernel

private def chunkLength (b : Fin 323) : ℕ :=
  min 128 (41294 - 128 * b.val)

private def chunkIndex (b : Fin 323) (k : Fin (chunkLength b)) : Fin 41294 :=
  ⟨128 * b.val + k.val, by
    have hb := b.isLt
    have hk := k.isLt
    dsimp [chunkLength] at hk
    omega⟩

private theorem mirror_blocks : ∀ b : Fin 323, ∀ k : Fin (chunkLength b),
    checkMirrorFields (chunkIndex b k) = true :=
  Fin.cases mirror_0 (Fin.cases mirror_128 (Fin.cases mirror_256 (Fin.cases mirror_384 (Fin.cases mirror_512 (Fin.cases mirror_640 (Fin.cases mirror_768 (Fin.cases mirror_896 (Fin.cases mirror_1024 (Fin.cases mirror_1152 (Fin.cases mirror_1280 (Fin.cases mirror_1408 (Fin.cases mirror_1536 (Fin.cases mirror_1664 (Fin.cases mirror_1792 (Fin.cases mirror_1920 (
  Fin.cases mirror_2048 (Fin.cases mirror_2176 (Fin.cases mirror_2304 (Fin.cases mirror_2432 (Fin.cases mirror_2560 (Fin.cases mirror_2688 (Fin.cases mirror_2816 (Fin.cases mirror_2944 (Fin.cases mirror_3072 (Fin.cases mirror_3200 (Fin.cases mirror_3328 (Fin.cases mirror_3456 (Fin.cases mirror_3584 (Fin.cases mirror_3712 (Fin.cases mirror_3840 (Fin.cases mirror_3968 (
  Fin.cases mirror_4096 (Fin.cases mirror_4224 (Fin.cases mirror_4352 (Fin.cases mirror_4480 (Fin.cases mirror_4608 (Fin.cases mirror_4736 (Fin.cases mirror_4864 (Fin.cases mirror_4992 (Fin.cases mirror_5120 (Fin.cases mirror_5248 (Fin.cases mirror_5376 (Fin.cases mirror_5504 (Fin.cases mirror_5632 (Fin.cases mirror_5760 (Fin.cases mirror_5888 (Fin.cases mirror_6016 (
  Fin.cases mirror_6144 (Fin.cases mirror_6272 (Fin.cases mirror_6400 (Fin.cases mirror_6528 (Fin.cases mirror_6656 (Fin.cases mirror_6784 (Fin.cases mirror_6912 (Fin.cases mirror_7040 (Fin.cases mirror_7168 (Fin.cases mirror_7296 (Fin.cases mirror_7424 (Fin.cases mirror_7552 (Fin.cases mirror_7680 (Fin.cases mirror_7808 (Fin.cases mirror_7936 (Fin.cases mirror_8064 (
  Fin.cases mirror_8192 (Fin.cases mirror_8320 (Fin.cases mirror_8448 (Fin.cases mirror_8576 (Fin.cases mirror_8704 (Fin.cases mirror_8832 (Fin.cases mirror_8960 (Fin.cases mirror_9088 (Fin.cases mirror_9216 (Fin.cases mirror_9344 (Fin.cases mirror_9472 (Fin.cases mirror_9600 (Fin.cases mirror_9728 (Fin.cases mirror_9856 (Fin.cases mirror_9984 (Fin.cases mirror_10112 (
  Fin.cases mirror_10240 (Fin.cases mirror_10368 (Fin.cases mirror_10496 (Fin.cases mirror_10624 (Fin.cases mirror_10752 (Fin.cases mirror_10880 (Fin.cases mirror_11008 (Fin.cases mirror_11136 (Fin.cases mirror_11264 (Fin.cases mirror_11392 (Fin.cases mirror_11520 (Fin.cases mirror_11648 (Fin.cases mirror_11776 (Fin.cases mirror_11904 (Fin.cases mirror_12032 (Fin.cases mirror_12160 (
  Fin.cases mirror_12288 (Fin.cases mirror_12416 (Fin.cases mirror_12544 (Fin.cases mirror_12672 (Fin.cases mirror_12800 (Fin.cases mirror_12928 (Fin.cases mirror_13056 (Fin.cases mirror_13184 (Fin.cases mirror_13312 (Fin.cases mirror_13440 (Fin.cases mirror_13568 (Fin.cases mirror_13696 (Fin.cases mirror_13824 (Fin.cases mirror_13952 (Fin.cases mirror_14080 (Fin.cases mirror_14208 (
  Fin.cases mirror_14336 (Fin.cases mirror_14464 (Fin.cases mirror_14592 (Fin.cases mirror_14720 (Fin.cases mirror_14848 (Fin.cases mirror_14976 (Fin.cases mirror_15104 (Fin.cases mirror_15232 (Fin.cases mirror_15360 (Fin.cases mirror_15488 (Fin.cases mirror_15616 (Fin.cases mirror_15744 (Fin.cases mirror_15872 (Fin.cases mirror_16000 (Fin.cases mirror_16128 (Fin.cases mirror_16256 (
  Fin.cases mirror_16384 (Fin.cases mirror_16512 (Fin.cases mirror_16640 (Fin.cases mirror_16768 (Fin.cases mirror_16896 (Fin.cases mirror_17024 (Fin.cases mirror_17152 (Fin.cases mirror_17280 (Fin.cases mirror_17408 (Fin.cases mirror_17536 (Fin.cases mirror_17664 (Fin.cases mirror_17792 (Fin.cases mirror_17920 (Fin.cases mirror_18048 (Fin.cases mirror_18176 (Fin.cases mirror_18304 (
  Fin.cases mirror_18432 (Fin.cases mirror_18560 (Fin.cases mirror_18688 (Fin.cases mirror_18816 (Fin.cases mirror_18944 (Fin.cases mirror_19072 (Fin.cases mirror_19200 (Fin.cases mirror_19328 (Fin.cases mirror_19456 (Fin.cases mirror_19584 (Fin.cases mirror_19712 (Fin.cases mirror_19840 (Fin.cases mirror_19968 (Fin.cases mirror_20096 (Fin.cases mirror_20224 (Fin.cases mirror_20352 (
  Fin.cases mirror_20480 (Fin.cases mirror_20608 (Fin.cases mirror_20736 (Fin.cases mirror_20864 (Fin.cases mirror_20992 (Fin.cases mirror_21120 (Fin.cases mirror_21248 (Fin.cases mirror_21376 (Fin.cases mirror_21504 (Fin.cases mirror_21632 (Fin.cases mirror_21760 (Fin.cases mirror_21888 (Fin.cases mirror_22016 (Fin.cases mirror_22144 (Fin.cases mirror_22272 (Fin.cases mirror_22400 (
  Fin.cases mirror_22528 (Fin.cases mirror_22656 (Fin.cases mirror_22784 (Fin.cases mirror_22912 (Fin.cases mirror_23040 (Fin.cases mirror_23168 (Fin.cases mirror_23296 (Fin.cases mirror_23424 (Fin.cases mirror_23552 (Fin.cases mirror_23680 (Fin.cases mirror_23808 (Fin.cases mirror_23936 (Fin.cases mirror_24064 (Fin.cases mirror_24192 (Fin.cases mirror_24320 (Fin.cases mirror_24448 (
  Fin.cases mirror_24576 (Fin.cases mirror_24704 (Fin.cases mirror_24832 (Fin.cases mirror_24960 (Fin.cases mirror_25088 (Fin.cases mirror_25216 (Fin.cases mirror_25344 (Fin.cases mirror_25472 (Fin.cases mirror_25600 (Fin.cases mirror_25728 (Fin.cases mirror_25856 (Fin.cases mirror_25984 (Fin.cases mirror_26112 (Fin.cases mirror_26240 (Fin.cases mirror_26368 (Fin.cases mirror_26496 (
  Fin.cases mirror_26624 (Fin.cases mirror_26752 (Fin.cases mirror_26880 (Fin.cases mirror_27008 (Fin.cases mirror_27136 (Fin.cases mirror_27264 (Fin.cases mirror_27392 (Fin.cases mirror_27520 (Fin.cases mirror_27648 (Fin.cases mirror_27776 (Fin.cases mirror_27904 (Fin.cases mirror_28032 (Fin.cases mirror_28160 (Fin.cases mirror_28288 (Fin.cases mirror_28416 (Fin.cases mirror_28544 (
  Fin.cases mirror_28672 (Fin.cases mirror_28800 (Fin.cases mirror_28928 (Fin.cases mirror_29056 (Fin.cases mirror_29184 (Fin.cases mirror_29312 (Fin.cases mirror_29440 (Fin.cases mirror_29568 (Fin.cases mirror_29696 (Fin.cases mirror_29824 (Fin.cases mirror_29952 (Fin.cases mirror_30080 (Fin.cases mirror_30208 (Fin.cases mirror_30336 (Fin.cases mirror_30464 (Fin.cases mirror_30592 (
  Fin.cases mirror_30720 (Fin.cases mirror_30848 (Fin.cases mirror_30976 (Fin.cases mirror_31104 (Fin.cases mirror_31232 (Fin.cases mirror_31360 (Fin.cases mirror_31488 (Fin.cases mirror_31616 (Fin.cases mirror_31744 (Fin.cases mirror_31872 (Fin.cases mirror_32000 (Fin.cases mirror_32128 (Fin.cases mirror_32256 (Fin.cases mirror_32384 (Fin.cases mirror_32512 (Fin.cases mirror_32640 (
  Fin.cases mirror_32768 (Fin.cases mirror_32896 (Fin.cases mirror_33024 (Fin.cases mirror_33152 (Fin.cases mirror_33280 (Fin.cases mirror_33408 (Fin.cases mirror_33536 (Fin.cases mirror_33664 (Fin.cases mirror_33792 (Fin.cases mirror_33920 (Fin.cases mirror_34048 (Fin.cases mirror_34176 (Fin.cases mirror_34304 (Fin.cases mirror_34432 (Fin.cases mirror_34560 (Fin.cases mirror_34688 (
  Fin.cases mirror_34816 (Fin.cases mirror_34944 (Fin.cases mirror_35072 (Fin.cases mirror_35200 (Fin.cases mirror_35328 (Fin.cases mirror_35456 (Fin.cases mirror_35584 (Fin.cases mirror_35712 (Fin.cases mirror_35840 (Fin.cases mirror_35968 (Fin.cases mirror_36096 (Fin.cases mirror_36224 (Fin.cases mirror_36352 (Fin.cases mirror_36480 (Fin.cases mirror_36608 (Fin.cases mirror_36736 (
  Fin.cases mirror_36864 (Fin.cases mirror_36992 (Fin.cases mirror_37120 (Fin.cases mirror_37248 (Fin.cases mirror_37376 (Fin.cases mirror_37504 (Fin.cases mirror_37632 (Fin.cases mirror_37760 (Fin.cases mirror_37888 (Fin.cases mirror_38016 (Fin.cases mirror_38144 (Fin.cases mirror_38272 (Fin.cases mirror_38400 (Fin.cases mirror_38528 (Fin.cases mirror_38656 (Fin.cases mirror_38784 (
  Fin.cases mirror_38912 (Fin.cases mirror_39040 (Fin.cases mirror_39168 (Fin.cases mirror_39296 (Fin.cases mirror_39424 (Fin.cases mirror_39552 (Fin.cases mirror_39680 (Fin.cases mirror_39808 (Fin.cases mirror_39936 (Fin.cases mirror_40064 (Fin.cases mirror_40192 (Fin.cases mirror_40320 (Fin.cases mirror_40448 (Fin.cases mirror_40576 (Fin.cases mirror_40704 (Fin.cases mirror_40832 (
  Fin.cases mirror_40960 (Fin.cases mirror_41088 (Fin.cases mirror_41216 (
  (fun b : Fin 0 => Fin.elim0 b))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

private theorem transitions_blocks : ∀ b : Fin 323, ∀ k : Fin (chunkLength b),
    checkTransitions (chunkIndex b k) = true :=
  Fin.cases transitions_0 (Fin.cases transitions_128 (Fin.cases transitions_256 (Fin.cases transitions_384 (Fin.cases transitions_512 (Fin.cases transitions_640 (Fin.cases transitions_768 (Fin.cases transitions_896 (Fin.cases transitions_1024 (Fin.cases transitions_1152 (Fin.cases transitions_1280 (Fin.cases transitions_1408 (Fin.cases transitions_1536 (Fin.cases transitions_1664 (Fin.cases transitions_1792 (Fin.cases transitions_1920 (
  Fin.cases transitions_2048 (Fin.cases transitions_2176 (Fin.cases transitions_2304 (Fin.cases transitions_2432 (Fin.cases transitions_2560 (Fin.cases transitions_2688 (Fin.cases transitions_2816 (Fin.cases transitions_2944 (Fin.cases transitions_3072 (Fin.cases transitions_3200 (Fin.cases transitions_3328 (Fin.cases transitions_3456 (Fin.cases transitions_3584 (Fin.cases transitions_3712 (Fin.cases transitions_3840 (Fin.cases transitions_3968 (
  Fin.cases transitions_4096 (Fin.cases transitions_4224 (Fin.cases transitions_4352 (Fin.cases transitions_4480 (Fin.cases transitions_4608 (Fin.cases transitions_4736 (Fin.cases transitions_4864 (Fin.cases transitions_4992 (Fin.cases transitions_5120 (Fin.cases transitions_5248 (Fin.cases transitions_5376 (Fin.cases transitions_5504 (Fin.cases transitions_5632 (Fin.cases transitions_5760 (Fin.cases transitions_5888 (Fin.cases transitions_6016 (
  Fin.cases transitions_6144 (Fin.cases transitions_6272 (Fin.cases transitions_6400 (Fin.cases transitions_6528 (Fin.cases transitions_6656 (Fin.cases transitions_6784 (Fin.cases transitions_6912 (Fin.cases transitions_7040 (Fin.cases transitions_7168 (Fin.cases transitions_7296 (Fin.cases transitions_7424 (Fin.cases transitions_7552 (Fin.cases transitions_7680 (Fin.cases transitions_7808 (Fin.cases transitions_7936 (Fin.cases transitions_8064 (
  Fin.cases transitions_8192 (Fin.cases transitions_8320 (Fin.cases transitions_8448 (Fin.cases transitions_8576 (Fin.cases transitions_8704 (Fin.cases transitions_8832 (Fin.cases transitions_8960 (Fin.cases transitions_9088 (Fin.cases transitions_9216 (Fin.cases transitions_9344 (Fin.cases transitions_9472 (Fin.cases transitions_9600 (Fin.cases transitions_9728 (Fin.cases transitions_9856 (Fin.cases transitions_9984 (Fin.cases transitions_10112 (
  Fin.cases transitions_10240 (Fin.cases transitions_10368 (Fin.cases transitions_10496 (Fin.cases transitions_10624 (Fin.cases transitions_10752 (Fin.cases transitions_10880 (Fin.cases transitions_11008 (Fin.cases transitions_11136 (Fin.cases transitions_11264 (Fin.cases transitions_11392 (Fin.cases transitions_11520 (Fin.cases transitions_11648 (Fin.cases transitions_11776 (Fin.cases transitions_11904 (Fin.cases transitions_12032 (Fin.cases transitions_12160 (
  Fin.cases transitions_12288 (Fin.cases transitions_12416 (Fin.cases transitions_12544 (Fin.cases transitions_12672 (Fin.cases transitions_12800 (Fin.cases transitions_12928 (Fin.cases transitions_13056 (Fin.cases transitions_13184 (Fin.cases transitions_13312 (Fin.cases transitions_13440 (Fin.cases transitions_13568 (Fin.cases transitions_13696 (Fin.cases transitions_13824 (Fin.cases transitions_13952 (Fin.cases transitions_14080 (Fin.cases transitions_14208 (
  Fin.cases transitions_14336 (Fin.cases transitions_14464 (Fin.cases transitions_14592 (Fin.cases transitions_14720 (Fin.cases transitions_14848 (Fin.cases transitions_14976 (Fin.cases transitions_15104 (Fin.cases transitions_15232 (Fin.cases transitions_15360 (Fin.cases transitions_15488 (Fin.cases transitions_15616 (Fin.cases transitions_15744 (Fin.cases transitions_15872 (Fin.cases transitions_16000 (Fin.cases transitions_16128 (Fin.cases transitions_16256 (
  Fin.cases transitions_16384 (Fin.cases transitions_16512 (Fin.cases transitions_16640 (Fin.cases transitions_16768 (Fin.cases transitions_16896 (Fin.cases transitions_17024 (Fin.cases transitions_17152 (Fin.cases transitions_17280 (Fin.cases transitions_17408 (Fin.cases transitions_17536 (Fin.cases transitions_17664 (Fin.cases transitions_17792 (Fin.cases transitions_17920 (Fin.cases transitions_18048 (Fin.cases transitions_18176 (Fin.cases transitions_18304 (
  Fin.cases transitions_18432 (Fin.cases transitions_18560 (Fin.cases transitions_18688 (Fin.cases transitions_18816 (Fin.cases transitions_18944 (Fin.cases transitions_19072 (Fin.cases transitions_19200 (Fin.cases transitions_19328 (Fin.cases transitions_19456 (Fin.cases transitions_19584 (Fin.cases transitions_19712 (Fin.cases transitions_19840 (Fin.cases transitions_19968 (Fin.cases transitions_20096 (Fin.cases transitions_20224 (Fin.cases transitions_20352 (
  Fin.cases transitions_20480 (Fin.cases transitions_20608 (Fin.cases transitions_20736 (Fin.cases transitions_20864 (Fin.cases transitions_20992 (Fin.cases transitions_21120 (Fin.cases transitions_21248 (Fin.cases transitions_21376 (Fin.cases transitions_21504 (Fin.cases transitions_21632 (Fin.cases transitions_21760 (Fin.cases transitions_21888 (Fin.cases transitions_22016 (Fin.cases transitions_22144 (Fin.cases transitions_22272 (Fin.cases transitions_22400 (
  Fin.cases transitions_22528 (Fin.cases transitions_22656 (Fin.cases transitions_22784 (Fin.cases transitions_22912 (Fin.cases transitions_23040 (Fin.cases transitions_23168 (Fin.cases transitions_23296 (Fin.cases transitions_23424 (Fin.cases transitions_23552 (Fin.cases transitions_23680 (Fin.cases transitions_23808 (Fin.cases transitions_23936 (Fin.cases transitions_24064 (Fin.cases transitions_24192 (Fin.cases transitions_24320 (Fin.cases transitions_24448 (
  Fin.cases transitions_24576 (Fin.cases transitions_24704 (Fin.cases transitions_24832 (Fin.cases transitions_24960 (Fin.cases transitions_25088 (Fin.cases transitions_25216 (Fin.cases transitions_25344 (Fin.cases transitions_25472 (Fin.cases transitions_25600 (Fin.cases transitions_25728 (Fin.cases transitions_25856 (Fin.cases transitions_25984 (Fin.cases transitions_26112 (Fin.cases transitions_26240 (Fin.cases transitions_26368 (Fin.cases transitions_26496 (
  Fin.cases transitions_26624 (Fin.cases transitions_26752 (Fin.cases transitions_26880 (Fin.cases transitions_27008 (Fin.cases transitions_27136 (Fin.cases transitions_27264 (Fin.cases transitions_27392 (Fin.cases transitions_27520 (Fin.cases transitions_27648 (Fin.cases transitions_27776 (Fin.cases transitions_27904 (Fin.cases transitions_28032 (Fin.cases transitions_28160 (Fin.cases transitions_28288 (Fin.cases transitions_28416 (Fin.cases transitions_28544 (
  Fin.cases transitions_28672 (Fin.cases transitions_28800 (Fin.cases transitions_28928 (Fin.cases transitions_29056 (Fin.cases transitions_29184 (Fin.cases transitions_29312 (Fin.cases transitions_29440 (Fin.cases transitions_29568 (Fin.cases transitions_29696 (Fin.cases transitions_29824 (Fin.cases transitions_29952 (Fin.cases transitions_30080 (Fin.cases transitions_30208 (Fin.cases transitions_30336 (Fin.cases transitions_30464 (Fin.cases transitions_30592 (
  Fin.cases transitions_30720 (Fin.cases transitions_30848 (Fin.cases transitions_30976 (Fin.cases transitions_31104 (Fin.cases transitions_31232 (Fin.cases transitions_31360 (Fin.cases transitions_31488 (Fin.cases transitions_31616 (Fin.cases transitions_31744 (Fin.cases transitions_31872 (Fin.cases transitions_32000 (Fin.cases transitions_32128 (Fin.cases transitions_32256 (Fin.cases transitions_32384 (Fin.cases transitions_32512 (Fin.cases transitions_32640 (
  Fin.cases transitions_32768 (Fin.cases transitions_32896 (Fin.cases transitions_33024 (Fin.cases transitions_33152 (Fin.cases transitions_33280 (Fin.cases transitions_33408 (Fin.cases transitions_33536 (Fin.cases transitions_33664 (Fin.cases transitions_33792 (Fin.cases transitions_33920 (Fin.cases transitions_34048 (Fin.cases transitions_34176 (Fin.cases transitions_34304 (Fin.cases transitions_34432 (Fin.cases transitions_34560 (Fin.cases transitions_34688 (
  Fin.cases transitions_34816 (Fin.cases transitions_34944 (Fin.cases transitions_35072 (Fin.cases transitions_35200 (Fin.cases transitions_35328 (Fin.cases transitions_35456 (Fin.cases transitions_35584 (Fin.cases transitions_35712 (Fin.cases transitions_35840 (Fin.cases transitions_35968 (Fin.cases transitions_36096 (Fin.cases transitions_36224 (Fin.cases transitions_36352 (Fin.cases transitions_36480 (Fin.cases transitions_36608 (Fin.cases transitions_36736 (
  Fin.cases transitions_36864 (Fin.cases transitions_36992 (Fin.cases transitions_37120 (Fin.cases transitions_37248 (Fin.cases transitions_37376 (Fin.cases transitions_37504 (Fin.cases transitions_37632 (Fin.cases transitions_37760 (Fin.cases transitions_37888 (Fin.cases transitions_38016 (Fin.cases transitions_38144 (Fin.cases transitions_38272 (Fin.cases transitions_38400 (Fin.cases transitions_38528 (Fin.cases transitions_38656 (Fin.cases transitions_38784 (
  Fin.cases transitions_38912 (Fin.cases transitions_39040 (Fin.cases transitions_39168 (Fin.cases transitions_39296 (Fin.cases transitions_39424 (Fin.cases transitions_39552 (Fin.cases transitions_39680 (Fin.cases transitions_39808 (Fin.cases transitions_39936 (Fin.cases transitions_40064 (Fin.cases transitions_40192 (Fin.cases transitions_40320 (Fin.cases transitions_40448 (Fin.cases transitions_40576 (Fin.cases transitions_40704 (Fin.cases transitions_40832 (
  Fin.cases transitions_40960 (Fin.cases transitions_41088 (Fin.cases transitions_41216 (
  (fun b : Fin 0 => Fin.elim0 b))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))

private theorem mirror_fields_accepted (j : Fin 41294) : checkMirrorFields j = true := by
  let b : Fin 323 := ⟨j.val / 128, by have hj := j.isLt; omega⟩
  let k : Fin (chunkLength b) := ⟨j.val % 128, by
    have hj := j.isLt
    dsimp [chunkLength, b]
    omega⟩
  have he : chunkIndex b k = j := by
    apply Fin.ext
    dsimp [chunkIndex, b, k]
    omega
  simpa only [he] using mirror_blocks b k

theorem mirror_accepted (j : Fin 41294) : checkMirror j = true := by
  apply decide_eq_true
  exact EJMLocalCacheData.checked_raw_mirror_eq (by decide) (by decide)
    (mirror j) (triangleTree.read 10 j.val) (mirror_fields_accepted j)


theorem transitions_accepted (j : Fin 41294) : checkTransitions j = true := by
  let b : Fin 323 := ⟨j.val / 128, by have hj := j.isLt; omega⟩
  let k : Fin (chunkLength b) := ⟨j.val % 128, by
    have hj := j.isLt
    dsimp [chunkLength, b]
    omega⟩
  have he : chunkIndex b k = j := by
    apply Fin.ext
    dsimp [chunkIndex, b, k]
    omega
  simpa only [he] using transitions_blocks b k

theorem assignment_eq (j : Fin 41294) :
    assignment j = (triangleRows j).assignment :=
  (of_decide_eq_true (mirror_accepted j)).1

theorem lower_eq (j : Fin 41294) : lower j = (triangleRows j).lower :=
  (of_decide_eq_true (mirror_accepted j)).2

/-- This consequence uses the original row assignments on both sides. -/
theorem transition_accepted (j : Fin 41294) (i : Fin 10) (a : Fin 4)
    (unknown : (triangleRows j).assignment i = none) :
    D5.S0.Certificates.FiniteLocalRowCache.checkTransition trianglePartialAct
      triangleRows transition j i a = true := by
  have hu : assignment j i = none := by rw [assignment_eq, unknown]
  have h := of_decide_eq_true (transitions_accepted j) i hu a
  rw [assignment_eq, assignment_eq] at h
  exact decide_eq_true h

theorem accepted : (∀ j, checkMirror j = true) ∧ (∀ j, checkTransitions j = true) :=
  ⟨mirror_accepted, transitions_accepted⟩

#print axioms assignment_eq
#print axioms lower_eq
#print axioms transition_accepted

end D5.S3.Quantum.Entanglement.EJM.EJMTriangleCacheReceipts
