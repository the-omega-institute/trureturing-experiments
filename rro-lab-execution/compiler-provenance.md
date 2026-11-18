# Compiler provenance

`constructor_compiler.py` ports the local instruction definitions of
`Complexity.StackSimulation.controller` and the `seekMachine`,
`insertMachine`, `rewindMachine`, `deleteMachine`, `moveLeftMachine`,
`insertRewindMachine` and `deleteRewindMachine` compositions that it uses.
The supplier is RBarish-UTokyo/NPCompleteness-Lean4 immutable revision
`3e603e795b44fedcf7d6dd65b7d1e30547ead925`, Apache-2.0; the complete license
is `RBarishBinaryMachine-LICENSE.txt`. The immutable upstream archive hash
is `1c0f46a8ccc32a190cf74266fa0ab81ed91263b7b7d2155cedfe677b7d614230`.
The archive has no NOTICE file. Original mathematical authorship belongs
to that supplier. This executable derivative uses symbolic phase labels
and enumerates their reachable closure in place of the original `Fin`
encoding; its transitions preserve the source definitions' local rules.
The port does not include a Lean proof of this client or a new universal
machine theorem. The two-bit adapter and evaluator network are separate
ordinary repository-derived constructions.

The controller starts at an already serialized `StackEncoding.homeTape`.
The executable therefore uses the `controller_simulates` contract, not
the `StackCompile.initMachine` input initializer or a polynomial-time
corollary. A terminating stack computation gives a finite room bound;
all-input termination of this particular evaluator is proved separately
in the owning theory volume. No timing bound from the generic supplier
is silently transferred to the Wang/LAB material representation.

The local content port was read from primary HEAD
`af2733014daf45edc42a2b72b54b176c5dda0714`. These content hashes pin the
exact definitions used:

[
  {
    "path": "D5/S0/Computability/BinaryMachine/StackProgram.lean",
    "sha256": "0526d36fe30fd698b057847a6a8331b2de2ee77bda658b2abe641a4adc8a271a",
    "bytes": 19645
  },
  {
    "path": "D5/S0/Computability/BinaryMachine/StackCompiler.lean",
    "sha256": "ac3e3fcaae5151b977cdf2c027baa0fe1f131be803e9453d9cb82093ba731097",
    "bytes": 23496
  },
  {
    "path": "D5/S0/Computability/BinaryMachine/StackPop.lean",
    "sha256": "ae1319b4911f23ddcad1256f9add990d27dcfcf2acabfe65383af96a8e7b8552",
    "bytes": 11953
  },
  {
    "path": "D5/S0/Computability/BinaryMachine/StackEncoding.lean",
    "sha256": "69a77033bfdf955d70e62a3b6c27846da41cc80611c2e77b9e8a2a64063661fe",
    "bytes": 8873
  },
  {
    "path": "D5/S0/Computability/BinaryMachine/StackSimulation.lean",
    "sha256": "37320357c62908e3f1f3a1397249df3ca532ab090bb1adbbb2b0c5796382a2ce",
    "bytes": 45501
  },
  {
    "path": "D5/S0/Computability/BinaryMachine/StackCompile.lean",
    "sha256": "383e67b4b7d1f73a46fc87c55299716b775e1d06ac10511270ba98c036b54b76",
    "bytes": 6962
  },
  {
    "path": "D5/S0/Computability/BinaryMachine/Machine.lean",
    "sha256": "6801eb5323880c1e617323390009b1dbd2fdc8436b7e5e0e47e147a936c2c3c5",
    "bytes": 9856
  }
]
