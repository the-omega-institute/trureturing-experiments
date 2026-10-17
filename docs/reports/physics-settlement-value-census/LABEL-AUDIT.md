# Per-settlement label and source audit

Original and checked values are separate columns. `verified` means the label or repository read has a checked source; it does not certify the semantic judgment. Full excerpts, hashes, module sets and importer edges are retained in [labels-audited-v2.json](labels-audited-v2.json).

The checked table contains 244 verified cells, eight missing tracks and no remaining inferred cells. The original inferred shape at #15993 and interpreted conclusion at #14159 are marked below. Live proof-path use is not audited.

## PR #11696

[Settlement](https://github.com/the-omega-institute/trureturing/pull/11696); merge `0f12cab1756a514a04b127f4a86d8695e4229978`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L88](https://github.com/the-omega-institute/trureturing/pull/11696) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L92](https://github.com/the-omega-institute/trureturing/pull/11696) |
| Conclusion | Proved | Proved | verified | [MabkSelfTestingPositivity.scribe.cs L32](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/QuantumBounds/MabkSelfTestingPositivity.scribe.cs#L32) |
| Utility | none | none | verified | [MabkSelfTestingPositivity.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/QuantumBounds/MabkSelfTestingPositivity.lean#L6) |
| Search track | unlabelled | missing (null) | missing | [PR missing-track search](https://github.com/the-omega-institute/trureturing/pull/11696) |
| Merge period | early (2026-10-01..09) | early (2026-10-01..09) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/11696) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #12188

[Settlement](https://github.com/the-omega-institute/trureturing/pull/12188); merge `b5ffd994ba41dabb6e62d94a7804ad2b091b0a4b`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L92](https://github.com/the-omega-institute/trureturing/pull/12188) |
| Admission basis | escape-witness + open-problem-resolution | escape-witness + open-problem-resolution | verified | [PR body L72](https://github.com/the-omega-institute/trureturing/pull/12188) |
| Conclusion | Proved | Proved | verified | [SumRules.scribe.cs L120](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Algebra/ZeitlinSixJ/SumRules.scribe.cs#L120) |
| Utility | none | none | verified | [Alternating.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Algebra/ZeitlinSixJ/Alternating.lean#L6) |
| Search track | unlabelled | missing (null) | missing | [PR missing-track search](https://github.com/the-omega-institute/trureturing/pull/12188) |
| Merge period | early (2026-10-01..09) | early (2026-10-01..09) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/12188) |
| External D5 importer task | has external D5 importer (1) | has external D5 importer (1) | verified | [RicciLimit.lean L9](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Algebra/ZeitlinSixJ/RicciLimit.lean#L9) |

## PR #12833

[Settlement](https://github.com/the-omega-institute/trureturing/pull/12833); merge `4680987ec71024c9ac4701ee27302301aee51955`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L88](https://github.com/the-omega-institute/trureturing/pull/12833) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L88](https://github.com/the-omega-institute/trureturing/pull/12833) |
| Conclusion | Proved | Proved | verified | [GHZMeasureBiseparableBound.scribe.cs L44](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Entanglement/GHZMeasureBiseparableBound.scribe.cs#L44) |
| Utility | none | none | verified | [GHZMeasureBiseparableBound.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Entanglement/GHZMeasureBiseparableBound.lean#L6) |
| Search track | unlabelled | missing (null) | missing | [PR missing-track search](https://github.com/the-omega-institute/trureturing/pull/12833) |
| Merge period | early (2026-10-01..09) | early (2026-10-01..09) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/12833) |
| External D5 importer task | has external D5 importer (5) | has external D5 importer (6) | verified | [PrecessionUniversalThresholdRefutation.lean L60](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Entanglement/PrecessionUniversalThresholdRefutation.lean#L60) |

## PR #13018

[Settlement](https://github.com/the-omega-institute/trureturing/pull/13018); merge `d4950dc44cd6adbb69d2ed0f229578dc49eaf90a`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L9](https://github.com/the-omega-institute/trureturing/pull/13018) |
| Admission basis | escape-witness only | escape-witness only | verified | [PR body L10](https://github.com/the-omega-institute/trureturing/pull/13018) |
| Conclusion | none (library result) | none (library result) | verified | [PR body L5](https://github.com/the-omega-institute/trureturing/pull/13018) |
| Utility | none | none | verified | [CanonicalChannelDPI.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Divergence/CanonicalChannelDPI.lean#L6) |
| Search track | unlabelled | missing (null) | missing | [PR missing-track search](https://github.com/the-omega-institute/trureturing/pull/13018) |
| Merge period | early (2026-10-01..09) | early (2026-10-01..09) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/13018) |
| External D5 importer task | has external D5 importer (3) | has external D5 importer (3) | verified | [CanonicalPetzEquality.lean L9](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Recovery/CanonicalPetzEquality.lean#L9) |

## PR #13245

[Settlement](https://github.com/the-omega-institute/trureturing/pull/13245); merge `4db169edaf01299cdfb137013444d4cbfee93e51`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L39](https://github.com/the-omega-institute/trureturing/pull/13245) |
| Admission basis | escape-witness + open-problem-resolution | escape-witness + open-problem-resolution | verified | [PR body L35](https://github.com/the-omega-institute/trureturing/pull/13245) |
| Conclusion | Proved | Proved | verified | [AlaiStaircase.scribe.cs L26](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/QuantumBounds/MerminMeasurementDependence/AlaiStaircase.scribe.cs#L26) |
| Utility | none | none | verified | [AlaiStaircase.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/QuantumBounds/MerminMeasurementDependence/AlaiStaircase.lean#L6) |
| Search track | unlabelled | missing (null) | missing | [PR missing-track search](https://github.com/the-omega-institute/trureturing/pull/13245) |
| Merge period | early (2026-10-01..09) | early (2026-10-01..09) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/13245) |
| External D5 importer task | has external D5 importer (2) | has external D5 importer (2) | verified | [CliffordPathRealizations.lean L71](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Measurements/CliffordJointMeasurability/CliffordPathRealizations.lean#L71) |

## PR #13409

[Settlement](https://github.com/the-omega-institute/trureturing/pull/13409); merge `782c7a079c3f3ff0e096b557b4774fc853f01d61`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L90](https://github.com/the-omega-institute/trureturing/pull/13409) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L78](https://github.com/the-omega-institute/trureturing/pull/13409) |
| Conclusion | Refuted | Refuted | verified | [XYPaddedGinibreRefutation.scribe.cs L49](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/StatisticalMechanics/XYPaddedGinibreRefutation.scribe.cs#L49) |
| Utility | certified-instance | certified-instance | verified | [XYPaddedGinibreRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/StatisticalMechanics/XYPaddedGinibreRefutation.lean#L6) |
| Search track | deep | deep | verified | [PR body L11](https://github.com/the-omega-institute/trureturing/pull/13409) |
| Merge period | early (2026-10-01..09) | early (2026-10-01..09) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/13409) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #13680

[Settlement](https://github.com/the-omega-institute/trureturing/pull/13680); merge `deb4a9ec28fe04a7f7e360683a88906e6f867901`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L93](https://github.com/the-omega-institute/trureturing/pull/13680) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L95](https://github.com/the-omega-institute/trureturing/pull/13680) |
| Conclusion | Refuted | Refuted | verified | [PrecessionUniversalThresholdRefutation.scribe.cs L39](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Entanglement/PrecessionUniversalThresholdRefutation.scribe.cs#L39) |
| Utility | certified-instance | certified-instance | verified | [PrecessionUniversalThresholdRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Entanglement/PrecessionUniversalThresholdRefutation.lean#L6) |
| Search track | unlabelled | missing (null) | missing | [PR missing-track search](https://github.com/the-omega-institute/trureturing/pull/13680) |
| Merge period | early (2026-10-01..09) | early (2026-10-01..09) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/13680) |
| External D5 importer task | has external D5 importer (1) | has external D5 importer (1) | verified | [GeometricTraceConvexity.lean L12](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Hiai/GeometricTraceConvexity.lean#L12) |

## PR #14064

[Settlement](https://github.com/the-omega-institute/trureturing/pull/14064); merge `f58a634fd626c700af9438262a63f9d453372f90`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L230](https://github.com/the-omega-institute/trureturing/pull/14064) |
| Admission basis | escape-witness only | escape-witness only | verified | [PR body L224](https://github.com/the-omega-institute/trureturing/pull/14064) |
| Conclusion | Proved | Proved | verified | [HiguchiSudberyEntropyMaximum.scribe.cs L39](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Entanglement/HiguchiSudbery/HiguchiSudberyEntropyMaximum.scribe.cs#L39) |
| Utility | checker | checker | verified | [PrimeHierarchyCertificate.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Entanglement/HiguchiSudbery/PrimeHierarchyCertificate.lean#L6) |
| Search track | deep | deep | verified | [PR body L10](https://github.com/the-omega-institute/trureturing/pull/14064) |
| Merge period | early (2026-10-01..09) | early (2026-10-01..09) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/14064) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #14159

[Settlement](https://github.com/the-omega-institute/trureturing/pull/14159); merge `330c5cb57e09dfeb96189920f9a5e6386b0b68b2`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L92](https://github.com/the-omega-institute/trureturing/pull/14159) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L64](https://github.com/the-omega-institute/trureturing/pull/14159) |
| Conclusion | Proved (interpreted from answer) | Proved | verified | [MayerZeroFidelityTightness.scribe.cs L83](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/QuantumChannels/MayerZeroFidelityTightness.scribe.cs#L83) |
| Utility | none | none | verified | [MayerZeroFidelityTightness.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/QuantumChannels/MayerZeroFidelityTightness.lean#L6) |
| Search track | deep | deep | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/14159) |
| Merge period | early (2026-10-01..09) | early (2026-10-01..09) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/14159) |
| External D5 importer task | has external D5 importer (1) | has external D5 importer (1) | verified | [LiJiangPolarPairOptimalityRefutation.lean L10](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Recovery/LiJiangPolarPairOptimalityRefutation.lean#L10) |

## PR #14605

[Settlement](https://github.com/the-omega-institute/trureturing/pull/14605); merge `934b24fe043ba4d5b63ffa610c516b153cc0e602`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L142](https://github.com/the-omega-institute/trureturing/pull/14605) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L142](https://github.com/the-omega-institute/trureturing/pull/14605) |
| Conclusion | Proved | Proved | verified | [SpinCoupledEvenHypercubeGroverFixedPoint.scribe.cs L45](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Dynamics/SpinCoupledEvenHypercubeGroverFixedPoint.scribe.cs#L45) |
| Utility | none | none | verified | [SpinCoupledEvenHypercubeGroverFixedPoint.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Dynamics/SpinCoupledEvenHypercubeGroverFixedPoint.lean#L6) |
| Search track | quick | quick | verified | [PR body L8](https://github.com/the-omega-institute/trureturing/pull/14605) |
| Merge period | early (2026-10-01..09) | early (2026-10-01..09) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/14605) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #15193

[Settlement](https://github.com/the-omega-institute/trureturing/pull/15193); merge `8799c189f58c905b82c725e127af6ca569781caf`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L74](https://github.com/the-omega-institute/trureturing/pull/15193) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L74](https://github.com/the-omega-institute/trureturing/pull/15193) |
| Conclusion | Proved | Proved | verified | [ThreeOutcomePermutationallyInvariantBell.scribe.cs L63](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Entanglement/ThreeOutcomePermutationallyInvariantBell.scribe.cs#L63) |
| Utility | none | none | verified | [ThreeOutcomePermutationallyInvariantBell.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Entanglement/ThreeOutcomePermutationallyInvariantBell.lean#L6) |
| Search track | quick | quick | verified | [PR body L10](https://github.com/the-omega-institute/trureturing/pull/15193) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/15193) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #15354

[Settlement](https://github.com/the-omega-institute/trureturing/pull/15354); merge `bcb25c259cf7d9e938443ad5f71ec7648520658f`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L146](https://github.com/the-omega-institute/trureturing/pull/15354) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L78](https://github.com/the-omega-institute/trureturing/pull/15354) |
| Conclusion | Proved | Proved | verified | [RenyiInformationCombiningOrderThree.scribe.cs L98](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Information/RenyiInformationCombiningOrderThree.scribe.cs#L98) |
| Utility | none | none | verified | [RenyiInformationCombiningOrderThree.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Information/RenyiInformationCombiningOrderThree.lean#L6) |
| Search track | deep | deep | verified | [PR body L9](https://github.com/the-omega-institute/trureturing/pull/15354) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/15354) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #15530

[Settlement](https://github.com/the-omega-institute/trureturing/pull/15530); merge `715e13e1ac0bffee9026ac384ad9817bb184c467`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L233](https://github.com/the-omega-institute/trureturing/pull/15530) |
| Admission basis | escape-witness + open-problem-resolution | escape-witness + open-problem-resolution | verified | [PR body L53](https://github.com/the-omega-institute/trureturing/pull/15530) |
| Conclusion | Refuted | Refuted | verified | [WassersteinComplexityNonconvex.scribe.cs L21](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/QuantumChannels/WassersteinComplexityNonconvex.scribe.cs#L21) |
| Utility | certified-instance | certified-instance | verified | [WassersteinComplexityNonconvex.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/QuantumChannels/WassersteinComplexityNonconvex.lean#L6) |
| Search track | quick | quick | verified | [PR body L10](https://github.com/the-omega-institute/trureturing/pull/15530) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/15530) |
| External D5 importer task | has external D5 importer (3) | has external D5 importer (3) | verified | [PurificationLowerBoundQutritRefutation.lean L9](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Dynamics/PurificationLowerBoundQutritRefutation.lean#L9) |

## PR #15565

[Settlement](https://github.com/the-omega-institute/trureturing/pull/15565); merge `6206db20cd068d82a97ac8988289b53106808f7f`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L105](https://github.com/the-omega-institute/trureturing/pull/15565) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L105](https://github.com/the-omega-institute/trureturing/pull/15565) |
| Conclusion | Refuted | Refuted | verified | [QuiverTreeIsingRefutation.scribe.cs L73](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Zeros/LeeYang/QuiverTreeIsingRefutation.scribe.cs#L73) |
| Utility | certified-instance | certified-instance | verified | [QuiverTreeIsingRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Zeros/LeeYang/QuiverTreeIsingRefutation.lean#L6) |
| Search track | quick | quick | verified | [PR body L10](https://github.com/the-omega-institute/trureturing/pull/15565) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/15565) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #15993

[Settlement](https://github.com/the-omega-institute/trureturing/pull/15993); merge `dfce518d4ed0176c6d24e7e321f5d8251e34f1b5`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only (inferred) | bind-only | verified | [QuantumTriangularJensenShannonRefutation.lean L27](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Divergence/QuantumTriangularJensenShannonRefutation.lean#L27) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L35](https://github.com/the-omega-institute/trureturing/pull/15993) |
| Conclusion | Refuted | Refuted | verified | [QuantumTriangularJensenShannonRefutation.scribe.cs L72](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Divergence/QuantumTriangularJensenShannonRefutation.scribe.cs#L72) |
| Utility | certified-instance | certified-instance | verified | [QuantumTriangularJensenShannonRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Divergence/QuantumTriangularJensenShannonRefutation.lean#L6) |
| Search track | deep | deep | verified | [PR body L4](https://github.com/the-omega-institute/trureturing/pull/15993) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/15993) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #16103

[Settlement](https://github.com/the-omega-institute/trureturing/pull/16103); merge `7a7b3d8aca1edfb066e44078674c32ab79a104e2`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L76](https://github.com/the-omega-institute/trureturing/pull/16103) |
| Admission basis | escape-witness + open-problem-resolution | escape-witness + open-problem-resolution | verified | [PR body L70](https://github.com/the-omega-institute/trureturing/pull/16103) |
| Conclusion | Refuted | Refuted | verified | [RenyiBoundaryMagicAdditivityRefutation.scribe.cs L127](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Magic/RenyiBoundaryMagicAdditivityRefutation.scribe.cs#L127) |
| Utility | certified-instance | certified-instance | verified | [RenyiBoundaryMagicAdditivityRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Magic/RenyiBoundaryMagicAdditivityRefutation.lean#L6) |
| Search track | deep | deep | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/16103) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/16103) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #16361

[Settlement](https://github.com/the-omega-institute/trureturing/pull/16361); merge `49b529ada15fe5a8c033d48bcf52c36f2e4d827f`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L29](https://github.com/the-omega-institute/trureturing/pull/16361) |
| Admission basis | escape-witness only | escape-witness only | verified | [PR body L18](https://github.com/the-omega-institute/trureturing/pull/16361) |
| Conclusion | none (library result) | none (library result) | verified | [PR body L18](https://github.com/the-omega-institute/trureturing/pull/16361) |
| Utility | none | none | verified | [CanonicalRelativeModularModel.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Divergence/CanonicalRelativeModularModel.lean#L6) |
| Search track | unlabelled | missing (null) | missing | [PR missing-track search](https://github.com/the-omega-institute/trureturing/pull/16361) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/16361) |
| External D5 importer task | has external D5 importer (2) | has external D5 importer (2) | verified | [CanonicalChannelDPI.lean L9](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Divergence/CanonicalChannelDPI.lean#L9) |

## PR #16702

[Settlement](https://github.com/the-omega-institute/trureturing/pull/16702); merge `6d5b8c368c02ce0568d1bcd2ab3c290229745d23`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L66](https://github.com/the-omega-institute/trureturing/pull/16702) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L65](https://github.com/the-omega-institute/trureturing/pull/16702) |
| Conclusion | Proved | Proved | verified | [SteeringEllipsoidCHSH.scribe.cs L50](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Entanglement/SteeringEllipsoidCHSH.scribe.cs#L50) |
| Utility | none | none | verified | [SteeringEllipsoidCHSH.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Entanglement/SteeringEllipsoidCHSH.lean#L6) |
| Search track | deep | deep | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/16702) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/16702) |
| External D5 importer task | has external D5 importer (1) | has external D5 importer (1) | verified | [UniversalTeleportationCriticalConcurrence.lean L9](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Entanglement/UniversalTeleportationCriticalConcurrence.lean#L9) |

## PR #16709

[Settlement](https://github.com/the-omega-institute/trureturing/pull/16709); merge `af9cb36df074c1e6245572f87b34aea33d534f6e`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L100](https://github.com/the-omega-institute/trureturing/pull/16709) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L102](https://github.com/the-omega-institute/trureturing/pull/16709) |
| Conclusion | Refuted | Refuted | verified | [FourMUBBoundEntanglementRefutation.scribe.cs L29](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Entanglement/FourMUBBoundEntanglementRefutation.scribe.cs#L29) |
| Utility | certified-instance | certified-instance | verified | [FourMUBBoundEntanglementRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Entanglement/FourMUBBoundEntanglementRefutation.lean#L6) |
| Search track | quick | quick | verified | [PR body L10](https://github.com/the-omega-institute/trureturing/pull/16709) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/16709) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #16727

[Settlement](https://github.com/the-omega-institute/trureturing/pull/16727); merge `5423f37863d8c9ac7bed0e2be8836d8bb1501562`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L67](https://github.com/the-omega-institute/trureturing/pull/16727) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L69](https://github.com/the-omega-institute/trureturing/pull/16727) |
| Conclusion | Refuted | Refuted | verified | [PositiveKrausSideEntanglementRefutation.scribe.cs L31](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/QuantumChannels/PositiveKrausSideEntanglementRefutation.scribe.cs#L31) |
| Utility | certified-instance | certified-instance | verified | [PositiveKrausSideEntanglementRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/QuantumChannels/PositiveKrausSideEntanglementRefutation.lean#L6) |
| Search track | quick | quick | verified | [PR body L8](https://github.com/the-omega-institute/trureturing/pull/16727) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/16727) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #16879

[Settlement](https://github.com/the-omega-institute/trureturing/pull/16879); merge `d5c6da8f60fc3f5e184760fa515e1f5c6d0839e3`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L99](https://github.com/the-omega-institute/trureturing/pull/16879) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L99](https://github.com/the-omega-institute/trureturing/pull/16879) |
| Conclusion | Proved | Proved | verified | [DetectorSensorDickeOptimality.scribe.cs L69](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Measurement/DetectorSensorDickeOptimality.scribe.cs#L69) |
| Utility | none | none | verified | [DetectorSensorDickeOptimality.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Measurement/DetectorSensorDickeOptimality.lean#L6) |
| Search track | unlabelled | missing (null) | missing | [PR missing-track search](https://github.com/the-omega-institute/trureturing/pull/16879) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/16879) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #16887

[Settlement](https://github.com/the-omega-institute/trureturing/pull/16887); merge `1f3e4c76bba7710950ab3ba52846d5981241da63`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L75](https://github.com/the-omega-institute/trureturing/pull/16887) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L93](https://github.com/the-omega-institute/trureturing/pull/16887) |
| Conclusion | Proved | Proved | verified | [RepetitionVarianceLinear.scribe.cs L43](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/StatisticalMechanics/NoisyVoter/RepetitionVarianceLinear.scribe.cs#L43) |
| Utility | none | none | verified | [RepetitionVarianceLinear.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/StatisticalMechanics/NoisyVoter/RepetitionVarianceLinear.lean#L6) |
| Search track | quick | quick | verified | [PR body L13](https://github.com/the-omega-institute/trureturing/pull/16887) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/16887) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #16902

[Settlement](https://github.com/the-omega-institute/trureturing/pull/16902); merge `ecabdecf933fd9af31463bb7381808e6b356fa31`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L63](https://github.com/the-omega-institute/trureturing/pull/16902) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L51](https://github.com/the-omega-institute/trureturing/pull/16902) |
| Conclusion | Proved | Proved | verified | [UniversalTeleportationCriticalConcurrence.scribe.cs L55](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Entanglement/UniversalTeleportationCriticalConcurrence.scribe.cs#L55) |
| Utility | none | none | verified | [UniversalTeleportationCriticalConcurrence.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Entanglement/UniversalTeleportationCriticalConcurrence.lean#L6) |
| Search track | deep | deep | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/16902) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/16902) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #16949

[Settlement](https://github.com/the-omega-institute/trureturing/pull/16949); merge `3d20b114e87b4c87ea7c3f14c71643c751341429`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L86](https://github.com/the-omega-institute/trureturing/pull/16949) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L88](https://github.com/the-omega-institute/trureturing/pull/16949) |
| Conclusion | Proved | Proved | verified | [SymmetricGHZResourceLowerBound.scribe.cs L74](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Contextuality/SymmetricGHZResourceLowerBound.scribe.cs#L74) |
| Utility | none | none | verified | [SymmetricGHZResourceLowerBound.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Contextuality/SymmetricGHZResourceLowerBound.lean#L6) |
| Search track | quick | quick | verified | [PR body L8](https://github.com/the-omega-institute/trureturing/pull/16949) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/16949) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #16971

[Settlement](https://github.com/the-omega-institute/trureturing/pull/16971); merge `373017ce510c2fb2c27a727ca1e426f7c5c2d3ff`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L66](https://github.com/the-omega-institute/trureturing/pull/16971) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L56](https://github.com/the-omega-institute/trureturing/pull/16971) |
| Conclusion | Proved | Proved | verified | [TwoQubitSVDEntropyLowerBound.scribe.cs L93](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Information/TwoQubitSVDEntropyLowerBound.scribe.cs#L93) |
| Utility | none | none | verified | [TwoQubitSVDEntropyLowerBound.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Information/TwoQubitSVDEntropyLowerBound.lean#L6) |
| Search track | quick | quick | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/16971) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/16971) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #17103

[Settlement](https://github.com/the-omega-institute/trureturing/pull/17103); merge `9b8a8715d523932eea2429bb953c4ac8201bec93`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L192](https://github.com/the-omega-institute/trureturing/pull/17103) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L196](https://github.com/the-omega-institute/trureturing/pull/17103) |
| Conclusion | Proved | Proved | verified | [SchurChannelMixedUnitaryDistance.scribe.cs L38](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/QuantumChannels/SchurChannelMixedUnitaryDistance.scribe.cs#L38) |
| Utility | none | none | verified | [SchurChannelMixedUnitaryDistance.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/QuantumChannels/SchurChannelMixedUnitaryDistance.lean#L6) |
| Search track | deep | deep | verified | [PR body L8](https://github.com/the-omega-institute/trureturing/pull/17103) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/17103) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #17145

[Settlement](https://github.com/the-omega-institute/trureturing/pull/17145); merge `1dfc887b8b9d972ce299308b93cea336fcb01ade`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L87](https://github.com/the-omega-institute/trureturing/pull/17145) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L76](https://github.com/the-omega-institute/trureturing/pull/17145) |
| Conclusion | Proved | Proved | verified | [QutritComplementaryEntropyFrontier.scribe.cs L131](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Information/QutritComplementaryEntropyFrontier.scribe.cs#L131) |
| Utility | none | none | verified | [QutritComplementaryEntropyFrontier.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Information/QutritComplementaryEntropyFrontier.lean#L6) |
| Search track | deep | deep | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/17145) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/17145) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #17178

[Settlement](https://github.com/the-omega-institute/trureturing/pull/17178); merge `f10f2b259f007d557ad22c203037e13059d8ed06`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L53](https://github.com/the-omega-institute/trureturing/pull/17178) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L51](https://github.com/the-omega-institute/trureturing/pull/17178) |
| Conclusion | Refuted | Refuted | verified | [GeigerKedemSpinEntropyPhaseRefutation.scribe.cs L74](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Entropy/GeigerKedemSpinEntropyPhaseRefutation.scribe.cs#L74) |
| Utility | certified-instance | certified-instance | verified | [GeigerKedemSpinEntropyPhaseRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Entropy/GeigerKedemSpinEntropyPhaseRefutation.lean#L6) |
| Search track | quick | quick | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/17178) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/17178) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #17487

[Settlement](https://github.com/the-omega-institute/trureturing/pull/17487); merge `0403ff6b9d7bf2159085feb909bc6b1c956f2ceb`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L66](https://github.com/the-omega-institute/trureturing/pull/17487) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L58](https://github.com/the-omega-institute/trureturing/pull/17487) |
| Conclusion | Refuted | Refuted | verified | [NegativityReductionInclusionRefutation.scribe.cs L59](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/QuantumBounds/NegativityReductionInclusionRefutation.scribe.cs#L59) |
| Utility | certified-instance | certified-instance | verified | [NegativityReductionInclusionRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/QuantumBounds/NegativityReductionInclusionRefutation.lean#L6) |
| Search track | quick | quick | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/17487) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/17487) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #17489

[Settlement](https://github.com/the-omega-institute/trureturing/pull/17489); merge `3f89a8fc073a859012fba13e652362e5aebc16b8`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L65](https://github.com/the-omega-institute/trureturing/pull/17489) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L57](https://github.com/the-omega-institute/trureturing/pull/17489) |
| Conclusion | Refuted | Refuted | verified | [SpacetimeDensitySchattenSubadditivityRefutation.scribe.cs L66](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Information/SpacetimeDensitySchattenSubadditivityRefutation.scribe.cs#L66) |
| Utility | certified-instance | certified-instance | verified | [SpacetimeDensitySchattenSubadditivityRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Information/SpacetimeDensitySchattenSubadditivityRefutation.lean#L6) |
| Search track | quick | quick | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/17489) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/17489) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #17605

[Settlement](https://github.com/the-omega-institute/trureturing/pull/17605); merge `6be5a5159fc72d7d6c717c461ec79c965b04680b`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L69](https://github.com/the-omega-institute/trureturing/pull/17605) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L59](https://github.com/the-omega-institute/trureturing/pull/17605) |
| Conclusion | Proved | Proved | verified | [SymmetricQubitWitnessSpectrum.scribe.cs L102](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Entanglement/AbsoluteSeparability/SymmetricQubitWitnessSpectrum.scribe.cs#L102) |
| Utility | none | none | verified | [SymmetricQubitWitnessSpectrum.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Entanglement/AbsoluteSeparability/SymmetricQubitWitnessSpectrum.lean#L6) |
| Search track | deep | deep | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/17605) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/17605) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #17628

[Settlement](https://github.com/the-omega-institute/trureturing/pull/17628); merge `efb7828e4720425c521698a021dad7525f8c47ec`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L69](https://github.com/the-omega-institute/trureturing/pull/17628) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L61](https://github.com/the-omega-institute/trureturing/pull/17628) |
| Conclusion | Refuted | Refuted | verified | [PeriodicSshExcessPseudoEntropyRefutation.scribe.cs L99](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Fermionic/PeriodicSshExcessPseudoEntropyRefutation.scribe.cs#L99) |
| Utility | certified-instance | certified-instance | verified | [PeriodicSshExcessPseudoEntropyRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Fermionic/PeriodicSshExcessPseudoEntropyRefutation.lean#L6) |
| Search track | quick | quick | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/17628) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/17628) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #17653

[Settlement](https://github.com/the-omega-institute/trureturing/pull/17653); merge `1382eb6d44503e5394acc2750f0f527605a4104b`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L93](https://github.com/the-omega-institute/trureturing/pull/17653) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L105](https://github.com/the-omega-institute/trureturing/pull/17653) |
| Conclusion | Refuted | Refuted | verified | [UniqueBasisFormationRefutation.scribe.cs L30](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Entanglement/GraphDensity/UniqueBasisFormationRefutation.scribe.cs#L30) |
| Utility | certified-instance | certified-instance | verified | [UniqueBasisFormationRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Entanglement/GraphDensity/UniqueBasisFormationRefutation.lean#L6) |
| Search track | quick | quick | verified | [PR body L10](https://github.com/the-omega-institute/trureturing/pull/17653) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/17653) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #17760

[Settlement](https://github.com/the-omega-institute/trureturing/pull/17760); merge `d39191b50a436b8ce6060d543a11567ce0301d3b`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L97](https://github.com/the-omega-institute/trureturing/pull/17760) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L101](https://github.com/the-omega-institute/trureturing/pull/17760) |
| Conclusion | Refuted | Refuted | verified | [SaleemNoisyCurvatureRefutation.scribe.cs L46](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Information/SaleemNoisyCurvatureRefutation.scribe.cs#L46) |
| Utility | certified-instance | certified-instance | verified | [SaleemNoisyCurvatureRefutation.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Information/SaleemNoisyCurvatureRefutation.lean#L6) |
| Search track | quick | quick | verified | [PR body L10](https://github.com/the-omega-institute/trureturing/pull/17760) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/17760) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |

## PR #17778

[Settlement](https://github.com/the-omega-institute/trureturing/pull/17778); merge `291aff96b947d2d07723a5b845a27ec17affa2cc`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | bind-only | bind-only | verified | [PR body L66](https://github.com/the-omega-institute/trureturing/pull/17778) |
| Admission basis | open-problem-resolution only | open-problem-resolution only | verified | [PR body L60](https://github.com/the-omega-institute/trureturing/pull/17778) |
| Conclusion | Proved | Proved | verified | [ConditionalHockeyStickOrthogonalReduction.scribe.cs L189](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Divergence/ConditionalHockeyStickOrthogonalReduction.scribe.cs#L189) |
| Utility | none | none | verified | [ConditionalHockeyStickOrthogonalReduction.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Divergence/ConditionalHockeyStickOrthogonalReduction.lean#L6) |
| Search track | quick | quick | verified | [PR body L3](https://github.com/the-omega-institute/trureturing/pull/17778) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/17778) |
| External D5 importer task | has external D5 importer (7) | has external D5 importer (7) | verified | [SubsystemLanczosPositivityRefutation.lean L31](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Dynamics/SubsystemLanczosPositivityRefutation.lean#L31) |

## PR #17791

[Settlement](https://github.com/the-omega-institute/trureturing/pull/17791); merge `9ebd3574b2f56569e23f6815ebabf0a95bf4faa3`.

| Field | Original value | Checked value | Status | Source |
| --- | --- | --- | --- | --- |
| Proof shape | content | content | verified | [PR body L181](https://github.com/the-omega-institute/trureturing/pull/17791) |
| Admission basis | escape-witness + open-problem-resolution | escape-witness + open-problem-resolution | verified | [PR body L311](https://github.com/the-omega-institute/trureturing/pull/17791) |
| Conclusion | Proved | Proved | verified | [StrongAbsorptionCumulants.scribe.cs L25](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/Blueprint/D5/S3/Quantum/Transport/StrongAbsorptionCumulants.scribe.cs#L25) |
| Utility | none | none | verified | [Density.lean L6](https://github.com/the-omega-institute/trureturing/blob/53ec26edbed43dc4f4862d2db303ad80cb3dbd92/D5/S3/Quantum/Transport/StrongAbsorption/Density.lean#L6) |
| Search track | quick | quick | verified | [PR body L10](https://github.com/the-omega-institute/trureturing/pull/17791) |
| Merge period | late (2026-10-10..16) | late (2026-10-10..16) | verified | [PR merge metadata](https://github.com/the-omega-institute/trureturing/pull/17791) |
| External D5 importer task | no external D5 importer (0) | no external D5 importer (0) | verified | [Full fixed-snapshot scan](audit_census.py) |
