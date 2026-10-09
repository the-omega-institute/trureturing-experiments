import D5.S3.Quantum.Measurement.FourQubitCompatibilityDegree
open Lean Elab Term
open scoped BigOperators
elab "lane_const " s:str : term => do
  let env ← getEnv
  let suffix := "D5.S3.Quantum.Measurement.FourQubitCompatibilityDegree." ++ s.getString
  let candidateNames := env.constants.toList.filter fun (n, _) =>
    n.toString == suffix || n.toString.endsWith ("." ++ suffix)
  let [(n, _)] := candidateNames | throwError "expected one compiled declaration for {suffix}"
  return mkConst n
open D5.S3.Geometry.FourVectorSignSumBound
open D5.S3.Quantum.Measurement.FourQubitCompatibilityDegree
open D5.S3.Quantum.Measurement.FourQubitParentConstruction
open scoped BigOperators
example : IsGreatest
    {q : ℝ | ∃ x : Fin 4 → EuclideanSpace ℝ (Fin 3),
      (∀ ε, ‖signedSum x ε‖ ≤ 1) ∧ q = ∑ i, ‖x i‖}
    (Real.sqrt 13 / 2) := by
  have hs : 0 < Real.sqrt 52 := Real.sqrt_pos.2 (by norm_num)
  have h13 : 0 < Real.sqrt 13 := Real.sqrt_pos.2 (by norm_num)
  constructor
  · refine ⟨fun i => (Real.sqrt 52)⁻¹ • (lane_const "optimum") i, ?_, ?_⟩
    · intro ε
      have heq : signedSum (fun i => (Real.sqrt 52)⁻¹ • (lane_const "optimum") i) ε =
          (Real.sqrt 52)⁻¹ • signedSum (lane_const "optimum") ε := by
        simp [signedSum, Finset.smul_sum, smul_smul, mul_comm]
      rw [heq, norm_smul, Real.norm_eq_abs, abs_of_pos (inv_pos.mpr hs)]
      calc
        _ ≤ (Real.sqrt 52)⁻¹ * Real.sqrt 52 :=
          mul_le_mul_of_nonneg_left ((signedSum_le_max (lane_const "optimum") ε).trans_eq (lane_const "optimum_max"))
            (inv_nonneg.mpr hs.le)
        _ = 1 := inv_mul_cancel₀ hs.ne'
    · simp_rw [norm_smul, Real.norm_eq_abs, abs_of_pos (inv_pos.mpr hs)]
      rw [← Finset.mul_sum, (lane_const "optimum_sum"), (lane_const "sqrt52_eq")]
      have hsq : (Real.sqrt 13)^2 = 13 := Real.sq_sqrt (by norm_num)
      field_simp
      nlinarith [hsq]
  · rintro q ⟨x, hx, rfl⟩
    have hm : maxNorm x ≤ 1 := Finset.sup'_le _ _ fun ε _ => hx ε
    exact (four_vector_inequality x).trans
      (by simpa using mul_le_mul_of_nonneg_left hm (by positivity : 0 ≤ Real.sqrt 13 / 2))

example : 1 / minCompatDegree = Real.sqrt 13 / 2 := by
  rw [D5.S3.Quantum.Measurement.FourQubitCompatibilityDegree.result]
  field_simp

example : IsLeast (compatDegree '' povmTuples) endpoint := lane_const "minimum_attained"
