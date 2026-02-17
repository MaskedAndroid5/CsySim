from Bio.Seq import Seq
import argparse
import json
import math
from typing import Dict, List, Tuple

class CsyBindingEngine:
    """ Type I-F (Csy complex) binding predictor with Gaussian probabilistic scoring.
    Grounded in P. aeruginosa cryo-EM and functional data:
    - Strict 5'-GG PAM on non-target strand
    - Positions 1–8 (seed) must be perfect match
    - Distal mismatches tolerated up to ~6–8 positions
    - Binding Kd ~nM for perfect matches, drops sharply with seed errors
    """

    def __init__(self, sigma: float = 3.0, seed_length: int = 8):
        """ sigma: Controls distal mismatch tolerance. Literature-inspired: σ ≈ 3 allows ~3 mismatches with moderate prob, 6 with low. """
        self.sigma = sigma
        self.seed_length = seed_length

    def gaussian_binding_prob(self, total_mismatches: int) -> float:
        """ Gaussian probability density (normalized to 1.0 at 0 mismatches).
        Math: P(x) = exp( -(x - μ)^2 / (2σ^2) ) with μ=0
        Reasoning:
        - Mismatch count x treated as deviation from ideal (0 mismatches)
        - Quadratic penalty models energy barriers in R-loop formation
        - σ=3 chosen so P(3) ≈ 0.606 (moderate), P(6) ≈ 0.135 (weak), P(9) ≈ 0.011 (near none)
        """
        if total_mismatches < 0:
            return 0.0
        exponent = - (total_mismatches ** 2) / (2 * self.sigma ** 2)
        prob = math.exp(exponent)
        return min(max(prob, 0.0), 1.0)  # Clamp [0,1]

    def evaluate_target(
        self, crrna_spacer: str, target_protospacer: str, pam_flank_5prime: str
    ) -> Dict:
        crrna_spacer = crrna_spacer.upper().strip()
        target_protospacer = target_protospacer.upper().strip()
        pam_flank_5prime = pam_flank_5prime.upper().strip()

        result: Dict = {
            "input": {
                "crrna_spacer": crrna_spacer,
                "target_protospacer": target_protospacer,
                "pam_5prime_non_target": pam_flank_5prime,
            },
            "mismatch_positions": [],
            "total_mismatches": 0,
            "binding": "none",
            "rloop_formation": "unknown",
            "interference": False,
            "binding_prob": 0.0,
            "reason": ""
        }

        # Check valid bases
        valid_bases = set("ACGT")
        if not all(base in valid_bases for base in crrna_spacer + target_protospacer + pam_flank_5prime):
            result["reason"] = "Invalid characters (only ACGT allowed in sequences)"
            return result

        # 1. PAM validation (strict GG on non-target strand)
        if len(pam_flank_5prime) < 2 or pam_flank_5prime[-2:] != "GG":
            result["reason"] = "Invalid PAM (requires 5'-GG on non-target strand)"
            return result

        # 2. Length check
        if len(crrna_spacer) != len(target_protospacer):
            result["reason"] = "Spacer and protospacer length mismatch"
            return result

        # 3. Compute mismatches (direct comparison — protospacer in same orientation as spacer)
        mismatches = [
            i + 1 for i in range(len(target_protospacer))
            if target_protospacer[i] != crrna_spacer[i]
        ]
        result["mismatch_positions"] = mismatches
        result["total_mismatches"] = len(mismatches)

        # 4. Seed region (positions 1–seed_length)
        seed_mismatches = [pos for pos in mismatches if pos <= self.seed_length]
        if seed_mismatches:
            result["reason"] = f"Seed mismatch at positions {', '.join(str(p) for p in seed_mismatches)}"
            return result

        # 5. Overall classification with Gaussian math
        prob = self.gaussian_binding_prob(len(mismatches))
        result["binding_prob"] = round(prob, 4)

        if len(mismatches) == 0:
            result["binding"] = "strong"
            result["rloop_formation"] = "full_stable"
            result["interference"] = True
            result["reason"] = "Perfect match → high-affinity binding expected"
        elif len(mismatches) <= 3:
            result["binding"] = "moderate"
            result["rloop_formation"] = "stable_with_bubbles"
            result["interference"] = True
            result["reason"] = "Few distal mismatches → stable R-loop likely"
        elif len(mismatches) <= 6:
            result["binding"] = "weak"
            result["rloop_formation"] = "partial_unstable"
            result["interference"] = False
            result["reason"] = "Moderate distal mismatches → weak/unstable binding"
        else:
            result["binding"] = "none"
            result["rloop_formation"] = "dissociated"
            result["reason"] = "Excessive mismatches → no stable binding"

        return result

# ----------------------
# Validation Test Suite
# ----------------------
def run_validation_tests(engine: CsyBindingEngine):
        tests: List[Tuple[str, str, str, Dict]] = [
        # Test 1: Perfect match
        ("ATGCAGTCGGTACGTACGTACGTACGTACGTA",
         "ATGCAGTCGGTACGTACGTACGTACGTACGTA",
         "CTGG",
         {"binding": "strong", "binding_prob": 1.0000}),

        # Test 2: Seed mismatch at position 5
        ("ATGCAGTCGGTACGTACGTACGTACGTACGTA",
         "ATGCTGTCGGTACGTACGTACGTACGTACGTA",  # mismatch at pos 5 (A→T)
         "CTGG",
         {"binding": "none", "binding_prob": 0.0000}),

        # Test 3: Distal mismatches (3) → moderate
        ("ATGCAGTCGGTACGTACGTACGTACGTACGTA",
         "ATGCAGTCAGTAGGTAGGTACGTACGTACGTA",  # mismatches at pos 9,13,17
         "CTGG",
         {"binding": "moderate", "binding_prob": 0.6065}),

        # Test 4: Invalid PAM
        ("ATGCAGTCGGTACGTACGTACGTACGTACGTA",
         "ATGCAGTCGGTACGTACGTACGTACGTACGTA",
         "CTAA",
         {"binding": "none", "binding_prob": 0.0000}),

        # Test 5: 7 distal mismatches → weak binding, low prob
        ("ATGCAGTCGGTACGTACGTACGTACGTACGTA",
         "ATGCAGTCGGTACGTACGTACGTACGTACGTACGTAAGGG",  # 7 changes at end
         "CTGG",
         {"binding": "weak", "binding_prob": 0.0657}),
    ]

    print("\n=== Validation Tests ===")
    for i, (spacer, target, pam, expected) in enumerate(tests, 1):
        result = engine.evaluate_target(spacer, target, pam)
        prob_match = abs(result["binding_prob"] - expected["binding_prob"]) < 0.001
        strength_match = result["binding"] == expected["binding"]
        passed = prob_match and strength_match
        print(f"Test {i}: { 'PASSED' if passed else 'FAILED' }")
        print(f" Expected: {expected['binding']} (prob ~{expected['binding_prob']})")
        print(f" Actual: {result['binding']} (prob {result['binding_prob']:.4f})")
        print(f" Reason: {result['reason']}\n")

# CLI Interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CsySim Type I-F Binding Predictor")
    parser.add_argument("--spacer", required=True, help="crRNA spacer sequence")
    parser.add_argument("--target", required=True, help="Target protospacer sequence")
    parser.add_argument("--pam", required=True, help="5' PAM flank on non-target strand (ends with GG)")
    parser.add_argument("--json", action="store_true", help="Output JSON only")
    parser.add_argument("--test", action="store_true", help="Run validation tests")

    args = parser.parse_args()

    engine = CsyBindingEngine(sigma=3.0)  # Adjustable tolerance

    if args.test:
        run_validation_tests(engine)
    else:
        result = engine.evaluate_target(args.spacer, args.target, args.pam)
        if args.json:
            print(json.dumps(result, indent=2))
        else:
            print(f"Binding: {result['binding']} (Prob: {result['binding_prob']:.4f})")
            print(f"R-loop: {result['rloop_formation']}")
            print(f"Interference: {result['interference']}")
            print(f"Mismatches: {result['total_mismatches']} at {result['mismatch_positions']}")
            print(f"Reason: {result['reason']}")
