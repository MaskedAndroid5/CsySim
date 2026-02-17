from typing import Dict, List
import argparse
import json
import math
import random
import csv
import sys

class CsyBindingEngine:
    """
    Core Type I-F binding logic: PAM + seed strictness + Gaussian distal penalty.
    Base class for both quick prediction and full simulation.
    """
    def __init__(self, sigma: float = 3.0):
        self.sigma = sigma

    def gaussian_binding_prob(self, mismatches: int) -> float:
        """Normalized Gaussian: P(0) = 1.0, falls off with mismatch count."""
        if mismatches < 0:
            return 0.0
        exponent = -(mismatches ** 2) / (2 * self.sigma ** 2)
        prob = math.exp(exponent)
        return min(max(prob, 0.0), 1.0)

    def evaluate_target(
        self,
        crrna_spacer: str,
        target_protospacer: str,
        pam_flank_5prime: str
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
            "mismatch_positions": [],  # now 1-based
            "total_mismatches": 0,
            "binding": "none",
            "rloop_formation": "unknown",
            "interference": False,
            "binding_prob": 0.0,
            "reason": ""
        }

        # Input validation
        valid_bases = set("ACGT")
        if not all(c in valid_bases for c in crrna_spacer + target_protospacer + pam_flank_5prime):
            result["reason"] = "Invalid characters (only ACGT allowed)"
            return result

        # 1. PAM check (strict 5'-GG on non-target strand)
        if len(pam_flank_5prime) < 2 or pam_flank_5prime[-2:] != "GG":
            result["reason"] = "Invalid PAM (requires 5'-..GG on non-target strand)"
            return result

        # 2. Length match
        if len(crrna_spacer) != len(target_protospacer):
            result["reason"] = "Spacer and protospacer length mismatch"
            return result

        # 3. Mismatches (report 1-based)
        mismatches_0 = [i for i in range(len(target_protospacer)) if target_protospacer[i] != crrna_spacer[i]]
        result["mismatch_positions"] = [i+1 for i in mismatches_0]
        result["total_mismatches"] = len(mismatches_0)

        # 4. Seed region (positions 1–8) must be perfect
        seed_mismatches = [pos+1 for pos in mismatches_0 if pos < 8]
        if seed_mismatches:
            result["reason"] = f"Seed mismatch at position(s): {', '.join(map(str, seed_mismatches))}"
            return result

        # 5. Distal penalty (Gaussian)
        prob = self.gaussian_binding_prob(result["total_mismatches"])
        result["binding_prob"] = round(prob, 4)

        if result["total_mismatches"] == 0:
            result.update({"binding": "strong", "rloop_formation": "full_stable", "interference": True,
                           "reason": "Perfect match → high-affinity binding expected"})
        elif result["total_mismatches"] <= 3:
            result.update({"binding": "moderate", "rloop_formation": "stable_with_bubbles", "interference": True,
                           "reason": "Few distal mismatches → stable R-loop likely"})
        elif result["total_mismatches"] <= 6:
            result.update({"binding": "weak", "rloop_formation": "partial_unstable", "interference": False,
                           "reason": "Moderate distal mismatches → weak/unstable binding"})
        else:
            result["reason"] = "Excessive mismatches → no stable binding"

        return result


# Simplified RNA:DNA NN ΔG parameters (kcal/mol, 37°C, approximate from hybrid literature)
NN_DG_HYBRID = {
    'AA': -0.9, 'AC': -2.3, 'AG': -2.1, 'AU': -1.1,
    'CA': -2.0, 'CC': -3.3, 'CG': -2.4, 'CU': -1.4,
    'GA': -1.6, 'GC': -3.4, 'GG': -2.9, 'GU': -1.5,
    'UA': -1.1, 'UC': -2.4, 'UG': -1.5, 'UU': -0.5,
}

MISMATCH_PENALTY = 3.2  # per mismatch stack penalty

def calculate_segment_dg(crrna_seg: str, target_seg: str) -> float:
    """Approximate ΔG for RNA:DNA hybrid segment (kcal/mol, 37°C)."""
    if len(crrna_seg) != len(target_seg) or len(crrna_seg) < 2:
        return 0.0

    dg = 4.1  # initiation penalty
    mm_count = 0

    for i in range(len(crrna_seg) - 1):
        r1, r2 = crrna_seg[i], crrna_seg[i+1]
        d1, d2 = target_seg[i], target_seg[i+1]
        key = r1 + r2
        if r1 == d1 and r2 == d2:
            dg += NN_DG_HYBRID.get(key, -1.5)
        else:
            dg += MISMATCH_PENALTY
            mm_count += 1

    # Terminal AU penalty
    if crrna_seg[-1] in 'AU' or target_seg[-1] in 'AU':
        dg += 0.45

    return dg / len(crrna_seg) if len(crrna_seg) > 0 else 0.0


class CsySimulator(CsyBindingEngine):
    """
    Monte Carlo simulator for directional, segmental R-loop unzipping along Csy3 backbone.
    Captures stochastic propagation, structural sensitivity, and Cas3 recruitment probability.
    """
    SEGMENT_DESCRIPTIONS = [
        "Csy3.1 – PAM-proximal (near Csy4/Csy2, Arg thumb loop)",
        "Csy3.2 – Helical core (Lys contact zone)",
        "Csy3.3 – Mid-backbone (potential thumb disruption)",
        "Csy3.4 – Distal curvature (Arg sensitive to bubbles)",
        "Csy3.5 – Pre-Cas3 zone (Glu loop region)",
        "Csy3.6 – PAM-distal tail (flexible His region)",
    ]

    def __init__(
        self,
        sigma_distal: float = 2.5,
        trials: int = 1000,
        segment_size: int = 6,
        dg_weight: float = 0.6,
        cas3_discount: float = 0.90
    ):
        super().__init__(sigma=sigma_distal)
        self.trials = trials
        self.segment_size = segment_size
        self.dg_weight = dg_weight
        self.cas3_discount = cas3_discount

    def simulate_unzipping(
        self,
        crrna_spacer: str,
        target_protospacer: str,
        pam_flank_5prime: str,
        use_dg: bool = False
    ) -> Dict:
        base = self.evaluate_target(crrna_spacer, target_protospacer, pam_flank_5prime)
        sim_result = {
            "trials": self.trials,
            "full_rloop_success_rate": 0.0,
            "avg_segments_unzipped": 0.0,
            "segment_mismatch_counts": [],
            "interference_probability": 0.0,
            "structural_notes": [],
            "detailed_reason": base["reason"]
        }

        if base["binding_prob"] == 0.0:
            sim_result["detailed_reason"] += " → simulation aborted (PAM or seed failure)"
            return {**base, "simulation": sim_result}

        # Segment the spacer (proximal → distal)
        spacer_len = len(crrna_spacer)
        num_segments = (spacer_len + self.segment_size - 1) // self.segment_size
        segment_mms = []
        segment_dgs = []

        for i in range(num_segments):
            start = i * self.segment_size
            end = min(start + self.segment_size, spacer_len)
            cr_seg = crrna_spacer[start:end]
            tg_seg = target_protospacer[start:end]
            mm_count = sum(a != b for a, b in zip(cr_seg, tg_seg))
            segment_mms.append(mm_count)
            segment_dgs.append(calculate_segment_dg(cr_seg, tg_seg) if use_dg else 0.0)

        sim_result["segment_mismatch_counts"] = segment_mms

        # Structural annotations
        for idx, mm in enumerate(segment_mms):
            if mm > 0 and idx < len(self.SEGMENT_DESCRIPTIONS):
                sim_result["structural_notes"].append(
                    f"Segment {idx+1} ({self.SEGMENT_DESCRIPTIONS[idx]}): {mm} mismatch(es) → potential local destabilization"
                )

        successes = 0
        total_opened = 0

        for _ in range(self.trials):
            opened = 0
            stalled = False
            for seg_idx in range(num_segments):
                mm = segment_mms[seg_idx]
                dg = segment_dgs[seg_idx]

                # Seed segments must be perfect
                if seg_idx * self.segment_size < 8:
                    if mm > 0:
                        stalled = True
                        break
                    prob = 1.0
                else:
                    gauss_p = self.gaussian_binding_prob(mm)
                    if use_dg and dg != 0:
                        # Sigmoid: more negative ΔG → higher prob (shifted/steepened)
                        dg_prob = 1 / (1 + math.exp((dg + 0.8) * 4))
                        prob = self.dg_weight * dg_prob + (1 - self.dg_weight) * gauss_p
                    else:
                        prob = gauss_p
                    prob = max(prob, 0.01)  # small escape floor

                if random.random() < prob:
                    opened += 1
                else:
                    stalled = True
                    break

            total_opened += opened
            if not stalled and opened == num_segments:
                successes += 1

        success_rate = successes / self.trials
        avg_opened = total_opened / self.trials

        sim_result.update({
            "full_rloop_success_rate": round(success_rate, 4),
            "avg_segments_unzipped": round(avg_opened, 2),
            "interference_probability": round(success_rate * self.cas3_discount, 4),
        })

        if success_rate > 0.80:
            level = "very high confidence full interference"
        elif success_rate > 0.50:
            level = "high probability interference"
        elif success_rate > 0.20:
            level = "moderate / stochastic interference"
        elif success_rate > 0.03:
            level = "low-probability escape possible"
        else:
            level = "unlikely to achieve interference"

        sim_result["detailed_reason"] = (
            f"Monte Carlo ({self.trials} trials): {level} "
            f"({success_rate*100:.1f}% full R-loops, avg {avg_opened:.1f}/{num_segments} segments unzipped). "
            + " ".join(sim_result["structural_notes"])
        )

        return {**base, "simulation": sim_result}


# ────────────────────────────────────────────────
# Tests + Batch
# ────────────────────────────────────────────────
def run_tests(sim: CsySimulator):
    spacer = "ATGCAGTCGGTACGTACGTACGTACGTACGTA"  # 32 nt
    cases = [
        ("Perfect match", spacer, spacer, "AAGG", False),
        ("Seed mismatch pos 5", spacer, "ATGCCGTCGGTACGTACGTACGTACGTACGTA", "AAGG", False),
        ("3 distal mismatches (pos 9,17,25)", spacer, "ATGCAGTCAGTACGTACGTACGTACGTACGTA", "AAGG", False),
        ("Invalid PAM", spacer, spacer, "AAAA", False),
        ("Same as above but with ΔG", spacer, "ATGCAGTCAGTACGTACGTACGTACGTACGTA", "AAGG", True),
    ]

    print("\n=== CsySim v0.2 Tests ===\n")
    for name, sp, tgt, pam, use_dg in cases:
        result = sim.simulate_unzipping(sp, tgt, pam, use_dg=use_dg)
        s = result["simulation"]
        print(f"→ {name}")
        print(f"  Legacy binding prob: {result['binding_prob']:.4f}")
        print(f"  Full R-loop success: {s['full_rloop_success_rate']*100:.1f}%")
        print(f"  Avg segments unzipped: {s['avg_segments_unzipped']}")
        print(f"  Interference prob:     {s['interference_probability']*100:.1f}%")
        print(f"  Structural notes:      {s['structural_notes']}")
        print(f"  Reason: {s['detailed_reason']}\n")


def batch_process(sim: CsySimulator, input_csv: str, output_csv: str, use_dg: bool):
    with open(input_csv, 'r', encoding='utf-8') as inf:
        reader = csv.DictReader(inf)
        if not {'spacer', 'target', 'pam'}.issubset(reader.fieldnames):
            print("Error: CSV must have columns: spacer, target, pam")
            sys.exit(1)
        rows = list(reader)

    out_rows = []
    for row in rows:
        spacer = row['spacer'].strip()
        target = row['target'].strip()
        pam    = row['pam'].strip()
        result = sim.simulate_unzipping(spacer, target, pam, use_dg=use_dg)
        s = result["simulation"]
        out_rows.append({
            'spacer': spacer,
            'target': target,
            'pam': pam,
            'total_mismatches': result['total_mismatches'],
            'mismatch_positions': ','.join(map(str, result['mismatch_positions'])),
            'full_rloop_success_rate': s['full_rloop_success_rate'],
            'avg_segments_unzipped': s['avg_segments_unzipped'],
            'interference_probability': s['interference_probability'],
            'reason': s['detailed_reason']
        })

    with open(output_csv, 'w', newline='', encoding='utf-8') as outf:
        writer = csv.DictWriter(outf, fieldnames=out_rows[0].keys())
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Batch processed {len(rows)} targets → {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CsySim v0.2 – Type I-F CRISPR-Cas Simulator")
    parser.add_argument("--spacer", help="crRNA spacer sequence")
    parser.add_argument("--target", help="Target protospacer sequence")
    parser.add_argument("--pam", help="5' flank on non-target strand (must end with GG)")
    parser.add_argument("--json", action="store_true", help="JSON output only")
    parser.add_argument("--test", action="store_true", help="Run built-in tests")
    parser.add_argument("--trials", type=int, default=1200, help="Monte Carlo trials")
    parser.add_argument("--seed", type=int, default=None, help="Random seed (for reproducibility)")
    parser.add_argument("--use_dg", action="store_true", help="Use biophysical ΔG scoring")
    parser.add_argument("--dg_weight", type=float, default=0.6, help="ΔG vs Gaussian weight (0–1)")
    parser.add_argument("--input_csv", help="Batch mode: input CSV (spacer,target,pam)")
    parser.add_argument("--output_csv", help="Batch mode: output CSV path")

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    sim = CsySimulator(
        trials=args.trials,
        sigma_distal=2.5,
        dg_weight=args.dg_weight
    )

    if args.test:
        run_tests(sim)
    elif args.input_csv and args.output_csv:
        batch_process(sim, args.input_csv, args.output_csv, args.use_dg)
    elif args.spacer and args.target and args.pam:
        result = sim.simulate_unzipping(args.spacer, args.target, args.pam, use_dg=args.use_dg)
        if args.json:
            print(json.dumps(result, indent=2))
        else:
            s = result["simulation"]
            print(f"CsySim v0.2 Simulation ({s['trials']} trials, ΔG={args.use_dg}, seed={args.seed})")
            print(f"Full R-loop success rate: {s['full_rloop_success_rate']*100:.1f}%")
            print(f"Avg segments unzipped:    {s['avg_segments_unzipped']}")
            print(f"Interference probability: {s['interference_probability']*100:.1f}%")
            print(f"Legacy binding prob:      {result['binding_prob']:.4f}")
            print(f"Structural notes:         {s['structural_notes']}")
            print(f"Reason:                   {s['detailed_reason']}")
            print(f"Segment mismatch counts:  {s['segment_mismatch_counts']}")
    else:
        parser.print_help()
