#!/usr/bin/env python3
"""
Filter and summarise IntaRNA predictions from all_predictions_combined.tsv

Thresholds:
  Strong   : dG < -15 kcal/mol   (high priority for experimental validation)
  Moderate : -15 <= dG < -10     (biologically plausible)
  Weak     : -10 <= dG < -7      (low confidence)
  Noise    : dG >= -7            (discarded)

Outputs (in results/sRNA_detection/07_intarna/filtered/):
  01_strong_hits.tsv              All strong predictions ranked by dG
  02_moderate_hits.tsv            All moderate predictions
  03_summary_per_srna.tsv         Per-sRNA stats (n strong, n moderate, best dG)
  04_top5_per_srna.tsv            Top 5 targets per sRNA
  05_P2_exp_all_ranked.tsv        P2 full ranked list (exp intergenic)

Run from project root:
    python3 scripts/07_intarna/04_filter_results.py
"""

import csv
import os
from pathlib import Path
from collections import defaultdict

# ============================================================================
# Config
# ============================================================================

INPUT_TSV  = "results/sRNA_detection/07_intarna/all_predictions_combined.tsv"
OUTPUT_DIR = "results/sRNA_detection/07_intarna/filtered"

STRONG_THRESHOLD   = -15.0
MODERATE_THRESHOLD = -10.0
WEAK_THRESHOLD     =  -7.0

# ============================================================================
# Load data
# ============================================================================

def load_predictions(tsv_file):
    rows = []
    with open(tsv_file, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            try:
                row['dG'] = float(row['E'])
            except (ValueError, KeyError):
                row['dG'] = 0.0
            rows.append(row)
    return rows


def classify(dg):
    if dg < STRONG_THRESHOLD:
        return 'Strong'
    elif dg < MODERATE_THRESHOLD:
        return 'Moderate'
    elif dg < WEAK_THRESHOLD:
        return 'Weak'
    else:
        return 'Noise'


def write_tsv(rows, output_file, fieldnames):
    with open(output_file, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote: {output_file}  ({len(rows)} rows)")


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 70)
    print("IntaRNA Results Filter")
    print("=" * 70)
    print()

    if not os.path.exists(INPUT_TSV):
        print(f"ERROR: {INPUT_TSV} not found.")
        print("Run 03_combine_results.sh first.")
        return

    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    rows = load_predictions(INPUT_TSV)
    print(f"Loaded {len(rows)} predictions from {INPUT_TSV}")
    print()

    # Add classification column
    for r in rows:
        r['confidence'] = classify(r['dG'])

    # Sort by dG ascending (most negative = strongest)
    rows.sort(key=lambda r: r['dG'])

    # =========================================================================
    # Count by priority and confidence
    # =========================================================================

    print("Distribution by priority:")
    priority_counts = defaultdict(lambda: defaultdict(int))
    for r in rows:
        priority_counts[r['priority']][r['confidence']] += 1

    for priority, counts in sorted(priority_counts.items()):
        total = sum(counts.values())
        print(f"  {priority}  (total={total})")
        for conf in ['Strong', 'Moderate', 'Weak', 'Noise']:
            if counts[conf]:
                print(f"    {conf:10s}: {counts[conf]:4d}")
    print()

    # =========================================================================
    # 01 Strong hits
    # =========================================================================

    strong = [r for r in rows if r['dG'] < STRONG_THRESHOLD]
    fields = ['priority', 'srna', 'target', 'dG', 'confidence',
              'start1', 'end1', 'start2', 'end2', 'hybridDP']
    write_tsv(strong, f"{OUTPUT_DIR}/01_strong_hits.tsv", fields)

    # =========================================================================
    # 02 Moderate hits
    # =========================================================================

    moderate = [r for r in rows
                if STRONG_THRESHOLD <= r['dG'] < MODERATE_THRESHOLD]
    write_tsv(moderate, f"{OUTPUT_DIR}/02_moderate_hits.tsv", fields)

    # =========================================================================
    # 03 Per-sRNA summary
    # =========================================================================

    srna_stats = defaultdict(lambda: {
        'priority': '', 'n_total': 0, 'n_strong': 0, 'n_moderate': 0,
        'best_dG': 0.0, 'best_target': '', 'targets_strong': []
    })

    for r in rows:
        s = srna_stats[r['srna']]
        s['priority']  = r['priority']
        s['n_total']  += 1
        s['n_strong']  += 1 if r['dG'] < STRONG_THRESHOLD   else 0
        s['n_moderate']+= 1 if STRONG_THRESHOLD <= r['dG'] < MODERATE_THRESHOLD else 0
        if r['dG'] < s['best_dG']:
            s['best_dG']     = r['dG']
            s['best_target'] = r['target']
        if r['dG'] < STRONG_THRESHOLD:
            s['targets_strong'].append(r['target'])

    summary_rows = []
    for srna, s in sorted(srna_stats.items(),
                          key=lambda x: x[1]['best_dG']):
        summary_rows.append({
            'srna':           srna,
            'priority':       s['priority'],
            'n_total':        s['n_total'],
            'n_strong':       s['n_strong'],
            'n_moderate':     s['n_moderate'],
            'best_dG':        round(s['best_dG'], 2),
            'best_target':    s['best_target'],
            'strong_targets': ';'.join(s['targets_strong'][:10]),
        })

    write_tsv(summary_rows,
              f"{OUTPUT_DIR}/03_summary_per_srna.tsv",
              ['srna', 'priority', 'n_total', 'n_strong', 'n_moderate',
               'best_dG', 'best_target', 'strong_targets'])

    # =========================================================================
    # 04 Top 5 per sRNA
    # =========================================================================

    top5_rows = []
    seen = defaultdict(int)
    for r in rows:
        if seen[r['srna']] < 5:
            top5_rows.append(r)
            seen[r['srna']] += 1

    write_tsv(top5_rows,
              f"{OUTPUT_DIR}/04_top5_per_srna.tsv",
              ['srna', 'priority', 'target', 'dG', 'confidence',
               'start1', 'end1', 'start2', 'end2'])

    # =========================================================================
    # 05 P2 exp full ranked list
    # =========================================================================

    p2 = [r for r in rows if r['priority'] == 'P2_intergenic_exp']
    write_tsv(p2,
              f"{OUTPUT_DIR}/05_P2_exp_all_ranked.tsv",
              ['srna', 'target', 'dG', 'confidence',
               'start1', 'end1', 'start2', 'end2'])

    # =========================================================================
    # Console summary
    # =========================================================================

    print()
    print("=" * 70)
    print("TOP 5 PER sRNA  (ranked by dG)")
    print("=" * 70)
    print()
    print(f"  {'sRNA':<30} {'Target':<20} {'dG':>12}  Confidence")
    print(f"  {'-'*30} {'-'*20} {'-'*12}  ----------")

    prev_srna = None
    for r in top5_rows:
        sep = "" if r['srna'] == prev_srna else "\n"
        print(f"{sep}  {r['srna']:<30} {r['target']:<20} {r['dG']:>12.2f}  {r['confidence']}")
        prev_srna = r['srna']

    print()
    print("=" * 70)
    print(f"Strong hits (dG < {STRONG_THRESHOLD}): {len(strong)}")
    print(f"  P1 leuC antisense : {sum(1 for r in strong if r['priority']=='P1_leuC_antisense')}")
    print(f"  P2 intergenic exp : {sum(1 for r in strong if r['priority']=='P2_intergenic_exp')}")
    print(f"  P3 intergenic sta : {sum(1 for r in strong if r['priority']=='P3_intergenic_sta')}")
    print()
    print(f"Outputs written to: {OUTPUT_DIR}/")
    print("=" * 70)


if __name__ == "__main__":
    main()
