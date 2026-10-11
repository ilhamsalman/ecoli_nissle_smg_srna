#!/usr/bin/env python3
"""Compare ANNOgesic cis-element Associated_CDS with phase-matched curated DEGs.

Run from the ecoli_nissle_smg_srna repository root:
  python3 scripts/09_cis_deg/compare_cis_elements_to_degs.py

Uses only Python 3 standard library. Does NOT calculate genomic overlap; matches
ANNOgesic's Associated_CDS locus tag to a DEG locus tag, optionally resolving
DEG symbols through the locally curated DEG mRNA FASTA headers.
"""
import csv
import re
from pathlib import Path
from collections import Counter, defaultdict

ROOT = Path.cwd()
DEG_ROOT = ROOT / 'results/sRNA_detection/06_intarna/curated_mrna_tables'
RNA_ROOT = ROOT / 'results/riboswitch_thermometer'
FASTA_ROOT = ROOT / 'results/sRNA_detection/06_intarna/sequences/mrnas'
OUTPUT = ROOT / 'results/riboswitch_thermometer/deg_comparison'

DEG_FILES = {'exp': DEG_ROOT / 'TableS1_EXP_curated.tsv',
             'sta': DEG_ROOT / 'TableS2_STA_curated.tsv'}
RNA_FILES = {
    ('exp', 'thermometer_module'): RNA_ROOT / 'exp/annotation_RNA_thermometer.csv',
    ('exp', 'riboswitch_module'): RNA_ROOT / 'exp/annotation_riboswitch.csv',
    ('sta', 'thermometer_module'): RNA_ROOT / 'sta/annotation_RNA_thermometer.csv',
    ('sta', 'riboswitch_module'): RNA_ROOT / 'sta/annotation_riboswitch.csv',
}
EXPECTED = {('exp', 'thermometer_module'): 43, ('exp', 'riboswitch_module'): 27,
            ('sta', 'thermometer_module'): 14, ('sta', 'riboswitch_module'): 11}
LOCI = re.compile(r'\bECOLIN[_-]?\d{4,6}\b', re.I)
FASTA_LOCUS = re.compile(r'\b(ECOLIN[_-]?\d{4,6})\|([^\s|]+)', re.I)


def clean_locus(s):
    m = LOCI.search(str(s))
    return m.group().upper().replace('-', '_') if m else None


def name_key(s):
    return re.sub(r'[^a-z0-9]', '', str(s).lower())


def read_tsv(path):
    if not path.is_file():
        raise FileNotFoundError(f'Missing input: {path}')
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError(f'No TSV header: {path}')
        rows = list(reader)
    return rows, reader.fieldnames


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def field_of(fields, candidates):
    lookup = {name_key(x): x for x in fields}
    for c in candidates:
        if name_key(c) in lookup:
            return lookup[name_key(c)]
    return None


def fasta_symbol_map(phase):
    # Externally prepared DEG mRNA FASTA header example:
    # >rcnA CP007799.1:2449147-2449847(+) ECOLIN_11935|rcnA logFC=-5.53
    fastadir = FASTA_ROOT / ('deg_exp' if phase == 'exp' else 'deg_sta')
    syms = defaultdict(set)
    if not fastadir.exists():
        return syms, fastadir
    for p in fastadir.glob('*.fasta'):
        with p.open(encoding='utf-8', errors='replace') as f:
            h = f.readline().strip()
        m = FASTA_LOCUS.search(h)
        if m:
            locus, symbol = clean_locus(m.group(1)), m.group(2)
            syms[symbol.casefold()].add(locus)
            syms[p.stem.casefold()].add(locus)
    return syms, fastadir


def deg_lookup(phase, path):
    rows, fields = read_tsv(path)
    print(f'\n{phase.upper()} DEG input: {path}\n  columns: {fields}\n  rows: {len(rows)}')
    symbol_col = field_of(fields, [
        'gene_symbol', 'gene_name', 'gene symbols', 'symbol', 'Gene', 'Gene Name', 'Gene.symbol',
        'Gene name', 'Name', 'gene'])
    fc_col = field_of(fields, ['log2FC', 'logFC', 'log2FoldChange',
                               'log2_fold_change', 'log2 fold change', 'log2 foldchange'])
    p_col = field_of(fields, ['p_value', 'pvalue', 'P.Value', 'P value', 'pval', 'P'])
    padj_col = field_of(fields, ['padj', 'adj.P.Val', 'adjusted_p_value', 'FDR', 'qvalue'])
    mapping, fastadir = fasta_symbol_map(phase)
    print(f'  symbol column: {symbol_col or "(none)"}; log2FC column: {fc_col or "(none)"}')
    print(f'  fallback FASTA directory: {fastadir} (exists={fastadir.is_dir()})')
    hits = defaultdict(list)
    audit = []
    methods = Counter()
    for i, row in enumerate(rows, 2):
        # Prefer explicit stable locus tags anywhere in the DEG row; only
        # fall back to exact gene symbol -> stable locus map if necessary.
        ids = {x for v in row.values() if v for x in [clean_locus(v)] if x}
        symbol = str(row.get(symbol_col, '')).strip() if symbol_col else ''
        method = ''
        if len(ids) == 1:
            method = 'direct_locus'
        elif not ids and symbol:
            ids = mapping.get(symbol.casefold(), set())
            if len(ids) == 1:
                method = 'unique_symbol_via_fasta'
        if len(ids) != 1:
            method = 'ambiguous' if len(ids) > 1 else 'unresolved'
        locus = next(iter(ids)) if len(ids) == 1 else ''
        methods[method] += 1
        audit.append({'phase': phase, 'deg_row': i, 'gene_symbol': symbol,
                      'resolved_locus': locus, 'resolution': method,
                      'possible_loci': ','.join(sorted(ids))})
        if locus:
            hits[locus].append({'gene_symbol': symbol,
                                'log2FC': row.get(fc_col, '') if fc_col else '',
                                'p_value': row.get(p_col, '') if p_col else '',
                                'adjusted_p': row.get(padj_col, '') if padj_col else '',
                                'deg_row': i, 'method': method})
    print(f'  DEG row resolution: {dict(methods)}')
    if not hits:
        raise ValueError(f'No DEG locus tags resolved for {phase}. Check the printed columns and '
                         f'FASTA directory; do not rely on gene symbols without a validated map.')
    return hits, audit, len(rows), methods


def main():
    lookups = {}
    audits = []
    totals = {}
    for phase, file in DEG_FILES.items():
        lookups[phase], audit, totals[phase], _ = deg_lookup(phase, file)
        audits.extend(audit)
    all_records = []
    matches = []
    summaries = []
    for (phase, module), path in RNA_FILES.items():
        ann, cols = read_tsv(path)
        if len(ann) != EXPECTED[(phase, module)]:
            raise ValueError(f'Unexpected number of rows in {path}: '
                             f'{len(ann)} vs expected {EXPECTED[(phase, module)]}. '
                             'Review source files and do not silently accept a phase swap.')
        for col in ['#ID', 'Associated_CDS', 'Rfam_name', 'E_value']:
            if col not in cols:
                raise ValueError(f'Missing {col} in {path}; columns={cols}')
        n_match = 0
        matched_genes = set()
        for a in ann:
            locus = clean_locus(a['Associated_CDS'])
            hitrows = lookups[phase].get(locus, [])
            if hitrows:
                n_match += 1
                matched_genes.add(locus)
            base = {'phase': 'exponential' if phase == 'exp' else 'stationary',
                    'prediction_module': module,
                    'prediction_id': a['#ID'],
                    'Rfam_name': a['Rfam_name'],
                    'Rfam_ID': a.get('Rfam_ID', ''),
                    'associated_CDS': a['Associated_CDS'],
                    'gene_locus': locus or '',
                    'chr': a.get('Genome', ''),
                    'element_start': a.get('Start_genome', ''),
                    'element_end': a.get('End_genome', ''),
                    'strand': a.get('Strand', ''),
                    'E_value': a['E_value'],
                    'same_phase_DEG': 'yes' if hitrows else 'no',
                    'deg_record_count': len(hitrows),
                    'matched_gene_symbols': ';'.join(sorted({str(h['gene_symbol']) for h in hitrows if h['gene_symbol']})),
                    'matched_log2FC': ';'.join(str(h['log2FC']) for h in hitrows if str(h['log2FC']).strip()),
                    'matched_p_values': ';'.join(str(h['p_value']) for h in hitrows if str(h['p_value']).strip()),
                    'matched_adjusted_p': ';'.join(str(h['adjusted_p']) for h in hitrows if str(h['adjusted_p']).strip())}
            all_records.append(base)
            if hitrows:
                matches.append(base)
        summaries.append({'phase': phase, 'prediction_module': module,
                          'total_predictions': len(ann),
                          'predictions_associated_with_same_phase_DEG': n_match,
                          'unique_matched_DEG_loci': len(matched_genes),
                          'total_curated_DEG_rows': totals[phase]})
        print(f'{phase} {module}: {len(ann)} total, {n_match} DEG-associated predictions, '
              f'{len(matched_genes)} distinct DEG loci')

    fields = list(all_records[0])
    write_tsv(OUTPUT / 'all_predictions_with_DEG_flags.tsv', all_records, fields)
    write_tsv(OUTPUT / 'same_phase_DEG_associations.tsv', matches, fields)
    write_tsv(OUTPUT / 'deg_mapping_audit.tsv', audits,
              ['phase', 'deg_row', 'gene_symbol', 'resolved_locus', 'resolution', 'possible_loci'])
    write_tsv(OUTPUT / 'summary.tsv', summaries,
              ['phase', 'prediction_module', 'total_predictions',
               'predictions_associated_with_same_phase_DEG',
               'unique_matched_DEG_loci', 'total_curated_DEG_rows'])
    print(f'\nResults saved in {OUTPUT}')
    print('IMPORTANT: These are ANNOgesic Associated_CDS matches, NOT coordinate overlaps.')
    print('Review deg_mapping_audit.tsv for unresolved or ambiguous DEG rows.')

if __name__ == '__main__':
    main()
