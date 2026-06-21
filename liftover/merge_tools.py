#!/usr/bin/env python3
"""
Collapse sorted VCF records with identical CHROM/POS/ID, joining their
LIFTOVER_TOOL INFO values into a comma-separated list.

For SVs the merge key is CHROM/POS/ID (the ID was stamped as
CHROM_POS_SVTYPE_END in hg19 coordinates before liftover).
Records at the same hg38 position but with different IDs — i.e. from
different hg19 variants — are kept separate.

Also:
- Carries tool-specific INFO fields (e.g. transanno MULTIMAP/REF_CHANGED).
- Optionally annotates records with Picard rejection metadata when
  --picard-rejected is supplied.

Reads from stdin (sorted VCF), writes to stdout.
"""
import sys
import argparse

_TOOL_FIELDS = {
    'transanno': ('MULTIMAP', 'REF_CHANGED'),
}

_PICARD_HEADERS = [
    '##INFO=<ID=PICARD_REJECT_REASON,Number=1,Type=String,Description="Picard LiftoverVcf rejection reason">',
    '##INFO=<ID=PICARD_ATTEMPTED_LOCUS,Number=1,Type=String,Description="hg38 locus Picard attempted before failing (MismatchedRefAllele only)">',
    '##INFO=<ID=PICARD_ATTEMPTED_ALLELES,Number=1,Type=String,Description="hg38 alleles Picard attempted before failing (MismatchedRefAllele only)">',
]


def parse_info(info_str):
    d = {}
    for field in info_str.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            d[k] = v
        elif field:
            d[field] = True
    return d


def format_info(d):
    parts = []
    for k, v in d.items():
        if v is True:
            parts.append(k)
        else:
            parts.append(f'{k}={v}')
    return ';'.join(parts)


def load_picard_rejected(path):
    """Parse Picard rejected VCF → {variant_id: {reason, locus, alleles}}."""
    rejected = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            vid    = parts[2]   # ID stamped as CHROM_POS_SVTYPE_END in hg19
            reason = parts[6]
            info   = parse_info(parts[7])
            rejected[vid] = {
                'reason':  reason,
                'locus':   info.get('AttemptedLocus', ''),
                'alleles': info.get('AttemptedAlleles', ''),
            }
    return rejected


def _apply_picard(variant_id, info, picard_rejected):
    if not picard_rejected or variant_id not in picard_rejected:
        return
    rec = picard_rejected[variant_id]
    info['PICARD_REJECT_REASON'] = rec['reason']
    if rec['locus']:
        info['PICARD_ATTEMPTED_LOCUS'] = rec['locus']
    if rec['alleles']:
        info['PICARD_ATTEMPTED_ALLELES'] = rec['alleles']


def flush_group(records, picard_rejected):
    if not records:
        return

    merged = list(records[0])

    if len(records) == 1:
        info = parse_info(merged[7])
        _apply_picard(merged[2], info, picard_rejected)
        merged[7] = format_info(info)
        sys.stdout.write('\t'.join(merged) + '\n')
        return

    seen = {}
    info = parse_info(merged[7])
    for i, parts in enumerate(records):
        rec_info = info if i == 0 else parse_info(parts[7])
        for t in rec_info.get('LIFTOVER_TOOL', '').split(','):
            t = t.strip()
            if t:
                seen[t] = True
        tool = rec_info.get('LIFTOVER_TOOL', '')
        for tool_name, fields in _TOOL_FIELDS.items():
            if tool_name in tool:
                for field in fields:
                    if field in rec_info:
                        info[field] = rec_info[field]

    info['LIFTOVER_TOOL'] = ','.join(seen.keys())
    _apply_picard(merged[2], info, picard_rejected)
    merged[7] = format_info(info)
    sys.stdout.write('\t'.join(merged) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--picard-rejected', metavar='VCF',
                        help='Picard rejected VCF; annotates records with rejection reason')
    args = parser.parse_args()

    picard_rejected = load_picard_rejected(args.picard_rejected) if args.picard_rejected else {}
    picard_headers_injected = not picard_rejected

    prev_key = None
    group = []

    for raw_line in sys.stdin:
        line = raw_line.rstrip('\n')
        if line.startswith('#'):
            if not picard_headers_injected and line.startswith('#CHROM'):
                for h in _PICARD_HEADERS:
                    sys.stdout.write(h + '\n')
                picard_headers_injected = True
            sys.stdout.write(line + '\n')
            continue

        parts = line.split('\t')
        if len(parts) < 8:
            sys.stdout.write(line + '\n')
            continue

        # Merge key: CHROM / POS / ID  (ID = CHROM_POS_SVTYPE_END from hg19)
        key = (parts[0], parts[1], parts[2])
        if key != prev_key:
            flush_group(group, picard_rejected)
            group = [parts]
            prev_key = key
        else:
            group.append(parts)

    flush_group(group, picard_rejected)


if __name__ == '__main__':
    main()
