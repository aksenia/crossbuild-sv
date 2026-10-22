#!/usr/bin/env python3
"""
Pool per-sample comparison metrics (compare_native.py metrics.tsv) into one
cohort table of the same long shape, plus n_samples.

Counts are summed per (tool, set, stratum, metric) and rates are recomputed
from the pooled counts (not averaged), with the formulas of compare_native.py.
median_abs_* are reported as the median of the per-sample medians.
n_samples = samples that have this metric row.

Usage: aggregate.py --output cohort.metrics.tsv runs/*/comparison/*.metrics.tsv
"""
import argparse
import csv
import statistics
import sys
from collections import defaultdict

NA = '.'
# rate: (numerator, denominator), as in compare_native.py
RATES = {
    'precision': ('matched_ids', 'lifted_ids'),
    'recall': ('native_matched', 'native_ids'),
    'exact_rate': ('exact_ids', 'matched_ids'),
    'within_ci_rate': ('within_ci_ids', 'matched_ids'),
    'gt_match_rate': ('gt_match_ids', 'matched_ids'),
    'yield': ('lifted_ids', 'source_ids'),
    'all_tools_agree_rate': ('all_tools_agree_ids', 'lifted_ids'),
}


def rate(num, den):
    return round(num / den, 4) if den else NA


def aggregate(paths):
    """{(tool, set, stratum): {metric: [values]}}, key order as first seen."""
    groups = defaultdict(lambda: defaultdict(list))
    samples = set()
    for path in paths:
        with open(path) as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        file_samples = {r['sample_id'] for r in rows}
        if file_samples & samples:   # would silently double the counts
            sys.exit(f'{path}: sample_id already seen: {", ".join(sorted(file_samples & samples))}')
        samples |= file_samples
        for r in rows:
            groups[(r['tool'], r['set'], r['stratum'])][r['metric']].append(r['value'])
    return groups


def pooled(metrics):
    """Yield (metric, value, n_samples) for one (tool, set, stratum)."""
    sums = {m: sum(int(v) for v in vs if v != NA) for m, vs in metrics.items()
            if m not in RATES and not m.startswith('median_')}
    for m, vs in metrics.items():
        if m in RATES:
            num, den = RATES[m]
            value = rate(sums.get(num, 0), sums.get(den, 0))
        elif m.startswith('median_'):
            vals = [float(v) for v in vs if v != NA]
            value = statistics.median(vals) if vals else NA
        else:
            value = sums[m]
        yield m, value, len(vs)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('metrics', nargs='+', help='per-sample metrics.tsv files')
    p.add_argument('--cohort', default='cohort', help='sample_id written to the output')
    p.add_argument('--output', required=True)
    a = p.parse_args()

    with open(a.output, 'w') as fh:
        fh.write('sample_id\ttool\tset\tstratum\tmetric\tvalue\tn_samples\n')
        for key, metrics in aggregate(a.metrics).items():
            for m, value, n in pooled(metrics):
                fh.write('\t'.join(map(str, (a.cohort, *key, m, value, n))) + '\n')


if __name__ == '__main__':
    main()
