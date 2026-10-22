#!/usr/bin/env python3
"""
Compare lifted SV calls with native hg38 calls of the same individual.

SVDB (run per tool and set, native first) only decides which native and
lifted calls belong together. Everything else is computed here from the
coordinates of the normalized VCFs, and written as two flat files keyed by
sample_id so a cohort is just their concatenation:

  links.tsv.gz   one row per (tool, set, link): matched native+lifted pair,
                 native_only call, or lifted_only placement. tool=union is the
                 merged set of all tools. Columns: see docs/comparison.md.
  metrics.tsv    long format (sample_id, tool, set, stratum, metric, value):
                 counts and rates per stratum (all, svtype=, size=,
                 svtype=|size=, region=, precision=). Counts are kept next
                 to rates so a cohort can pool numerators/denominators
                 instead of averaging rates.
"""
import argparse
import bisect
import gzip
import re
import statistics
from collections import Counter, defaultdict

from finalize_sv import ChainIndex, bare_contig   # sibling script

SRC_FLAGS = ('CUPs', 'DISCREPs', 'CENTEL', 'SEGDUP')
HG38_FLAGS = ('CUPs', 'DISCREPs')
INTERVAL_TYPES = {'DEL', 'DUP', 'INV'}
SETS = ('main', 'refchanged')
SIZE_BINS = [(50, '<50'), (100, '50-100'), (500, '100-500'), (1000, '500-1k'),
             (10000, '1k-10k'), (float('inf'), '>=10k')]
NA = '.'
# Unmatched calls: what the other build has within NEAR_WINDOW bp, first
# matching reason wins (see OtherBuild.explain).
NEAR_WINDOW = 1000
REASONS = ('near_same_type', 'filtered_other_build', 'other_type_other_build',
           'no_hg19_sequence', 'none')
OFFSET_BINS = [(0, '0'), (10, '1-10'), (100, '11-100'), (1000, '101-1k'), (float('inf'), '>1k')]

COLUMNS = [
    'sample_id', 'tool', 'set', 'match_class', 'native_id', 'source_id', 'n_links',
    'nat_chrom', 'nat_pos', 'nat_end', 'nat_svtype', 'nat_svlen', 'nat_gt', 'nat_precise',
    'nat_cipos', 'nat_ciend', *(f'nat_{f}' for f in HG38_FLAGS),
    'src_chrom', 'src_pos', 'src_end', 'src_svtype', 'src_svlen', 'src_gt',
    *(f'src_{f}' for f in SRC_FLAGS),
    'lift_chrom', 'lift_pos', 'lift_end', 'lift_svlen', 'tools_same_placement',
    'n_tools_lifted', 'n_placements', 'tools_rejected', *(f'lift_{f}' for f in HG38_FLAGS),
    'exact', 'd_pos', 'd_end', 'd_svlen', 'recip_overlap', 'within_ci',
    'chrom_match', 'svtype_match', 'gt_match', 'closest_tool',
    'diff_reason', 'nearest_dist', 'nearest_size_ratio',
    'svtype', 'size_bin', 'in_region', 'precision',
]


# ── VCF reading ──────────────────────────────────────────────────────────────
def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def parse_info(info_str):
    d = {}
    for field in info_str.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            d[k] = v
        elif field and field != '.':
            d[field] = True
    return d


def records(path):
    """Yield dicts for each record: core columns, parsed INFO, first-sample GT."""
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            gt = NA
            if len(p) > 9:
                keys, vals = p[8].split(':'), p[9].split(':')
                gt = vals[keys.index('GT')] if 'GT' in keys and keys.index('GT') < len(vals) else NA
            yield {'chrom': p[0], 'pos': int(p[1]), 'id': p[2], 'ref': p[3], 'alt': p[4],
                   'filter': p[6], 'info': parse_info(p[7]), 'info_str': p[7], 'gt': gt, 'last': p[-1]}


def end_of(r):
    try:
        return int(r['info']['END'])
    except (KeyError, ValueError):
        return r['pos']


def svlen_of(r):
    """Event size in bp (positive) from the record itself."""
    alt = r['alt'].split(',')[0]
    if not alt.startswith('<'):
        return abs(len(alt) - len(r['ref'])) or len(r['ref']) - 1
    if r['info'].get('SVTYPE') == 'INS':
        try:
            return abs(int(str(r['info']['SVLEN']).split(',')[0]))
        except (KeyError, ValueError):
            return 0
    return end_of(r) - r['pos']


def norm_gt(gt):
    return '/'.join(sorted(gt.replace('|', '/').split('/'))) if gt != NA else NA


def ci(r, key):
    try:
        lo, hi = (int(x) for x in r['info'][key].split(','))
        return lo, hi
    except (KeyError, ValueError):
        return 0, 0


def flags(r, names):
    return {f: int(f in r['info']) for f in names}


# ── Inputs ───────────────────────────────────────────────────────────────────
def kv_args(items):
    """['crossmap=path', ...] -> {'crossmap': 'path'}"""
    return dict(item.split('=', 1) for item in items or [])


def load_rejections(paths, qc_paths):
    """{source_id: {tool: reason}} from tool-native and post-liftover QC rejects."""
    rej = defaultdict(dict)
    for tool, path in paths.items():
        for r in records(path):
            # transanno appends FAILED_REASON without a ';' separator
            failed = re.search(r'FAILED_REASON=([^;]+)', r['info_str'])
            reason = ((failed and failed.group(1))              # transanno
                      or (r['last'] if r['last'].startswith('Fail') else None)  # CrossMap .unmap
                      or r['filter'])                           # bcftools, Picard
            rej[r['id']][tool] = str(reason).replace(',', '|')   # tools_rejected is comma-joined
    for tool, path in qc_paths.items():
        for r in records(path):
            rej[r['id']].setdefault(tool, 'qc:' + str(r['info'].get('QC_REJECT_REASON', NA)))
    return rej


def load_placements(tool_vcfs):
    """{source_id: {(chrom, pos, end): [tools]}} and per-tool {(id, chrom, pos): svlen}."""
    placements = defaultdict(lambda: defaultdict(list))
    svlens = {}
    for tool, path in tool_vcfs.items():
        for r in records(path):
            key = (r['chrom'], r['pos'], end_of(r))
            placements[r['id']][key].append(tool)
            svlens[(r['id'], r['chrom'], r['pos'])] = svlen_of(r)
    return placements, svlens


def svdb_clusters(path):
    """Yield {'native': [(id, chrom, pos)], 'lifted': [...]} for each SVDB record."""
    for r in records(path):
        members = {}
        for tag in ('native', 'lifted'):
            pos_f, chrom_f = r['info'].get(f'{tag}_POS'), r['info'].get(f'{tag}_CHROM')
            if pos_f and chrom_f:
                chroms = dict(x.rsplit('|', 1) for x in chrom_f.split(','))
                members[tag] = [(vid, chroms[vid], int(p))
                                for vid, p in (x.rsplit('|', 1) for x in pos_f.split(','))]
            elif r['info'].get('set') == tag:          # singleton without member fields
                members[tag] = [(r['id'], r['chrom'], r['pos'])]
            else:
                members[tag] = []
        yield members


def offset_bin(d):
    return next(label for limit, label in OFFSET_BINS if d <= limit)


# ── Other build: why is a call unmatched? ────────────────────────────────────
class OtherBuild:
    """
    Raw calls of the other build (any FILTER, any SVTYPE incl. BND), on hg38
    coordinates. A call missing from the comparison may still be there as a
    filtered call, a BND, or a nearby call SVDB did not link.
    """

    def __init__(self, calls):
        """calls: iterable of (chrom, pos, filter, svtype, svlen)."""
        self.calls = defaultdict(list)
        for chrom, *rest in calls:
            self.calls[bare_contig(chrom)].append(tuple(rest))
        for v in self.calls.values():
            v.sort()
        self.starts = {c: [x[0] for x in v] for c, v in self.calls.items()}

    def near(self, chrom, pos):
        c = bare_contig(chrom)
        starts = self.starts.get(c, [])
        return self.calls[c][bisect.bisect_left(starts, pos - NEAR_WINDOW):
                             bisect.bisect_right(starts, pos + NEAR_WINDOW)] if starts else []

    def explain(self, chrom, pos, svtype, svlen):
        """(reason, nearest_dist, nearest_size_ratio); the last two only for near_same_type."""
        hits = self.near(chrom, pos)
        same = [h for h in hits if h[2] == svtype and h[1] in ('PASS', '.')]
        if same:
            p, _, _, n = min(same, key=lambda h: abs(h[0] - pos))
            ratio = round(min(svlen, n) / max(svlen, n), 3) if svlen and n else NA
            return 'near_same_type', abs(p - pos), ratio
        if any(h[2] == svtype for h in hits):
            return 'filtered_other_build', NA, NA
        if hits:
            return 'other_type_other_build', NA, NA
        return 'none', NA, NA


def raw_calls(path, chains=None):
    """(chrom, pos, filter, svtype, svlen) per record; hg19 records lifted by POS through chains."""
    for r in records(path):
        chrom, pos = r['chrom'], r['pos']
        if chains:
            hit = chains.lift_point(chrom, pos)
            if not hit:
                continue
            chrom, pos = hit
        yield chrom, pos, r['filter'], r['info'].get('SVTYPE', NA), svlen_of(r)


class Hg19Coverage:
    """hg38 positions some chain block maps from hg19 (else the locus has no hg19 sequence)."""

    def __init__(self, chains):
        spans = defaultdict(list)
        for cs in chains.by_source.values():
            for c in cs:
                for b in c.blocks:
                    n = b.source_end - b.source_start
                    s = b.target_start if c.target_strand == '+' else c.target_size - b.target_start - n
                    spans[bare_contig(c.target_name)].append((s + 1, s + n))
        self.spans = {}
        for chrom, v in spans.items():
            merged = []
            for s, e in sorted(v):
                if merged and s <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self.spans[chrom] = ([s for s, _ in merged], [e for _, e in merged])

    def __contains__(self, chrom_pos):
        chrom, pos = chrom_pos
        starts, ends = self.spans.get(bare_contig(chrom), ([], []))
        i = bisect.bisect_right(starts, pos) - 1
        return i >= 0 and pos <= ends[i]


# ── Rows ─────────────────────────────────────────────────────────────────────
def size_bin(n):
    return next(label for limit, label in SIZE_BINS if n < limit)


def recip_overlap(svtype, a_pos, a_end, b_pos, b_end):
    if svtype not in INTERVAL_TYPES:
        return NA
    ov = min(a_end, b_end) - max(a_pos, b_pos)
    la, lb = max(a_end - a_pos, 1), max(b_end - b_pos, 1)
    return round(max(ov, 0) / max(la, lb), 3)


class Builder:
    def __init__(self, a, tools):
        self.sample, self.tools = a.sample, tools
        self.src = {r['id']: r for r in records(a.source)}
        self.nat = {r['id']: r for r in records(a.native)}
        self.placements, self.svlens = load_placements(kv_args(a.tool_vcf))
        self.rejected = load_rejections(kv_args(a.rejected), kv_args(a.qc_rejected))
        self.lift_flags = {}
        for path in a.lifted_annotated:
            for r in records(path):
                self.lift_flags[(r['id'], r['chrom'], r['pos'], end_of(r))] = flags(r, HG38_FLAGS)
        chains = ChainIndex(a.chain)
        self.hg19 = OtherBuild(raw_calls(a.raw_source, chains))
        self.hg38 = OtherBuild(raw_calls(a.raw_native))
        self.hg19_coverage = Hg19Coverage(chains)

    def placement(self, sid, chrom, pos, tool):
        """The placement SVDB saw; prefer this tool's own if several share CHROM/POS."""
        keys = [(k, ts) for k, ts in self.placements.get(sid, {}).items() if k[:2] == (chrom, pos)]
        keys.sort(key=lambda kt: tool not in kt[1])
        return keys[0] if keys else ((chrom, pos, pos), [])

    def closest_tool(self, sid, nat):
        dist = {}
        for (chrom, pos, end), tools in self.placements.get(sid, {}).items():
            if chrom == nat['chrom']:
                for t in tools:
                    dist[t] = abs(pos - nat['pos']) + abs(end - end_of(nat))
        if not dist:
            return NA
        best = min(dist.values())
        return ','.join(t for t in self.tools if dist.get(t) == best)

    def row(self, tool, set_, cls, nid, lifted):
        row = dict.fromkeys(COLUMNS, NA)
        row.update(sample_id=self.sample, tool=tool, set=set_, match_class=cls)
        nat = self.nat.get(nid) if nid else None
        regions = []

        if nat:
            nlen = svlen_of(nat)
            row.update(native_id=nid, nat_chrom=nat['chrom'], nat_pos=nat['pos'], nat_end=end_of(nat),
                       nat_svtype=nat['info'].get('SVTYPE', NA), nat_svlen=nlen, nat_gt=nat['gt'],
                       nat_precise=int('IMPRECISE' not in nat['info']),
                       nat_cipos='{},{}'.format(*ci(nat, 'CIPOS')), nat_ciend='{},{}'.format(*ci(nat, 'CIEND')))
            nf = flags(nat, HG38_FLAGS)
            row.update({f'nat_{k}': v for k, v in nf.items()})
            regions += nf.values()

        if lifted:
            sid, chrom, pos = lifted
            src = self.src.get(sid)
            if src:
                sf = flags(src, SRC_FLAGS)
                row.update(source_id=sid, src_chrom=src['chrom'], src_pos=src['pos'], src_end=end_of(src),
                           src_svtype=src['info'].get('SVTYPE', NA), src_svlen=svlen_of(src), src_gt=src['gt'],
                           **{f'src_{k}': v for k, v in sf.items()})
                regions += sf.values()
            (chrom, pos, end), here = self.placement(sid, chrom, pos, tool)
            lf = self.lift_flags.get((sid, chrom, pos, end), dict.fromkeys(HG38_FLAGS, 0))
            regions += lf.values()
            lifted_tools = {t for ts in self.placements.get(sid, {}).values() for t in ts}
            row.update(lift_chrom=chrom, lift_pos=pos, lift_end=end,
                       lift_svlen=self.svlens.get((sid, chrom, pos), NA),
                       tools_same_placement=','.join(t for t in self.tools if t in here) or NA,
                       n_tools_lifted=len(lifted_tools), n_placements=len(self.placements.get(sid, {})),
                       tools_rejected=','.join(f'{t}:{r}' for t, r in sorted(self.rejected.get(sid, {}).items())
                                               if t not in lifted_tools) or NA,
                       **{f'lift_{k}': v for k, v in lf.items()})

        if nat and lifted:
            same = row['lift_chrom'] == nat['chrom']
            row.update(chrom_match=int(same),
                       svtype_match=int(row['src_svtype'] == row['nat_svtype']),
                       gt_match=int(norm_gt(row['src_gt']) == norm_gt(nat['gt'])),
                       closest_tool=self.closest_tool(row['source_id'], nat))
            if same:
                d_pos, d_end = row['lift_pos'] - nat['pos'], row['lift_end'] - end_of(nat)
                (plo, phi), (elo, ehi) = ci(nat, 'CIPOS'), ci(nat, 'CIEND')
                row.update(exact=int(d_pos == 0 and d_end == 0), d_pos=d_pos, d_end=d_end,
                           within_ci=int(plo <= d_pos <= phi and elo <= d_end <= ehi),
                           recip_overlap=recip_overlap(row['nat_svtype'], row['lift_pos'], row['lift_end'],
                                                       nat['pos'], end_of(nat)))
            else:
                row.update(exact=0, within_ci=0)
            if row['lift_svlen'] != NA:
                row['d_svlen'] = row['lift_svlen'] - row['nat_svlen']

        # Strata come from the native call when there is one (the reference),
        # else from the source call.
        ref = nat or self.src.get(row['source_id'])
        if ref:
            row['svtype'] = ref['info'].get('SVTYPE', NA)
            row['size_bin'] = size_bin(svlen_of(ref))
            row['precision'] = 'imprecise' if 'IMPRECISE' in ref['info'] else 'precise'
        row['in_region'] = int(any(regions))

        if cls == 'native_only':
            why = self.hg19.explain(nat['chrom'], nat['pos'], row['nat_svtype'], row['nat_svlen'])
            if why[0] == 'none' and (nat['chrom'], nat['pos']) not in self.hg19_coverage:
                why = ('no_hg19_sequence', NA, NA)
        elif cls == 'lifted_only':
            svlen = row['lift_svlen'] if row['lift_svlen'] != NA else row['src_svlen']
            why = self.hg38.explain(row['lift_chrom'], row['lift_pos'], row['src_svtype'], svlen)
        if cls != 'matched':
            row.update(zip(('diff_reason', 'nearest_dist', 'nearest_size_ratio'), why))
        return row

    def links(self, tool, set_, path):
        rows = []
        for m in svdb_clusters(path):
            if m['native'] and m['lifted']:
                rows += [self.row(tool, set_, 'matched', n[0], l) for n in m['native'] for l in m['lifted']]
            elif m['lifted']:
                rows += [self.row(tool, set_, 'lifted_only', None, l) for l in m['lifted']]
            elif set_ == 'main':   # REF-changed set: native-only is not defined
                rows += [self.row(tool, set_, 'native_only', n[0], None) for n in m['native']]
        return self.one_to_one(rows)

    def union(self, rows):
        """
        tool=union rows from per-tool rows: every distinct matched placement;
        lifted placements of source IDs no tool matched; native calls no tool
        matched. (One SVDB run over the merged set would let a native call join
        only one of a source ID's placements.)
        """
        out = []
        for set_ in SETS:
            rs = [r for r in rows if r['set'] == set_]
            matched_src = {r['source_id'] for r in rs if r['match_class'] == 'matched'}
            matched_nat = {r['native_id'] for r in rs if r['match_class'] == 'matched'}
            seen = set()
            for r in rs:
                cls = r['match_class']
                if cls == 'lifted_only' and r['source_id'] in matched_src:
                    continue
                if cls == 'native_only' and r['native_id'] in matched_nat:
                    continue
                key = (cls, r['native_id'], r['source_id'], r['lift_chrom'], r['lift_pos'], r['lift_end'])
                if key not in seen:
                    seen.add(key)
                    out.append(dict(r, tool='union'))
        return self.one_to_one(out)

    def one_to_one(self, rows):
        """
        Keep one source ID per native call. SVDB's greedy clustering can put
        several lifted calls into one native cluster (seen: near-duplicate
        source calls a few bp apart); the closest stays matched, the others
        become lifted_only unless matched to another native call. n_links
        keeps the count before this step, so n_links > 1 marks such clusters.
        """
        self.count_links(rows)
        by_native = defaultdict(list)
        for r in rows:
            if r['match_class'] == 'matched':
                by_native[r['native_id']].append(r)
        demoted = set()
        for rs in by_native.values():
            if len({r['source_id'] for r in rs}) > 1:
                best = min(rs, key=self.link_distance)['source_id']
                demoted |= {id(r) for r in rs if r['source_id'] != best}
        kept = [r for r in rows if id(r) not in demoted]
        still_matched = {r['source_id'] for r in kept if r['match_class'] == 'matched'}
        seen = set()
        for r in rows:
            key = (r['source_id'], r['lift_chrom'], r['lift_pos'])
            if id(r) in demoted and r['source_id'] not in still_matched and key not in seen:
                seen.add(key)
                kept.append(self.row(r['tool'], r['set'], 'lifted_only', None, key))
        return kept

    @staticmethod
    def link_distance(r):
        """Sort key: same chrom first, then |d_pos| + |d_end|; source_id breaks ties."""
        if r['d_pos'] == NA:
            return (1, 0, r['source_id'])
        return (0, abs(r['d_pos']) + abs(r['d_end']), r['source_id'])

    @staticmethod
    def count_links(rows):
        per_native = defaultdict(set)
        for r in rows:
            if r['match_class'] == 'matched':
                per_native[r['native_id']].add(r['source_id'])
        for r in rows:
            if r['match_class'] == 'matched':
                r['n_links'] = len(per_native[r['native_id']])
        return rows


# ── Metrics ──────────────────────────────────────────────────────────────────
def rate(num, den):
    return round(num / den, 4) if den else NA


def strata(row):
    yield 'all'
    for col, name in (('svtype', 'svtype'), ('size_bin', 'size'), ('precision', 'precision')):
        if row[col] != NA:
            yield f'{name}={row[col]}'
    # crossed stratum: size bins mean different things per SV type
    if row['svtype'] != NA and row['size_bin'] != NA:
        yield f"svtype={row['svtype']}|size={row['size_bin']}"
    yield 'region=any' if row['in_region'] else 'region=none'


def offset_metrics(rs, matched):
    """
    Matched IDs by their smallest max(|d_pos|, |d_end|); unmatched calls by
    diff_reason (first in REASONS order if an ID has several rows), and
    near_same_type ones by distance to the nearest same-type call. Large
    offsets mostly show up there: SVDB links only close calls.
    """
    m = {}
    best = {}
    for r in matched:   # other chrom ranks after any offset
        off = float('inf') if r['d_pos'] == NA else max(abs(r['d_pos']), abs(r['d_end']))
        best[r['source_id']] = min(best.get(r['source_id'], off), off)
    counts = Counter('other_chrom' if o == float('inf') else offset_bin(o) for o in best.values())
    for label in [label for _, label in OFFSET_BINS] + ['other_chrom']:
        m[f'matched_offset_{label}'] = counts[label]

    for cls, key in (('native_only', 'native_id'), ('lifted_only', 'source_id')):
        linked = {r[key] for r in matched}
        per_id = {}
        for r in rs:
            if r['match_class'] == cls and r[key] not in linked:
                why = (REASONS.index(r['diff_reason']), r['nearest_dist'] if r['nearest_dist'] != NA else 0)
                per_id[r[key]] = min(per_id.get(r[key], why), why)
        reasons = Counter(REASONS[i] for i, _ in per_id.values())
        for reason in REASONS:
            if cls == 'native_only' or reason != 'no_hg19_sequence':
                m[f'{cls}_{reason}'] = reasons[reason]
        near = Counter(offset_bin(d) for i, d in per_id.values() if REASONS[i] == 'near_same_type')
        for _, label in OFFSET_BINS[:-1]:   # nothing beyond NEAR_WINDOW
            m[f'{cls}_near_offset_{label}'] = near[label]
    return m


def comparison_metrics(rows):
    groups = defaultdict(list)
    for r in rows:
        for s in strata(r):
            groups[(r['tool'], r['set'], s)].append(r)
    out = []
    for (tool, set_, stratum), rs in sorted(groups.items()):
        matched = [r for r in rs if r['match_class'] == 'matched']
        ids = lambda cls: {r['source_id'] for r in rs if r['match_class'] in cls}
        m = {
            'lifted_ids': len(ids(('matched', 'lifted_only'))),
            'matched_ids': len(ids(('matched',))),
            'lifted_only_ids': len(ids(('lifted_only',)) - ids(('matched',))),
            'native_ids': len({r['native_id'] for r in rs if r['match_class'] in ('matched', 'native_only')}),
            'native_matched': len({r['native_id'] for r in matched}),
            'exact_ids': len({r['source_id'] for r in matched if r['exact'] == 1}),
            'within_ci_ids': len({r['source_id'] for r in matched if r['within_ci'] == 1}),
            'gt_match_ids': len({r['source_id'] for r in matched if r['gt_match'] == 1}),
        }
        m['native_only'] = m['native_ids'] - m['native_matched']
        m['precision'] = rate(m['matched_ids'], m['lifted_ids'])
        m['recall'] = rate(m['native_matched'], m['native_ids'])
        if set_ == 'refchanged':   # no native-only calls here, so no recall
            for k in ('native_ids', 'native_only', 'recall'):
                del m[k]
        m['exact_rate'] = rate(m['exact_ids'], m['matched_ids'])
        m['within_ci_rate'] = rate(m['within_ci_ids'], m['matched_ids'])
        m['gt_match_rate'] = rate(m['gt_match_ids'], m['matched_ids'])
        m.update({k: v for k, v in offset_metrics(rs, matched).items()
                  if set_ == 'main' or not k.startswith('native_only')})
        for d in ('d_pos', 'd_end'):
            vals = [abs(r[d]) for r in matched if r[d] != NA]
            m[f'median_abs_{d}'] = statistics.median(vals) if vals else NA
        out += [(tool, set_, stratum, k, v) for k, v in m.items()]
    return out


def yield_metrics(b, refchanged):
    """Per-tool liftover accounting over all source IDs (set=all)."""
    out = []
    by_type = defaultdict(set)
    for sid, r in b.src.items():
        by_type['all'].add(sid)
        by_type[f"svtype={r['info'].get('SVTYPE', NA)}"].add(sid)
    for stratum, sids in sorted(by_type.items()):
        for tool in b.tools:
            lifted = {s for s in sids if any(tool in ts for ts in b.placements.get(s, {}).values())}
            rej = {s for s in sids if tool in b.rejected.get(s, {}) and s not in lifted}
            qc = {s for s in rej if b.rejected[s][tool].startswith('qc:')}
            m = {'source_ids': len(sids), 'lifted_ids': len(lifted),
                 'tool_rejected_ids': len(rej - qc), 'qc_rejected_ids': len(qc),
                 'unaccounted_ids': len(sids - lifted - rej),
                 'refchanged_flagged_ids': len({s for s in sids if tool in refchanged.get(s, ())}),
                 'yield': rate(len(lifted), len(sids))}
            out += [(tool, 'all', stratum, k, v) for k, v in m.items()]
        any_tool = {s for s in sids if b.placements.get(s)}
        agree = {s for s in any_tool if len(b.placements[s]) == 1
                 and len(next(iter(b.placements[s].values()))) == len(b.tools)}
        m = {'source_ids': len(sids), 'lifted_ids': len(any_tool),
             'all_tools_agree_ids': len(agree),
             'multi_placement_ids': len({s for s in any_tool if len(b.placements[s]) > 1}),
             'refchanged_ids': len(sids & set(refchanged)),
             'yield': rate(len(any_tool), len(sids)),
             'all_tools_agree_rate': rate(len(agree), len(any_tool))}
        out += [('union', 'all', stratum, k, v) for k, v in m.items()]
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--sample', required=True)
    p.add_argument('--tools', required=True, help='comma-separated tool order')
    p.add_argument('--source', required=True, help='annotated hg19 source VCF')
    p.add_argument('--native', required=True, help='normalized, annotated native hg38 VCF')
    p.add_argument('--tool-vcf', nargs='+', required=True, metavar='TOOL=VCF')
    p.add_argument('--rejected', nargs='+', required=True, metavar='TOOL=VCF')
    p.add_argument('--qc-rejected', nargs='+', required=True, metavar='TOOL=VCF')
    p.add_argument('--lifted-annotated', nargs='+', required=True, help='annotated merged hg38 VCFs')
    p.add_argument('--refchanged-ids', required=True)
    p.add_argument('--svdb', nargs='+', required=True, metavar='TOOL:SET=VCF')
    p.add_argument('--raw-source', required=True, help='unfiltered hg19 input VCF (diff_reason)')
    p.add_argument('--raw-native', required=True, help='unfiltered native hg38 VCF (diff_reason)')
    p.add_argument('--chain', required=True, help='hg19->hg38 UCSC chain (diff_reason)')
    p.add_argument('--links', required=True)
    p.add_argument('--metrics', required=True)
    a = p.parse_args()

    b = Builder(a, a.tools.split(','))
    refchanged = defaultdict(set)
    with open(a.refchanged_ids) as fh:
        for line in fh:
            sid, tool = line.split()
            refchanged[sid].add(tool)

    rows = []
    for key, path in kv_args(a.svdb).items():
        tool, set_ = key.split(':')
        rows += b.links(tool, set_, path)
    rows += b.union(rows)

    with gzip.open(a.links, 'wt') as fh:
        fh.write('\t'.join(COLUMNS) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in COLUMNS) + '\n')

    with open(a.metrics, 'w') as fh:
        fh.write('sample_id\ttool\tset\tstratum\tmetric\tvalue\n')
        for rec in yield_metrics(b, refchanged) + comparison_metrics(rows):
            fh.write('\t'.join(map(str, (a.sample, *rec))) + '\n')


if __name__ == '__main__':
    main()
