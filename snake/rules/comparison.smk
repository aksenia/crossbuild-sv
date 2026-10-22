import os

# ─── Compare lifted calls with native hg38 calls of the same individual ──────
# Native goes through the same filter, ID stamp and normalization as the
# lifted calls, so both sides share one representation. SVDB links native and
# lifted calls per tool and set (main / refchanged); compare_native.py
# computes all metrics from the coordinates and derives tool=union.

COMP_DIR   = os.path.join(RESULTS_DIR, "comparison")
SETS       = ["main", "refchanged"]

# Same contig subset as the source, in hg38 (chr) naming.
_NATIVE_TARGETS = ("-t " + ",".join(c if c.startswith("chr") else ("chrM" if c == "MT" else f"chr{c}")
                                    for c in map(str, _CONTIGS))) if _CONTIGS else ""

_REJECTED = {t: (f"{LIFTOVER_DIR}/{t}/{SAMPLE}.hg38.unmap.vcf" if t == "crossmap"
                 else f"{LIFTOVER_DIR}/{t}/{SAMPLE}.rejected.vcf") for t in TOOLS}


rule prep_native_vcf:
    """Filter, stamp IDs, normalize and region-annotate the native hg38 calls."""
    input:
        vcf = config["native_vcf"]
    output:
        vcf      = f"{COMP_DIR}/native/{SAMPLE}.hg38.vcf.gz",
        tbi      = f"{COMP_DIR}/native/{SAMPLE}.hg38.vcf.gz.tbi",
        rejected = f"{COMP_DIR}/native/{SAMPLE}.qc_rejected.vcf"
    log:
        f"{COMP_DIR}/native/{SAMPLE}.normalize.log"
    params:
        pass_flag = _PASS_FLAG,
        sv_filter = _SV_FILTER,
        targets   = _NATIVE_TARGETS,
        stamp_id  = _STAMP_ID,
        fasta     = config["ref"]["hg38_fasta"],
        script    = f"{SCRIPTS}/finalize_sv.py",
        cups      = config["regions"]["cups_hg38"],
        discreps  = config["regions"]["discreps_hg38"],
        header    = config["regions"]["header"],
        tmp       = f"{COMP_DIR}/native/{SAMPLE}.tmp.vcf"
    shell:
        """
        mkdir -p {COMP_DIR}/native
        HDR=$(mktemp)
        grep -E 'ID=(CUPs|DISCREPs),' {params.header} > $HDR
        bcftools view {params.pass_flag} {params.targets} -i '{params.sv_filter}' {input.vcf} \\
            | {params.stamp_id} \\
            | bcftools norm -f {params.fasta} -c w 2> {log} \\
            | python3 {params.script} --fasta {params.fasta} \\
                --output {params.tmp} --reject {output.rejected} 2>> {log}
        bcftools sort {params.tmp} \\
            | bcftools annotate -a {params.cups} -c CHROM,FROM,TO -m +CUPs -h $HDR \\
            | bcftools annotate -a {params.discreps} -c CHROM,FROM,TO -m +DISCREPs -Oz -o {output.vcf}
        tabix -f {output.vcf}
        rm -f {params.tmp} $HDR
        """


rule lifted_set:
    """One tool's normalized calls restricted to the main or REF-changed set."""
    input:
        vcf = f"{LIFTOVER_DIR}/{{tool}}/{SAMPLE}.hg38.vcf.gz",
        ids = rules.refchanged_ids.output.ids
    output:
        vcf = f"{COMP_DIR}/sets/{{tool}}.{{set}}.vcf.gz"
    wildcard_constraints:
        tool = "|".join(TOOLS),
        set  = "|".join(SETS)
    params:
        expr = lambda w: "-e" if w.set == "main" else "-i"
    shell:
        """
        mkdir -p {COMP_DIR}/sets
        IDS=$(mktemp)
        cut -f1 {input.ids} | sort -u > $IDS
        bcftools view {params.expr} "ID=@$IDS" -Oz -o {output.vcf} {input.vcf}
        rm -f $IDS
        """


rule svdb_link:
    """
    Link native and lifted calls. Native is listed first (it is the
    reference, so merged clusters carry native coordinates); --no_intra stops
    SVDB merging calls within one file.
    """
    input:
        native = rules.prep_native_vcf.output.vcf,
        lifted = rules.lifted_set.output.vcf
    output:
        vcf = f"{COMP_DIR}/svdb/{{tool}}.{{set}}.vcf"
    log:
        f"{COMP_DIR}/svdb/{{tool}}.{{set}}.log"
    params:
        overlap = config.get("comparison", {}).get("svdb_overlap", 0.8),
        profile = config.get("comparison", {}).get("svdb_profile", "sample")
    shell:
        """
        mkdir -p {COMP_DIR}/svdb
        svdb --merge --no_intra --overlap {params.overlap} --data_profile {params.profile} \\
            --vcf {input.native}:native {input.lifted}:lifted --priority native,lifted \\
            > {output.vcf} 2> {log}
        """


rule compare_native:
    """Build the per-sample links table and long-format metrics."""
    input:
        svdb      = expand(f"{COMP_DIR}/svdb/{{tool}}.{{set}}.vcf", tool=TOOLS, set=SETS),
        native    = rules.prep_native_vcf.output.vcf,
        source    = f"{RESULTS_DIR}/annotated/source/{SAMPLE}.hg19.annotated.vcf.gz",
        tool_vcfs = expand(f"{LIFTOVER_DIR}/{{tool}}/{SAMPLE}.hg38.vcf.gz", tool=TOOLS),
        qc        = expand(f"{LIFTOVER_DIR}/{{tool}}/{SAMPLE}.qc_rejected.vcf", tool=TOOLS),
        rejected  = list(_REJECTED.values()),
        annotated = expand(f"{RESULTS_DIR}/annotated/merged/{SAMPLE}.{{part}}.hg38.annotated.vcf.gz",
                           part=["merged", "refchanged"]),
        ids       = rules.refchanged_ids.output.ids
    output:
        links   = f"{COMP_DIR}/{SAMPLE}.links.tsv.gz",
        metrics = f"{COMP_DIR}/{SAMPLE}.metrics.tsv"
    params:
        script    = f"{SCRIPTS}/compare_native.py",
        tools     = ",".join(TOOLS),
        tool_vcfs = " ".join(f"{t}={LIFTOVER_DIR}/{t}/{SAMPLE}.hg38.vcf.gz" for t in TOOLS),
        qc        = " ".join(f"{t}={LIFTOVER_DIR}/{t}/{SAMPLE}.qc_rejected.vcf" for t in TOOLS),
        rejected  = " ".join(f"{t}={p}" for t, p in _REJECTED.items()),
        svdb      = " ".join(f"{t}:{s}={COMP_DIR}/svdb/{t}.{s}.vcf" for t in TOOLS for s in SETS)
    shell:
        """
        python3 {params.script} --sample {SAMPLE} --tools {params.tools} \\
            --source {input.source} --native {input.native} \\
            --tool-vcf {params.tool_vcfs} --rejected {params.rejected} --qc-rejected {params.qc} \\
            --lifted-annotated {input.annotated} --refchanged-ids {input.ids} \\
            --svdb {params.svdb} --links {output.links} --metrics {output.metrics}
        """
