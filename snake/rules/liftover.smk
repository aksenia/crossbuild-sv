import os

RESULTS_DIR  = config['results_dir']
SAMPLE       = config['sample']
LIFTOVER_DIR = os.path.join(RESULTS_DIR, 'liftover')
TOOLS        = config.get('liftover_tools', ['crossmap', 'bcftools', 'picard', 'transanno'])

# Build SVTYPE include-filter expression for bcftools
_SV_TYPES    = config.get('sv_types', ['DEL', 'DUP', 'INS', 'INV'])
_SV_FILTER   = ' || '.join([f'INFO/SVTYPE="{t}"' for t in _SV_TYPES])
_PASS_FLAG   = '-f PASS,.' if config.get('pass_only', True) else ''
_CONTIGS     = config.get('contigs') or []
_TARGETS     = f"-t {','.join(map(str, _CONTIGS))}" if _CONTIGS else ''

# awk: print "<old>\t<new>" for every ##contig ID, from header lines on stdin.
# Uses the full header (not only contigs seen in records) so renamed headers
# never mix numeric and chr names.
# (params are inserted verbatim, so single braces here.)
_ADD_CHR = r'''awk -F'[=,>]' '/^##contig=<ID=/ && $3 !~ /^chr/ {print $3 "\t" ($3=="MT" ? "chrM" : "chr" $3)}' '''
_DEL_CHR = r'''awk -F'[=,>]' '/^##contig=<ID=chr/ {c=substr($3,4); print $3 "\t" (c=="M" ? "MT" : c)}' '''

# ─── 1. Filter to target SV types, optionally require PASS, stamp hg19 IDs ───
rule prep_source_vcf:
    """
    Keep only target SVTYPE values from config sv_types (optionally PASS only,
    optionally a contig subset). Stamp each record with a stable ID
    CHROM_POS_SVTYPE_END (hg19 coords). Distinct records sharing that key get
    _2, _3, ... in file order, so IDs stay unique and do not depend on which
    contigs were selected.
    """
    input:
        vcf = config["input_vcf"]
    output:
        vcf = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.vcf.gz",
        tbi = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.vcf.gz.tbi"
    params:
        pass_flag = _PASS_FLAG,
        sv_filter = _SV_FILTER,
        targets   = _TARGETS
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/source
        bcftools view \
            {params.pass_flag} {params.targets} \
            -i '{params.sv_filter}' \
            {input.vcf} \
        | awk 'BEGIN{{OFS="\\t"}} /^#/{{print; next}} {{
            split($8,a,";"); svtype=""; end_=$2;
            for(i in a){{
                if(a[i]~/^SVTYPE=/) svtype=substr(a[i],8);
                if(a[i]~/^END=/) end_=substr(a[i],5)
            }}
            id=$1"_"$2"_"svtype"_"end_;
            $3=(++seen[id]>1) ? id"_"seen[id] : id; print
        }}' \
        | bcftools sort -Oz -o {output.vcf}
        tabix -f {output.vcf}
        """

# ─── 2. Normalize contig names to chr-prefix for UCSC chain ──────────────────
rule normalize_contigs:
    """
    Ensure contig names use chr-prefix (e.g. 22 -> chr22, MT -> chrM) to match
    the UCSC chain file. If the VCF already uses chr-prefix this is a
    pass-through. All tools except transanno use this output.
    """
    input:
        vcf = rules.prep_source_vcf.output.vcf
    output:
        vcf     = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.chrnamed.vcf.gz",
        tbi     = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.chrnamed.vcf.gz.tbi",
        chr_map = f"{LIFTOVER_DIR}/source/chr_map.txt"
    params:
        add_chr = _ADD_CHR
    shell:
        """
        bcftools view -h {input.vcf} | {params.add_chr} > {output.chr_map}
        if [ -s "{output.chr_map}" ]; then
            bcftools annotate --rename-chrs {output.chr_map} {input.vcf} -Oz -o {output.vcf}
        else
            bcftools view -Oz -o {output.vcf} {input.vcf}
        fi
        tabix -f {output.vcf}
        """

# ─── 3. CrossMap liftover ─────────────────────────────────────────────────────
rule crossmap_liftover:
    """Lift to hg38 with CrossMap vcf; --chromid l keeps chr-prefixed names."""
    input:
        vcf = rules.normalize_contigs.output.vcf
    output:
        vcf   = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.raw.vcf.gz",
        unmap = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.unmap.vcf"
    params:
        chain = config["ref"]["liftover_chain"],
        fasta = config["ref"]["hg38_fasta"],
        raw   = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.tmp.vcf"
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/crossmap
        CrossMap vcf --chromid l {params.chain} {input.vcf} {params.fasta} {params.raw}
        [ -f {params.raw}.unmap ] \
            && mv {params.raw}.unmap {output.unmap} \
            || touch {output.unmap}
        bcftools sort -Oz -o {output.vcf} {params.raw}
        rm -f {params.raw}
        """

# ─── 4. bcftools +liftover ────────────────────────────────────────────────────
rule bcftools_liftover:
    """
    Lift to hg38 with bcftools +liftover.
    --no-left-align: keep the tool's own placement; left-alignment is applied
    to every tool alike in normalize_liftover, which also recomputes END.
    """
    input:
        vcf = rules.normalize_contigs.output.vcf
    output:
        vcf      = f"{LIFTOVER_DIR}/bcftools/{SAMPLE}.hg38.raw.vcf.gz",
        rejected = f"{LIFTOVER_DIR}/bcftools/{SAMPLE}.rejected.vcf"
    params:
        chain = config["ref"]["liftover_chain"],
        srcfa = config["ref"]["hg19_fasta_chr"],
        dstfa = config["ref"]["hg38_fasta"]
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/bcftools
        bcftools +liftover -Ov {input.vcf} \
            -- -c {params.chain} \
               --src-fasta-ref {params.srcfa} \
               --fasta-ref {params.dstfa} \
               --no-left-align \
               --reject {output.rejected} \
               --write-reject \
            | bcftools sort -Oz -o {output.vcf}
        """

# ─── 5. Picard LiftoverVcf ────────────────────────────────────────────────────
rule picard_liftover:
    """
    Lift to hg38 with Picard LiftoverVcf.
    LIFTOVER_MIN_MATCH=0.0: prevents rejection of large SVs whose sequence
    spans a chain break with low fraction of bases mapping to the same interval.
    Header is patched to VCFv4.3 as Picard rejects v4.4+. No ALT pre-filter:
    an IUPAC regex here used to drop every symbolic allele (<DUP:TANDEM> ...).
    """
    input:
        vcf = rules.normalize_contigs.output.vcf
    output:
        vcf      = f"{LIFTOVER_DIR}/picard/{SAMPLE}.hg38.raw.vcf.gz",
        rejected = f"{LIFTOVER_DIR}/picard/{SAMPLE}.rejected.vcf"
    params:
        chain     = config["ref"]["liftover_chain"],
        fasta     = config["ref"]["hg38_fasta"],
        java_opts = config.get("picard_java_opts", "-Xmx6g")
    resources:
        mem_mb = 7000
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/picard
        RAW={LIFTOVER_DIR}/picard/{SAMPLE}.hg38.tmp.vcf
        PATCHED={LIFTOVER_DIR}/picard/{SAMPLE}.hg19.patched.vcf.gz

        bcftools view {input.vcf} \
            | sed 's/^##fileformat=VCFv4\\.[4-9]/##fileformat=VCFv4.3/' \
            | bgzip > $PATCHED
        tabix -f $PATCHED

        JAVA_OPTS="{params.java_opts}" picard LiftoverVcf \
            -INPUT $PATCHED \
            -OUTPUT $RAW \
            -CHAIN {params.chain} \
            -REJECT {output.rejected} \
            -REFERENCE_SEQUENCE {params.fasta} \
            -LIFTOVER_MIN_MATCH 0.0 \
            -WARN_ON_MISSING_CONTIG true \
            -VALIDATION_STRINGENCY LENIENT \
            -CREATE_INDEX false
        bcftools sort -Oz -o {output.vcf} $RAW
        rm -f $RAW $PATCHED $PATCHED.tbi
        """

# ─── 6. Transanno liftover ────────────────────────────────────────────────────
rule transanno_liftover:
    """
    Lift to hg38 with Transanno liftvcf.
    Takes the non-chr-prefixed source VCF directly (transanno chain uses numeric
    contig names). Strips chr prefix if present (chrM -> MT).
    Note: transanno may reject large sequence-resolved SVs.
    """
    input:
        vcf = rules.prep_source_vcf.output.vcf
    output:
        vcf      = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.hg38.raw.vcf.gz",
        rejected = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.rejected.vcf"
    log:
        f"{LIFTOVER_DIR}/transanno/{SAMPLE}.log"
    params:
        chain   = config["ref"]["transanno_chain"],
        srcfa   = config["ref"]["hg19_fasta"],
        dstfa   = config["ref"]["hg38_fasta"],
        raw     = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.hg38.tmp.vcf",
        del_chr = _DEL_CHR
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/transanno
        exec > {log} 2>&1
        set -x
        INPUT={LIFTOVER_DIR}/transanno/{SAMPLE}.hg19.nochr.vcf.gz

        CHR_MAP=$(mktemp)
        bcftools view -h {input.vcf} | {params.del_chr} > $CHR_MAP
        if [ -s "$CHR_MAP" ]; then
            bcftools annotate --rename-chrs $CHR_MAP {input.vcf} -Oz -o $INPUT
        else
            bcftools view -Oz -o $INPUT {input.vcf}
        fi
        tabix -f $INPUT
        rm -f $CHR_MAP

        transanno liftvcf \
            --original-assembly {params.srcfa} \
            --new-assembly      {params.dstfa} \
            --chain             {params.chain} \
            --vcf               $INPUT \
            --output            {params.raw} \
            --fail              {output.rejected}
        bcftools sort -Oz -o {output.vcf} {params.raw}
        rm -f {params.raw} $INPUT $INPUT.tbi
        """

# ─── 7. Normalize each tool's output the same way ───────────────────────────
rule normalize_liftover:
    """
    Put every tool's hg38 records into one canonical representation so the
    merge compares coordinates, not allele spelling: bcftools pads REF/ALT,
    transanno drops a shared padding base (POS+1). `bcftools norm`
    left-aligns and trims; finalize_sv.py then checks REF against hg38, sets
    END from the normalized record (chain-mapped for symbolic intervals),
    and tags LIFTOVER_TOOL. Records failing those checks go to qc_rejected.
    """
    input:
        raw    = f"{LIFTOVER_DIR}/{{tool}}/{SAMPLE}.hg38.raw.vcf.gz",
        source = rules.prep_source_vcf.output.vcf
    output:
        vcf      = f"{LIFTOVER_DIR}/{{tool}}/{SAMPLE}.hg38.vcf.gz",
        tbi      = f"{LIFTOVER_DIR}/{{tool}}/{SAMPLE}.hg38.vcf.gz.tbi",
        rejected = f"{LIFTOVER_DIR}/{{tool}}/{SAMPLE}.qc_rejected.vcf"
    log:
        f"{LIFTOVER_DIR}/{{tool}}/{SAMPLE}.normalize.log"
    wildcard_constraints:
        tool = "|".join(TOOLS)
    params:
        fasta  = config["ref"]["hg38_fasta"],
        chain  = lambda w: config["ref"]["transanno_chain" if w.tool == "transanno" else "liftover_chain"],
        tmp    = f"{LIFTOVER_DIR}/{{tool}}/{SAMPLE}.hg38.norm.tmp.vcf",
        script = f"{SCRIPTS}/finalize_sv.py"
    shell:
        """
        # -c w: warn on REF mismatch; finalize_sv.py rejects those records.
        bcftools norm -f {params.fasta} -c w {input.raw} 2> {log} \
        | python3 {params.script} --tool {wildcards.tool} --source {input.source} \
            --chain {params.chain} --fasta {params.fasta} \
            --output {params.tmp} --reject {output.rejected} 2>> {log}
        bcftools sort -Oz -o {output.vcf} {params.tmp}
        tabix -f {output.vcf}
        rm -f {params.tmp}
        """

# ─── 8. Merge all normalized VCFs into one hg38 VCF ─────────────────────────
rule merge_liftovers:
    """
    Concatenate per-tool hg38 VCFs, sort, then collapse records with identical
    CHROM/POS/END/ID by joining their LIFTOVER_TOOL values. After
    normalization a remaining POS or END difference is a real coordinate
    disagreement and stays on its own line with its own LIFTOVER_TOOL tag.
    """
    input:
        vcfs = expand(f"{LIFTOVER_DIR}/{{tool}}/{SAMPLE}.hg38.vcf.gz", tool=TOOLS),
        tbis = expand(f"{LIFTOVER_DIR}/{{tool}}/{SAMPLE}.hg38.vcf.gz.tbi", tool=TOOLS)
    output:
        vcf = f"{LIFTOVER_DIR}/merged/{SAMPLE}.merged.hg38.vcf.gz",
        tbi = f"{LIFTOVER_DIR}/merged/{SAMPLE}.merged.hg38.vcf.gz.tbi"
    params:
        script          = f"{SCRIPTS}/merge_tools.py",
        picard_rejected = (f"{LIFTOVER_DIR}/picard/{SAMPLE}.rejected.vcf"
                           if 'picard' in TOOLS else "")
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/merged
        PICARD_FLAG=""
        [ -n "{params.picard_rejected}" ] && PICARD_FLAG="--picard-rejected {params.picard_rejected}"
        bcftools concat --allow-overlaps {input.vcfs} \
            | bcftools sort \
            | python3 {params.script} $PICARD_FLAG \
            | bgzip > {output.vcf}
        tabix -f {output.vcf}
        """
