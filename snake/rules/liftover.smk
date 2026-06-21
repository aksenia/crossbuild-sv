import os

RESULTS_DIR  = config['results_dir']
SAMPLE       = config['sample']
LIFTOVER_DIR = os.path.join(RESULTS_DIR, 'liftover')
TOOLS        = config.get('liftover_tools', ['crossmap', 'bcftools', 'picard', 'transanno'])

# Build SVTYPE include-filter expression for bcftools
_SV_TYPES    = config.get('sv_types', ['DEL', 'DUP', 'DUP:TANDEM', 'INS', 'INV'])
_SV_FILTER   = ' | '.join([f'INFO/SVTYPE="{t}"' for t in _SV_TYPES])
_PASS_FLAG   = '-f PASS,.' if config.get('pass_only', True) else ''

# ─── 1. Filter to target SV types, optionally require PASS, stamp hg19 IDs ───
rule prep_source_vcf:
    """
    Keep only target SVTYPE values (DEL/DUP/INS/INV etc.) from config sv_types.
    Optionally restrict to PASS records (pass_only: true).
    Stamp each record with a stable ID: CHROM_POS_SVTYPE_END in hg19 coords.
    """
    input:
        vcf = config["input_vcf"]
    output:
        vcf = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.vcf.gz",
        tbi = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.vcf.gz.tbi"
    params:
        pass_flag = _PASS_FLAG,
        sv_filter = _SV_FILTER
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/source
        bcftools view \
            {params.pass_flag} \
            -i '{params.sv_filter}' \
            {input.vcf} \
        | awk 'BEGIN{{OFS="\\t"}} /^#/{{print; next}} {{
            split($8,a,";"); svtype=""; end_=$2;
            for(i in a){{
                if(a[i]~/^SVTYPE=/) svtype=substr(a[i],8);
                if(a[i]~/^END=/) end_=substr(a[i],5)
            }}
            $3=$1"_"$2"_"svtype"_"end_; print
        }}' \
        | bcftools sort -Oz -o {output.vcf}
        tabix {output.vcf}
        """

# ─── 2. Normalize contig names to chr-prefix for UCSC chain ──────────────────
rule normalize_contigs:
    """
    Ensure contig names use chr-prefix (e.g. 22 -> chr22) to match the UCSC
    chain file. If the VCF already uses chr-prefix this is a pass-through.
    All tools except transanno use this output.
    """
    input:
        vcf = rules.prep_source_vcf.output.vcf
    output:
        vcf     = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.chrnamed.vcf.gz",
        tbi     = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.chrnamed.vcf.gz.tbi",
        chr_map = f"{LIFTOVER_DIR}/source/chr_map.txt"
    shell:
        """
        bcftools view -h {input.vcf} \
            | grep '^##contig' | grep -oP 'ID=\K[^,>]+' \
            | grep -v '^chr' \
            | awk '{{print $0"\tchr"$0}}' > {output.chr_map} || true

        if [ -s "{output.chr_map}" ]; then
            bcftools annotate --rename-chrs {output.chr_map} {input.vcf} -Oz -o {output.vcf}
        else
            bcftools view -Oz -o {output.vcf} {input.vcf}
        fi
        tabix {output.vcf}
        """

# ─── 3. CrossMap liftover ─────────────────────────────────────────────────────
rule crossmap_liftover:
    """Lift to hg38 with CrossMap vcf."""
    input:
        vcf = rules.normalize_contigs.output.vcf
    output:
        vcf   = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.vcf.gz",
        tbi   = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.vcf.gz.tbi",
        unmap = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.unmap.vcf"
    params:
        chain = config["ref"]["liftover_chain"],
        fasta = config["ref"]["hg38_fasta"],
        raw   = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.raw.vcf"
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/crossmap
        CrossMap vcf {params.chain} {input.vcf} {params.fasta} {params.raw}
        [ -f {params.raw}.unmap ] \
            && mv {params.raw}.unmap {output.unmap} \
            || touch {output.unmap}
        bcftools sort -Oz -o {output.vcf} {params.raw}
        tabix {output.vcf}
        rm -f {params.raw}
        """

# ─── 4. bcftools +liftover ────────────────────────────────────────────────────
rule bcftools_liftover:
    """
    Lift to hg38 with bcftools +liftover.
    --no-left-align: prevents re-anchoring large indel/SV sequences after lift,
    which would shift POS and corrupt the INFO/END coordinate.
    """
    input:
        vcf = rules.normalize_contigs.output.vcf
    output:
        vcf      = f"{LIFTOVER_DIR}/bcftools/{SAMPLE}.hg38.vcf.gz",
        tbi      = f"{LIFTOVER_DIR}/bcftools/{SAMPLE}.hg38.vcf.gz.tbi",
        rejected = f"{LIFTOVER_DIR}/bcftools/{SAMPLE}.rejected.vcf"
    params:
        chain = config["ref"]["liftover_chain"],
        srcfa = config["ref"]["hg19_fasta_chr"],
        dstfa = config["ref"]["hg38_fasta"]
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/bcftools
        BCFTOOLS_PLUGINS=/usr/local/lib/bcftools_plugins \
        bcftools +liftover -Ov {input.vcf} \
            -- -c {params.chain} \
               --src-fasta-ref {params.srcfa} \
               --fasta-ref {params.dstfa} \
               --no-left-align \
               --reject {output.rejected} \
               --write-reject \
            | bcftools sort -Oz -o {output.vcf}
        tabix {output.vcf}
        """

# ─── 5. Picard LiftoverVcf ────────────────────────────────────────────────────
rule picard_liftover:
    """
    Lift to hg38 with Picard LiftoverVcf.
    LIFTOVER_MIN_MATCH=0.0: prevents rejection of large SVs whose sequence
    spans a chain break with low fraction of bases mapping to the same interval.
    Header is patched to VCFv4.3 as Picard rejects v4.4+.
    """
    input:
        vcf = rules.normalize_contigs.output.vcf
    output:
        vcf      = f"{LIFTOVER_DIR}/picard/{SAMPLE}.hg38.vcf.gz",
        tbi      = f"{LIFTOVER_DIR}/picard/{SAMPLE}.hg38.vcf.gz.tbi",
        rejected = f"{LIFTOVER_DIR}/picard/{SAMPLE}.rejected.vcf"
    params:
        chain = config["ref"]["liftover_chain"],
        fasta = config["ref"]["hg38_fasta"]
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/picard
        RAW={LIFTOVER_DIR}/picard/{SAMPLE}.hg38.raw.vcf
        PATCHED={LIFTOVER_DIR}/picard/{SAMPLE}.hg19.patched.vcf.gz

        bcftools view -e 'ALT~"[RYSWKMBDHV]"' {input.vcf} \
            | sed 's/^##fileformat=VCFv4\.[4-9]/##fileformat=VCFv4.3/' \
            | bgzip > $PATCHED
        tabix $PATCHED

        picard LiftoverVcf \
            -INPUT $PATCHED \
            -OUTPUT $RAW \
            -CHAIN {params.chain} \
            -REJECT {output.rejected} \
            -REFERENCE_SEQUENCE {params.fasta} \
            -LIFTOVER_MIN_MATCH 0.0 \
            -WARN_ON_MISSING_CONTIG true \
            -VALIDATION_STRINGENCY LENIENT
        bcftools sort -Oz -o {output.vcf} $RAW
        tabix {output.vcf}
        rm -f $RAW $PATCHED $PATCHED.tbi
        """

# ─── 6. Transanno liftover ────────────────────────────────────────────────────
rule transanno_liftover:
    """
    Lift to hg38 with Transanno liftvcf.
    Takes the non-chr-prefixed source VCF directly (transanno chain uses numeric
    contig names). Strips chr prefix if present.
    Note: transanno may reject large sequence-resolved SVs; rejected records
    simply won't carry the transanno LIFTOVER_TOOL tag.
    """
    input:
        vcf = rules.prep_source_vcf.output.vcf
    output:
        vcf      = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.hg38.vcf.gz",
        tbi      = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.hg38.vcf.gz.tbi",
        rejected = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.rejected.vcf"
    log:
        f"{LIFTOVER_DIR}/transanno/{SAMPLE}.log"
    params:
        chain = config["ref"]["transanno_chain"],
        srcfa = config["ref"]["hg19_fasta"],
        dstfa = config["ref"]["hg38_fasta"],
        raw   = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.hg38.raw.vcf"
    shell:
        """
        exec > {log} 2>&1
        set -x
        mkdir -p {LIFTOVER_DIR}/transanno
        INPUT={LIFTOVER_DIR}/transanno/{SAMPLE}.hg19.nochr.vcf.gz

        CHR_MAP=$(mktemp)
        bcftools view -h {input.vcf} \
            | grep '^##contig' | grep -oP 'ID=\K[^,>]+' \
            | grep '^chr' \
            | awk '{{print $0"\t"substr($0,4)}}' > $CHR_MAP || true
        if [ -s "$CHR_MAP" ]; then
            bcftools annotate --rename-chrs $CHR_MAP {input.vcf} -Oz -o $INPUT
        else
            bcftools view -Oz -o $INPUT {input.vcf}
        fi
        tabix $INPUT
        rm -f $CHR_MAP

        transanno liftvcf \
            --original-assembly {params.srcfa} \
            --new-assembly      {params.dstfa} \
            --chain             {params.chain} \
            --vcf               $INPUT \
            --output            {params.raw} \
            --fail              {output.rejected}
        bcftools sort -Oz -o {output.vcf} {params.raw}
        tabix {output.vcf}
        rm -f {params.raw} $INPUT $INPUT.tbi
        """

# ─── 7. Tag each VCF with LIFTOVER_TOOL INFO field ───────────────────────────
def _tool_vcf(wildcards):
    return f"{LIFTOVER_DIR}/{wildcards.tool}/{SAMPLE}.hg38.vcf.gz"

rule tag_liftover_tool:
    """
    Annotate each per-tool hg38 VCF with LIFTOVER_TOOL=<tool> in INFO.
    Matches on CHROM/POS/ID (ID = CHROM_POS_SVTYPE_END from hg19) rather than
    CHROM/POS/REF/ALT, since SV REF/ALT sequences can be thousands of bases.
    """
    input:
        vcf = _tool_vcf
    output:
        vcf = f"{LIFTOVER_DIR}/tagged/{{tool}}.hg38.tagged.vcf.gz",
        tbi = f"{LIFTOVER_DIR}/tagged/{{tool}}.hg38.tagged.vcf.gz.tbi"
    wildcard_constraints:
        tool = "|".join(TOOLS)
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/tagged
        ANN={LIFTOVER_DIR}/tagged/{wildcards.tool}.ann.tsv.gz

        bcftools query -f '%CHROM\t%POS\t%ID\t{wildcards.tool}\n' \
            {input.vcf} \
        | sort -k1,1V -k2,2n \
        | bgzip > $ANN
        tabix -s1 -b2 -e2 $ANN

        HDR=$(mktemp)
        echo '##INFO=<ID=LIFTOVER_TOOL,Number=1,Type=String,Description="Liftover tool that produced this record">' > $HDR
        bcftools annotate \
            -a $ANN \
            -c CHROM,POS,ID,INFO/LIFTOVER_TOOL \
            -h $HDR \
            -Oz -o {output.vcf} \
            {input.vcf}
        rm -f $HDR
        tabix {output.vcf}
        """

# ─── 8. Merge all tagged VCFs into one hg38 VCF ──────────────────────────────
rule merge_liftovers:
    """
    Concatenate per-tool hg38 VCFs, sort, then collapse records with identical
    CHROM/POS/ID by joining their LIFTOVER_TOOL values.
    Records where tools disagree on position (same hg19 ID, different hg38 POS)
    are kept as separate lines — each carrying its own LIFTOVER_TOOL tag.
    """
    input:
        vcfs = expand(
            f"{LIFTOVER_DIR}/tagged/{{tool}}.hg38.tagged.vcf.gz",
            tool=TOOLS
        )
    output:
        vcf = f"{LIFTOVER_DIR}/merged/{SAMPLE}.merged.hg38.vcf.gz",
        tbi = f"{LIFTOVER_DIR}/merged/{SAMPLE}.merged.hg38.vcf.gz.tbi"
    params:
        script          = config.get("tools", {}).get("merge_script", "/app/liftover/merge_tools.py"),
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
        tabix {output.vcf}
        """
