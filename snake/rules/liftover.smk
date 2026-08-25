import os

RESULTS_DIR  = config["results_dir"]
SAMPLE       = config["sample"]
LIFTOVER_DIR = os.path.join(RESULTS_DIR, "liftover")

TOOLS = config.get(
    "liftover_tools",
    ["crossmap", "bcftools", "picard", "transanno"],
)

_SV_TYPES  = config.get(
    "sv_types",
    ["DEL", "DUP", "DUP:TANDEM", "INS", "INV"],
)
_SV_FILTER = " || ".join([f'INFO/SVTYPE="{t}"' for t in _SV_TYPES])
_PASS_FLAG = "-f PASS,." if config.get("pass_only", True) else ""

_NORMALIZE_SV = config.get("tools", {}).get(
    "normalize_sv_script",
    "/app/liftover/normalize_sv_coordinates.py",
)


# ─── 1. Filter target SVs and stamp a unique source-derived ID ────────────────

rule prep_source_vcf:
    """
    Filter the source hg19 VCF and stamp every retained record with a
    source-coordinate ID plus a deterministic ordinal:

        CHROM_POS_SVTYPE_END_rN
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
        | awk 'BEGIN{{OFS="\\t"; n=0}} /^#/{{print; next}} {{
            split($8,a,";");
            svtype="";
            end_=$2;

            for(i in a){{
                if(a[i]~/^SVTYPE=/)
                    svtype=substr(a[i],8);

                if(a[i]~/^END=/)
                    end_=substr(a[i],5);
            }}

            n++;
            $3=$1 "_" $2 "_" svtype "_" end_ "_r" n;

            print;
        }}' \
        | bcftools sort -Oz -o {output.vcf}

        bcftools index -f -t {output.vcf}
        """


# ─── 2. Add chr prefix for UCSC-chain tools ──────────────────────────────────

rule normalize_contigs:
    input:
        vcf = rules.prep_source_vcf.output.vcf
    output:
        vcf     = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.chrnamed.vcf.gz",
        tbi     = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.chrnamed.vcf.gz.tbi",
        chr_map = f"{LIFTOVER_DIR}/source/chr_map.txt"
    shell:
        """
        bcftools view -h {input.vcf} \
        | awk 'BEGIN{{OFS="\\t"}} /^##contig=<ID=/ {{
            contig=$0;
            sub(/^##contig=<ID=/,"",contig);
            sub(/[,>].*$/,"",contig);
            if(contig !~ /^chr/){{
                target=(contig=="MT" ? "chrM" : "chr" contig);
                print contig,target;
            }}
        }}' > {output.chr_map}

        if [ -s "{output.chr_map}" ]; then
            bcftools annotate \
                --rename-chrs {output.chr_map} \
                {input.vcf} \
                -Oz -o {output.vcf}
        else
            bcftools view -Oz -o {output.vcf} {input.vcf}
        fi

        bcftools index -f -t {output.vcf}
        """


# ─── 3. CrossMap ──────────────────────────────────────────────────────────────

rule crossmap_liftover:
    input:
        vcf = rules.normalize_contigs.output.vcf,
        source = rules.prep_source_vcf.output.vcf
    output:
        vcf          = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.vcf.gz",
        tbi          = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.vcf.gz.tbi",
        unmap        = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.unmap.vcf",
        end_rejected = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.end_rejected.vcf"
    params:
        chain      = config["ref"]["liftover_chain"],
        fasta      = config["ref"]["hg38_fasta"],
        raw        = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.raw.vcf",
        normalized = f"{LIFTOVER_DIR}/crossmap/{SAMPLE}.hg38.normalized.vcf",
        normalizer = _NORMALIZE_SV
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/crossmap

        CrossMap vcf \
            --chromid l \
            {params.chain} \
            {input.vcf} \
            {params.fasta} \
            {params.raw}

        if [ -f {params.raw}.unmap ]; then
            mv {params.raw}.unmap {output.unmap}
        else
            touch {output.unmap}
        fi

        python3 {params.normalizer} \
            --source {input.source} \
            --lifted {params.raw} \
            --chain {params.chain} \
            --target-fasta {params.fasta} \
            --tool crossmap \
            --output {params.normalized} \
            --reject {output.end_rejected}

        bcftools sort \
            -Oz -o {output.vcf} \
            {params.normalized}

        bcftools index -f -t {output.vcf}

        rm -f \
            {params.raw} \
            {params.normalized}
        """


# ─── 4. bcftools +liftover ───────────────────────────────────────────────────

rule bcftools_liftover:
    input:
        vcf = rules.normalize_contigs.output.vcf,
        source = rules.prep_source_vcf.output.vcf
    output:
        vcf          = f"{LIFTOVER_DIR}/bcftools/{SAMPLE}.hg38.vcf.gz",
        tbi          = f"{LIFTOVER_DIR}/bcftools/{SAMPLE}.hg38.vcf.gz.tbi",
        rejected     = f"{LIFTOVER_DIR}/bcftools/{SAMPLE}.rejected.vcf",
        end_rejected = f"{LIFTOVER_DIR}/bcftools/{SAMPLE}.end_rejected.vcf"
    params:
        chain      = config["ref"]["liftover_chain"],
        srcfa      = config["ref"]["hg19_fasta_chr"],
        dstfa      = config["ref"]["hg38_fasta"],
        raw        = f"{LIFTOVER_DIR}/bcftools/{SAMPLE}.hg38.raw.vcf",
        normalized = f"{LIFTOVER_DIR}/bcftools/{SAMPLE}.hg38.normalized.vcf",
        normalizer = _NORMALIZE_SV
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/bcftools

        BCFTOOLS_PLUGINS=/usr/local/lib/bcftools_plugins \
        bcftools +liftover -Ov {input.vcf} \
            -- \
            -c {params.chain} \
            --src-fasta-ref {params.srcfa} \
            --fasta-ref {params.dstfa} \
            --no-left-align \
            --reject {output.rejected} \
            --write-reject \
            > {params.raw}

        python3 {params.normalizer} \
            --source {input.source} \
            --lifted {params.raw} \
            --chain {params.chain} \
            --target-fasta {params.dstfa} \
            --tool bcftools \
            --output {params.normalized} \
            --reject {output.end_rejected}

        bcftools sort \
            -Oz -o {output.vcf} \
            {params.normalized}

        bcftools index -f -t {output.vcf}

        rm -f \
            {params.raw} \
            {params.normalized}
        """


# ─── 5. Picard LiftoverVcf ───────────────────────────────────────────────────

rule picard_liftover:
    input:
        vcf = rules.normalize_contigs.output.vcf,
        source = rules.prep_source_vcf.output.vcf
    output:
        vcf          = f"{LIFTOVER_DIR}/picard/{SAMPLE}.hg38.vcf.gz",
        tbi          = f"{LIFTOVER_DIR}/picard/{SAMPLE}.hg38.vcf.gz.tbi",
        rejected     = f"{LIFTOVER_DIR}/picard/{SAMPLE}.rejected.vcf",
        end_rejected = f"{LIFTOVER_DIR}/picard/{SAMPLE}.end_rejected.vcf"
    params:
        chain      = config["ref"]["liftover_chain"],
        fasta      = config["ref"]["hg38_fasta"],
        raw        = f"{LIFTOVER_DIR}/picard/{SAMPLE}.hg38.raw.vcf",
        patched    = f"{LIFTOVER_DIR}/picard/{SAMPLE}.hg19.patched.vcf.gz",
        normalized = f"{LIFTOVER_DIR}/picard/{SAMPLE}.hg38.normalized.vcf",
        normalizer = _NORMALIZE_SV
    resources:
        mem_mb = 7000
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/picard

        bcftools view {input.vcf} \
            | sed 's/^##fileformat=VCFv4\\.[4-9]/##fileformat=VCFv4.3/' \
            | bgzip -c > {params.patched}

        bcftools index -f -t {params.patched}

        JAVA_TOOL_OPTIONS="-Xmx6g" \
        picard LiftoverVcf \
            -INPUT {params.patched} \
            -OUTPUT {params.raw} \
            -CHAIN {params.chain} \
            -REJECT {output.rejected} \
            -REFERENCE_SEQUENCE {params.fasta} \
            -LIFTOVER_MIN_MATCH 0.0 \
            -WARN_ON_MISSING_CONTIG true \
            -VALIDATION_STRINGENCY LENIENT

        python3 {params.normalizer} \
            --source {input.source} \
            --lifted {params.raw} \
            --chain {params.chain} \
            --target-fasta {params.fasta} \
            --tool picard \
            --output {params.normalized} \
            --reject {output.end_rejected}

        bcftools sort \
            -Oz -o {output.vcf} \
            {params.normalized}

        bcftools index -f -t {output.vcf}

        rm -f \
            {params.raw} \
            {params.patched} \
            {params.patched}.tbi \
            {params.normalized}
        """


# ─── 6. Transanno ─────────────────────────────────────────────────────────────

rule transanno_liftover:
    input:
        vcf = rules.prep_source_vcf.output.vcf
    output:
        vcf          = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.hg38.vcf.gz",
        tbi          = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.hg38.vcf.gz.tbi",
        rejected     = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.rejected.vcf",
        end_rejected = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.end_rejected.vcf"
    log:
        f"{LIFTOVER_DIR}/transanno/{SAMPLE}.log"
    params:
        chain      = config["ref"]["transanno_chain"],
        srcfa      = config["ref"]["hg19_fasta"],
        dstfa      = config["ref"]["hg38_fasta"],
        raw        = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.hg38.raw.vcf",
        normalized = f"{LIFTOVER_DIR}/transanno/{SAMPLE}.hg38.normalized.vcf",
        normalizer = _NORMALIZE_SV
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/transanno

        exec > {log} 2>&1
        set -x

        INPUT={LIFTOVER_DIR}/transanno/{SAMPLE}.hg19.nochr.vcf.gz
        CHR_MAP=$(mktemp)

        bcftools view -h {input.vcf} \
        | awk 'BEGIN{{OFS="\\t"}} /^##contig=<ID=/ {{
            contig=$0;
            sub(/^##contig=<ID=/,"",contig);
            sub(/[,>].*$/,"",contig);
            if(contig ~ /^chr/){{
                target=substr(contig,4);
                if(target=="M") target="MT";
                print contig,target;
            }}
        }}' > $CHR_MAP

        if [ -s "$CHR_MAP" ]; then
            bcftools annotate \
                --rename-chrs $CHR_MAP \
                {input.vcf} \
                -Oz -o $INPUT
        else
            bcftools view -Oz -o $INPUT {input.vcf}
        fi

        bcftools index -f -t $INPUT
        rm -f $CHR_MAP

        transanno liftvcf \
            --original-assembly {params.srcfa} \
            --new-assembly {params.dstfa} \
            --chain {params.chain} \
            --vcf $INPUT \
            --output {params.raw} \
            --fail {output.rejected}

        python3 {params.normalizer} \
            --source {input.vcf} \
            --lifted {params.raw} \
            --chain {params.chain} \
            --target-fasta {params.dstfa} \
            --tool transanno \
            --output {params.normalized} \
            --reject {output.end_rejected}

        bcftools sort \
            -Oz -o {output.vcf} \
            {params.normalized}

        bcftools index -f -t {output.vcf}

        rm -f \
            {params.raw} \
            {params.normalized} \
            $INPUT \
            $INPUT.tbi
        """


# ─── 7. Add LIFTOVER_TOOL ─────────────────────────────────────────────────────

def _tool_vcf(wildcards):
    return f"{LIFTOVER_DIR}/{wildcards.tool}/{SAMPLE}.hg38.vcf.gz"


rule tag_liftover_tool:
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

        bcftools query \
            -f '%CHROM\\t%POS\\t%ID\\t{wildcards.tool}\\n' \
            {input.vcf} \
        | bgzip -c > $ANN

        tabix -f -s1 -b2 -e2 $ANN

        HDR=$(mktemp)
        echo '##INFO=<ID=LIFTOVER_TOOL,Number=.,Type=String,Description="Liftover tool(s) that produced this exact target representation">' \
            > $HDR

        bcftools annotate \
            -a $ANN \
            -c CHROM,POS,~ID,INFO/LIFTOVER_TOOL \
            -h $HDR \
            -Oz -o {output.vcf} \
            {input.vcf}

        bcftools index -f -t {output.vcf}

        rm -f \
            $HDR \
            $ANN \
            $ANN.tbi
        """


# ─── 8. Merge matching target intervals ──────────────────────────────────────

rule merge_liftovers:
    """
    merge_tools.py buffers each CHROM/POS bucket and only collapses records
    agreeing on CHROM/POS/END/source-derived ID.
    """
    input:
        vcfs = expand(
            f"{LIFTOVER_DIR}/tagged/{{tool}}.hg38.tagged.vcf.gz",
            tool=TOOLS,
        )
    output:
        vcf = f"{LIFTOVER_DIR}/merged/{SAMPLE}.merged.hg38.vcf.gz",
        tbi = f"{LIFTOVER_DIR}/merged/{SAMPLE}.merged.hg38.vcf.gz.tbi"
    params:
        script = config.get("tools", {}).get(
            "merge_script",
            "/app/liftover/merge_tools.py",
        ),
        picard_rejected = (
            f"{LIFTOVER_DIR}/picard/{SAMPLE}.rejected.vcf"
            if "picard" in TOOLS
            else ""
        )
    shell:
        """
        mkdir -p {LIFTOVER_DIR}/merged

        PICARD_FLAG=""
        if [ -n "{params.picard_rejected}" ]; then
            PICARD_FLAG="--picard-rejected {params.picard_rejected}"
        fi

        bcftools concat \
            --allow-overlaps \
            -Ou \
            {input.vcfs} \
        | bcftools sort -Ov \
        | python3 {params.script} $PICARD_FLAG \
        | bgzip -c > {output.vcf}

        bcftools index -f -t {output.vcf}
        """
