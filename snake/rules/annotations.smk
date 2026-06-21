import os

# ─── Region annotation ────────────────────────────────────────────────────────
# Annotate source (hg19) and merged (hg38) VCFs with CUP and DISCREP flags.
#
# Chr-prefix notes:
#   hg19 source VCF  →  no chr prefix (e.g. "1")
#   hg19 BEDs        →  no chr prefix  → direct match, no renaming needed
#   hg38 merged VCF  →  chr prefix (e.g. "chr1")
#   hg38 BEDs        →  chr prefix     → direct match, no renaming needed
#
# hg38 BEDs are pre-processed to 3 columns (CHROM/FROM/TO) so bcftools
# annotate can use them directly without a runtime temp-file step.
#
# Region files live inside the container at /app/regions/ and are listed as
# params (not input) so Snakemake does not try to validate them on the host.

LIFTOVER_DIR = os.path.join(config["results_dir"], "liftover")
ANNOT_DIR    = os.path.join(config["results_dir"], "annotated")


rule annotate_regions_hg19:
    """Tag source hg19 VCF with CUPs, DISCREPs, and CENTEL flags."""
    input:
        vcf = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.vcf.gz",
        tbi = f"{LIFTOVER_DIR}/source/{SAMPLE}.hg19.vcf.gz.tbi",
    output:
        vcf = f"{ANNOT_DIR}/source/{SAMPLE}.hg19.annotated.vcf.gz",
        tbi = f"{ANNOT_DIR}/source/{SAMPLE}.hg19.annotated.vcf.gz.tbi"
    params:
        cups     = config["regions"]["cups_hg19"],
        cups_tbi = config["regions"]["cups_hg19"] + ".tbi",
        discreps = config["regions"]["discreps_hg19"],
        disc_tbi = config["regions"]["discreps_hg19"] + ".tbi",
        centel   = config["regions"]["centel_hg19"],
        cent_tbi = config["regions"]["centel_hg19"] + ".tbi",
        segdups  = config["regions"]["segdups_hg19"],
        segd_tbi = config["regions"]["segdups_hg19"] + ".tbi",
        header   = config["regions"]["header"]
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        TMP=$(mktemp -d)

        bcftools view --no-version -h {input.vcf} | grep -v '/ess/' > $TMP/clean_header.txt
        bcftools reheader -h $TMP/clean_header.txt {input.vcf} \
        | bcftools annotate \
            --no-version \
            -a {params.cups} \
            -c CHROM,FROM,TO \
            -m +CUPs \
            -h {params.header} \
        | bcftools annotate \
            --no-version \
            -a {params.discreps} \
            -c CHROM,FROM,TO \
            -m +DISCREPs \
        | bcftools annotate \
            --no-version \
            -a {params.centel} \
            -c CHROM,FROM,TO \
            -m +CENTEL \
        | bcftools annotate \
            --no-version \
            -a {params.segdups} \
            -c CHROM,FROM,TO \
            -m +SEGDUP \
            -Oz -o {output.vcf}

        tabix {output.vcf}
        rm -rf $TMP
        """


rule annotate_regions_hg38:
    """Tag merged hg38 VCF with CUPs and DISCREPs flags (Ormond 2021, Li 2021)."""
    input:
        vcf = f"{LIFTOVER_DIR}/merged/{SAMPLE}.merged.hg38.vcf.gz",
        tbi = f"{LIFTOVER_DIR}/merged/{SAMPLE}.merged.hg38.vcf.gz.tbi",
    output:
        vcf = f"{ANNOT_DIR}/merged/{SAMPLE}.merged.hg38.annotated.vcf.gz",
        tbi = f"{ANNOT_DIR}/merged/{SAMPLE}.merged.hg38.annotated.vcf.gz.tbi"
    params:
        cups    = config["regions"]["cups_hg38"],
        cups_tbi= config["regions"]["cups_hg38"] + ".tbi",
        discreps = config["regions"]["discreps_hg38"],
        disc_tbi= config["regions"]["discreps_hg38"] + ".tbi",
        header  = config["regions"]["header"]
    shell:
        """
        mkdir -p $(dirname {output.vcf})
        TMP=$(mktemp -d)

        bcftools view --no-version -h {input.vcf} | grep -v '/ess/' > $TMP/clean_header.txt
        bcftools reheader -h $TMP/clean_header.txt {input.vcf} \
        | bcftools annotate \
            --no-version \
            -a {params.cups} \
            -c CHROM,FROM,TO \
            -m +CUPs \
            -h {params.header} \
        | bcftools annotate \
            --no-version \
            -a {params.discreps} \
            -c CHROM,FROM,TO \
            -m +DISCREPs \
            -Oz -o {output.vcf}

        tabix {output.vcf}
        rm -rf $TMP
        """
