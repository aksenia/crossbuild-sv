


"""
rules:
an id is not unique to one row in the id tsv file, one if can appear on mutliple lines because one row per tool,
so cannot say number of lines = nr of variants

an id in refchanged does not automatically mean a swap! the pipeline flags any changes in swap, refchanged or multialleic alt
therefore need to distinguish between the allele states

we join using the original source IDs. coordinates and allele description can change during liftover

we keep the original liftover - need a swap and original homozygous-alternative genotype to be able to distinguish between the allele states

conflicts need to stay visible, so if tools diagree we want to keep them

have a shared classification oer source ID


maybe SVTYPE correction



"""


# Inspect evidence for REF-changed source variants before classification.

RUN=$(cd "$1" && pwd)
SAMPLE=HG002

# input files
REFCHANGED_IDS="$RUN/liftover/merged/$SAMPLE.refchanged_ids.tsv"
SOURCE_VCF="$RUN/liftover/source/$SAMPLE.hg19.vcf.gz"
LINKS="$RUN/comparison/$SAMPLE.links.tsv.gz"
METRICS="$RUN/comparison/$SAMPLE.metrics.tsv"

# Require each shared input to be readable.
for file in "$REFCHANGED_IDS" "$SOURCE_VCF" "$LINKS" "$METRICS"; do
    if [[ ! -r "$file" ]]; then
        printf 'Error: missing or unreadable file: %s\n' "$file" >&2
        exit 1
    fi
done



for tool in "${TOOLS[@]}"; do
    for suffix in hg38.raw.vcf.gz hg38.vcf.gz; do
        file="$RUN/liftover/$tool/$SAMPLE.$suffix"

        if [[ ! -r "$file" ]]; then
            printf 'Error: missing or unreadable file: %s\n' "$file" >&2
            exit 1
        fi
    done
done


# Count unique source IDs
N_REFCHANGED=$(
    cut -f1 "$REFCHANGED_IDS" \
        | LC_ALL=C sort -u \
        | wc -l \
        | tr -d '[:space:]'
)

printf '\nUnique REF-changed source IDs: %s\n' "$N_REFCHANGED"
printf 'The step 2 plan records 171 IDs for the wg4 run.\n'

printf '\nFirst 10 rows of the REF-changed ID/tool table:\n'
head -n 10 "$REFCHANGED_IDS"



# Inspect allele-state tag definitions in each tool's raw output.
for tool in "${TOOLS[@]}"; do
    vcf="$RUN/liftover/$tool/$SAMPLE.hg38.raw.vcf.gz"

    printf '\n=== %s: allele-state tag definitions ===\n' "$tool"

    bcftools view -h "$vcf" \
        | awk '
            /^##INFO=<ID=(SWAP|REF_CHANGED),/ {
                print
                found = 1
            }
            END {
                if (!found)
                    print "No SWAP or REF_CHANGED definition in this header."
            }
        '

    printf '\nSample names in this VCF:\n'
    bcftools query -l "$vcf"
done

# Check the source VCF sample names too.
printf '\n=== Original hg19 source: sample names ===\n'
bcftools query -l "$SOURCE_VCF"


# Store inspection results alongside the run data.
OUT="$RUN/step2_refstate"
mkdir -p "$OUT"

# Save one entry per source ID.
cut -f1 "$REFCHANGED_IDS" \
    | LC_ALL=C sort -u \
    > "$OUT/refchanged.ids.txt"

# Extract records belonging to the REF-changed source IDs.
# Print one row per VCF record and sample.
extract_evidence() {
    local label="$1"
    local vcf="$2"

    bcftools query \
        -f '[%ID\t%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO\t%SAMPLE\t%GT\n]' \
        "$vcf" \
        | awk -F '\t' -v OFS='\t' -v label="$label" '
            FILENAME == ARGV[1] {
                ids[$1] = 1
                next
            }
            $1 in ids {
                print label, $0
            }
        ' "$OUT/refchanged.ids.txt" -
}

# Collect original alleles/genotypes and each tool's raw evidence.
{
    printf 'origin\tsource_id\tchrom\tpos\tref\talt\tfilter\tinfo\tsample\tgt\n'

    extract_evidence "source_hg19" "$SOURCE_VCF"

    for tool in "${TOOLS[@]}"; do
        extract_evidence \
            "$tool" \
            "$RUN/liftover/$tool/$SAMPLE.hg38.raw.vcf.gz"
    done
} > "$OUT/refchanged.evidence.tsv"

printf '\nEvidence table saved to:\n%s\n' \
    "$OUT/refchanged.evidence.tsv"




    # Collect the IDs successfully extracted from the original source VCF.
awk -F '\t' '
    NR > 1 && $1 == "source_hg19" {
        print $2
    }
' "$OUT/refchanged.evidence.tsv" \
    | LC_ALL=C sort -u \
    > "$OUT/source_found.ids.txt"

# Find requested IDs that were not extracted from the source.
LC_ALL=C comm -23 \
    "$OUT/refchanged.ids.txt" \
    "$OUT/source_found.ids.txt" \
    > "$OUT/source_missing.ids.txt"

if [[ -s "$OUT/source_missing.ids.txt" ]]; then
    printf '\nError: some REF-changed IDs have no source evidence:\n' >&2
    cat "$OUT/source_missing.ids.txt" >&2
    exit 1
fi

printf '\nAll REF-changed IDs have source evidence.\n'

# Summarize unique IDs represented by each origin.
printf '\nUnique source IDs represented in the evidence table:\n'

awk -F '\t' '
    NR > 1 {
        print $1 "\t" $2
    }
' "$OUT/refchanged.evidence.tsv" \
    | LC_ALL=C sort -u \
    | awk -F '\t' '
        { counts[$1]++ }
        END {
            for (origin in counts)
                print origin "\t" counts[origin]
        }
    ' \
    | LC_ALL=C sort




    