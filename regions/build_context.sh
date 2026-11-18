#!/usr/bin/env bash
set -euo pipefail

# gives the location of where TSV file is
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
FLAGS_FILE="${SCRIPT_DIR}/giab_context_flags.tsv"


#three inputs, only run when exactly three inputs and in what order
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <giab_version> <hg38.fa.fai> <outdir>" >&2
    exit 1
fi

GIAB_VERSION="$1"
FAI_FILE="$2"
OUT_DIR="$3"

# check the two input files exist and are readable

if [[ ! -r "$FLAGS_FILE" ]]; then
    echo "Error: cannot read flags table: $FLAGS_FILE" >&2
    exit 1
fi

if [[ ! -r "$FAI_FILE" ]]; then
    echo "Error: cannot read FASTA index: $FAI_FILE" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

# checks that the commands: wget, md5sum, bedtools, gzip, bgzip, tabix, awk, sort, cut and grep are available 
for tool in wget bedtools bgzip tabix awk sort cut grep gzip md5sum; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Error: required tool is not installed: $tool" >&2
        exit 1
    fi
done

# where to download the giab files from and where to store them
BASE_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/${GIAB_VERSION}/GRCh38@all"

DOWNLOAD_DIR="${OUT_DIR}/downloads"
FLAGS_DIR="${OUT_DIR}/flags"

mkdir -p "$DOWNLOAD_DIR" "$FLAGS_DIR"

# loop to read tsv file, one row at a time, donwloading the files
while IFS=$'\t' read -r priority flag file; do
    if [[ "$priority" == "priority" ]]; then
        continue
    fi

    url="${BASE_URL}/${file}"
    destination="${DOWNLOAD_DIR}/${file}"

    mkdir -p "$(dirname -- "$destination")"

    printf 'Downloading: %s\n' "$file"
    wget -O "$destination" "$url"
done < "$FLAGS_FILE"


# a checksum to check for incomplete downloads
CHECKSUM_FILE="${DOWNLOAD_DIR}/GRCh38-genome-stratifications-md5s.txt"

printf 'Downloading the checksum list\n'
wget -O "$CHECKSUM_FILE" \
    "${BASE_URL}/GRCh38-genome-stratifications-md5s.txt"

SELECTED_CHECKSUMS="${DOWNLOAD_DIR}/selected-md5s.txt"
: > "$SELECTED_CHECKSUMS"

while IFS=$'\t' read -r priority flag file; do
    if [[ "$priority" == "priority" ]]; then
        continue
    fi

    awk -v target="$file" '
        $2 == target {
            print
            found++
        }
        END {
            if (found != 1) {
                print "Expected one checksum for: " target > "/dev/stderr"
                exit 1
            }
        }
    ' "$CHECKSUM_FILE" >> "$SELECTED_CHECKSUMS"

done < "$FLAGS_FILE"

# actual verification of the downloaded files
printf 'Verifying the 19 downloaded files\n'

(
    cd -- "$DOWNLOAD_DIR"
    md5sum -c selected-md5s.txt
)


# now we are creatubg the 10 flag priority order
FLAG_ORDER="${OUT_DIR}/flag-order.tsv"

awk -F '\t' 'NR > 1 { print $1 "\t" $2 }' "$FLAGS_FILE" \
    | LC_ALL=C sort -k1,1n -k2,2 -u \
    > "$FLAG_ORDER"


# combine downloaded files into one BED file per flag
while IFS=$'\t' read -r priority flag; do
    printf 'Building flag: %s\n' "$flag"

    {
        while IFS=$'\t' read -r file_priority file_flag file; do
            if [[ "$file_flag" == "$flag" ]]; then
                gzip -cd -- "${DOWNLOAD_DIR}/${file}"
            fi
        done < "$FLAGS_FILE"
    } \
        | awk 'BEGIN { OFS="\t" }
            $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {
                print $1, $2, $3
            }' \
        | LC_ALL=C sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        > "${FLAGS_DIR}/${flag}.bed"

done < "$FLAG_ORDER"


# now we are applying the priority order
CLAIMED_BED="${OUT_DIR}/claimed.bed"
CONTEXT_BED="${OUT_DIR}/context.unsorted.bed"

: > "$CLAIMED_BED"
: > "$CONTEXT_BED"

while IFS=$'\t' read -r priority flag; do
    printf 'Applying priority %s: %s\n' "$priority" "$flag"

    {
        if [[ -s "$CLAIMED_BED" ]]; then
            bedtools subtract \
                -a "${FLAGS_DIR}/${flag}.bed" \
                -b "$CLAIMED_BED"
        else
            cat "${FLAGS_DIR}/${flag}.bed"
        fi
    } \
        | awk -v label="$flag" 'BEGIN { OFS="\t" }
            { print $1, $2, $3, label }' \
        >> "$CONTEXT_BED"

    cat "$CLAIMED_BED" "${FLAGS_DIR}/${flag}.bed" \
        | LC_ALL=C sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        > "${OUT_DIR}/claimed.tmp.bed"

    mv -- "${OUT_DIR}/claimed.tmp.bed" "$CLAIMED_BED"

done < "$FLAG_ORDER"


# now we sort, compress and index the final context file
CONTEXT_GZ="${OUT_DIR}/giab_${GIAB_VERSION}_context.hg38.bed.gz"

LC_ALL=C sort -k1,1 -k2,2n "$CONTEXT_BED" \
    | bgzip -c \
    > "$CONTEXT_GZ"

tabix -f -p bed "$CONTEXT_GZ"

printf 'Checking the context map for overlaps\n'

gzip -cd -- "$CONTEXT_GZ" \
    | awk '
        {
            if ($1 == previous_chrom && $2 < previous_end) {
                print "Error: overlapping context intervals at " \
                    $1 ":" $2 "-" $3 > "/dev/stderr"
                exit 1
            }

            previous_chrom = $1
            previous_end = $3
        }
        END {
            if (NR == 0) {
                print "Error: the context map is empty" > "/dev/stderr"
                exit 1
            }
        }
    '

printf 'Overlap check passed\n'

# calculation of total reference length
PRIMARY_FAI="${OUT_DIR}/primary-contigs.tsv"

awk 'BEGIN { OFS="\t" }
    $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {
        print $1, $2
    }
' "$FAI_FILE" > "$PRIMARY_FAI"

# check ref index contains all 24 primary chromosomes
PRIMARY_COUNT="$(awk 'END { print NR }' "$PRIMARY_FAI")"

if [[ "$PRIMARY_COUNT" -ne 24 ]]; then
    echo "Error: expected 24 primary chromosomes, found $PRIMARY_COUNT" >&2
    echo "Use the full chr-prefixed hg38 FASTA index." >&2
    exit 1
fi

GENOME_BP="$(
    awk '{ total += $2 } END { printf "%.0f\n", total }' "$PRIMARY_FAI"
)"

if [[ "$GENOME_BP" -le 0 ]]; then
    echo "Error: no primary hg38 chromosome lengths found" >&2
    exit 1
fi

printf 'Primary reference length: %s bp\n' "$GENOME_BP"


# calculation of base pairs and genome fraction for each context, incl unique
SIZES_FILE="${OUT_DIR}/giab_${GIAB_VERSION}_context_sizes.tsv"

awk -v genome="$GENOME_BP" '
    FILENAME == ARGV[1] {
        labels[++n] = $2
        next
    }

    {
        bp[$4] += $3 - $2
    }

    END {
        for (i = 1; i <= n; i++) {
            assigned += bp[labels[i]]
        }

        if (assigned > genome) {
            print "Error: context bases exceed reference length" > "/dev/stderr"
            exit 1
        }

        print "context\tbp\tfraction"

        for (i = 1; i <= n; i++) {
            label = labels[i]
            printf "%s\t%.0f\t%.10f\n", \
                label, bp[label], bp[label] / genome
        }

        remaining = genome - assigned
        printf "unique\t%.0f\t%.10f\n", \
            remaining, remaining / genome
    }
' "$FLAG_ORDER" "$CONTEXT_BED" > "${SIZES_FILE}.tmp"

mv -- "${SIZES_FILE}.tmp" "$SIZES_FILE"

# check if applying priority hasnt lost any covered bases
UNION_BP="$(
    awk '{ total += $3 - $2 } END { printf "%.0f\n", total }' \
        "$CLAIMED_BED"
)"

LABELLED_BP="$(
    awk '{ total += $3 - $2 } END { printf "%.0f\n", total }' \
        "$CONTEXT_BED"
)"

if [[ "$UNION_BP" -ne "$LABELLED_BP" ]]; then
    echo "Error: labelled bases do not equal the combined region coverage" >&2
    exit 1
fi

printf 'Coverage check passed: %s labelled bases\n' "$LABELLED_BP"

