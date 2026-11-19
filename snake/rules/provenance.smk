# ─── Run provenance ──────────────────────────────────────────────────────────
# Every file the results depend on, by sha256, and the tool versions.
_PROV_FILES = {
    "input_vcf": config["input_vcf"],
    **({"native_vcf": config["native_vcf"]} if COMPARE else {}),
    **{f"ref.{k}": v for k, v in config["ref"].items()},
    **{f"regions.{k}": v for k, v in config["regions"].items()},
    **{f"config.{i}": p for i, p in enumerate(workflow.configfiles)},
}

rule provenance:
    """
    provenance.tsv: kind (sha256 | version), name, value, path.
    Python tools by package version: CrossMap 0.7.3 reports 0.7.0 itself.
    """
    output:
        f"{RESULTS}/provenance.tsv"
    params:
        files = " ".join(f"{k}={v}" for k, v in _PROV_FILES.items())
    shell:
        r"""
        (
        printf 'kind\tname\tvalue\tpath\n'
        for kv in {params.files}; do
            f=${{kv#*=}}
            printf 'sha256\t%s\t%s\t%s\n' "${{kv%%=*}}" "$(sha256sum "$f" | cut -d' ' -f1)" "$f"
        done
        printf 'version\t%s\t%s\t.\n' \
            bcftools  "$(bcftools --version | head -1 | cut -d' ' -f2)" \
            picard    "$(picard LiftoverVcf --version 2>&1 | tail -1 | cut -d: -f2)" \
            transanno "$(transanno --version | cut -d' ' -f2)" \
            snakemake "$(snakemake --version)"
        python3 -c 'import importlib.metadata as m
for p in ("CrossMap", "svdb", "pysam"): print(f"version\t{{p}}\t{{m.version(p)}}\t.")'
        ) > {output}
        """
