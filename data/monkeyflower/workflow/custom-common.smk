def custom_multiqc_roi_input(wildcards):
    return expand(
        f"{wildcards.roi}/vcftools/redyellow.allsites.subset.{{stat}}.{{ext}}",
        stat=VCFTOOLS_STATS.keys(),
        ext=["tab", "png"],
    )
