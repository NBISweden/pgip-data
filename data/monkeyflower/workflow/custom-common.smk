def custom_multiqc_roi_input(wildcards, callmode, callsets):
    return expand(
        f"{wildcards.roi}/vcftools-{{callmode}}/{{callsets}}.allsites.subset.{{stat}}.{{ext}}",
        stat=VCFTOOLS_STATS.keys(),
        ext=["tab", "png"],
        callmode=callmode
    )
