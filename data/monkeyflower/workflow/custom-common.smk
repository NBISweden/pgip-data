def custom_gatk_combine_gvcfs_input(wildcards):
    yellow = ["PUN-Y-BCRD", "PUN-Y-INJ", "PUN-Y-LO",
              "PUN-Y-PCT", "PUN-Y-POTR"]
    red = ["PUN-R-ELF", "PUN-R-JMC", "PUN-R-LH",
           "PUN-R-MT", "PUN-R-UCSD"]

    fmt = (
        f"{wildcards.roi}{wildcards.sep}{{samplealias}}/{{srrun}}"
        ".sort.dup.recal.hc.g.vcf.gz"
    )
    if wildcards.label == "redyellow.":
        i = sampleinfo.ScientificName.str.find("puniceus") > 0
        df = sampleinfo[i]
        vcf = expand(fmt, zip, samplealias=df.SampleAlias.values,
                     srrun=df.Run.values)
    return vcf

