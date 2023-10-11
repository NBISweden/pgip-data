def get_read_group(wildcards):
    run=wildcards.srrid
    sample = sraruninfo[sraruninfo["Run"] == run].Sample.values[0]
    rg = f'-R "@RG\\tID:{wildcards.srrid}\\tSM:{sample}\\tPL:ILLUMINA"',
    return rg


def gatk_raw_or_bqsr_variant_filtration_options(wildcards):
    if wildcards.raw == ".raw":
        options = [
            "--filter-name", "FisherStrand",
            "--filter", "'FS > 50.0'",
            "--filter-name", "QualByDepth",
            "--filter", "'QD < 4.0'",
            "--filter-name", "MappingQuality",
            "--filter", "'MQ < 50.0'"
        ]
    else:
        options = [
            "--filter-name", "FisherStrand",
            "--filter", "'FS > 60.0'",
            "--filter-name", "QualByDepth",
            "--filter", "'QD < 2.0'",
            "--filter-name", "MappingQuality",
            "--filter", "'MQ < 40.0'"
        ]

    return " ".join(options)

def make_roi_save_script_input(wildcards):
    roi = config["output"][wildcards.roi]
    data = [os.path.join(wildcards.roi, "allsites.vcf.gz")]
    if roi["ubam"]:
        data.extend(
            expand(f"{wildcards.roi}/{{srrid}}/{{srrid}}.{{sampleid}}.unmapped.bam",
                   zip,
                   srrid=sraruninfo.Run.values,
                   sampleid=sraruninfo.SampleName.values)
        )
        data.extend(
            expand(f"{wildcards.roi}/{{srrid}}/{{srrid}}.{{sampleid}}.unmapped.bam.bai",
                   zip,
                   srrid=sraruninfo.Run.values,
                   sampleid=sraruninfo.SampleName.values)
        )
    if roi["fastq"]:
        data.extend(
            expand(f"{wildcards.roi}/{{srrid}}/{{srrid}}.{{sampleid}}_R1.fastq.gz",
                   zip,
                   srrid=sraruninfo.Run.values,
                   sampleid=sraruninfo.SampleName.values)
        )
        data.extend(
            expand(f"{wildcards.roi}/{{srrid}}/{{srrid}}.{{sampleid}}_R2.fastq.gz",
                   zip,
                   srrid=sraruninfo.Run.values,
                   sampleid=sraruninfo.SampleName.values)
        )
    if roi["reference"]:
        data.append(os.path.join(wildcards.roi, config["reference"]))
    return data

