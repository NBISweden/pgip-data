import re
import itertools

def get_read_group(wildcards):
    info = get_srrun_dict(wildcards)
    run = info["srrun"]
    sample = info["samplename"]
    rg = f'-R "@RG\\tID:{run}\\tSM:{sample}\\tPL:ILLUMINA"'
    return rg


def get_roi_read_group(wildcards):
    info = get_srrun_dict(wildcards)
    run = info["srrun"]
    sample = info["sample"]
    rg = f'-R "@RG\\tID:{run}\\tSM:{sample}\\tPL:ILLUMINA"'
    return rg


def get_srrun_dict(wildcards):
    srrun = wildcards.srrun
    samplename = sampleinfo[sampleinfo["Run"] == srrun].SampleName.values[0]
    sample = sampleinfo[sampleinfo["Run"] == srrun].SampleAlias.values[0]

    info = {'srrun': srrun, 'samplename': samplename, 'sample': sample}
    return info


def get_sequence_fai(wildcards):
    if wildcards.roi == "":
        path = config["reference"]
    else:
        path = os.path.join(wildcards.roi, config["reference"])
    return path + ".fai"


def get_sequence_dictionary(wildcards):
    if wildcards.roi == "":
        path = config["reference"]
    else:
        path = os.path.join(wildcards.roi, config["reference"])
    return re.sub(".fasta$", ".dict", path)


def subset_bam_for_roi_input(wildcards, *, ext=""):
    srsamplename = sampleinfo[sampleinfo["SampleAlias"] == wildcards.samplealias].SampleName.values[0]
    return f"{srsamplename}/{wildcards.srrun}.sort.md.bam{ext}"


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
            expand(f"{wildcards.roi}/{{srrun}}/{{srrun}}.{{samplealias}}.unmapped.bam",
                   zip,
                   srrun=sampleinfo.Run.values,
                   samplealias=sampleinfo.SampleAlias.values)
        )
        data.extend(
            expand(f"{wildcards.roi}/{{srrun}}/{{srrun}}.{{samplealias}}.unmapped.bam.bai",
                   zip,
                   srrun=sampleinfo.Run.values,
                   samplealias=sampleinfo.SampleAlias.values)
        )
    if roi["fastq"]:
        data.extend(
            expand(f"{wildcards.roi}/{{srrun}}/{{srrun}}.{{samplealias}}_R1.fastq.gz",
                   zip,
                   srrun=sampleinfo.Run.values,
                   samplealias=sampleinfo.SampleAlias.values)
        )
        data.extend(
            expand(f"{wildcards.roi}/{{srrun}}/{{srrun}}.{{samplealias}}_R2.fastq.gz",
                   zip,
                   srrun=sampleinfo.Run.values,
                   samplealias=sampleinfo.SampleAlias.values)
        )
    if roi["reference"]:
        data.append(os.path.join(wildcards.roi, config["reference"]))
    return data


def gatk_combine_gvcfs_input(wildcards):
    fmt = (
        f"{wildcards.roi}{wildcards.sep}{{samplealias}}/{{srrun}}"
        ".sort.dup.recal.hc.g.vcf.gz"
    )
    vcf = expand(fmt, zip,
                 samplealias=sampleinfo.SampleAlias.values,
                 srrun=sampleinfo.Run.values)
    return vcf


def multiqc_roi_input(wildcards):
    def _expand_fmt(fmt, qc):
        fmt = fmt.format(qc=qc)
        return expand(fmt, zip,
                      samplealias=sampleinfo.SampleAlias.values,
                      srrun=sampleinfo.Run.values)

    fmt = f"{wildcards.roi}/{{qc}}/{{{{samplealias}}}}/{{{{srrun}}}}"
    results = dict()
    results["fastqc"] = _expand_fmt(fmt + "_R1_fastqc/summary.txt",
                                    "fastqc")
    results["qualimap"] = _expand_fmt(fmt + ".sort_stats/genome_results.txt",
                                      "qualimap")
    results["markdups"] = _expand_fmt(fmt + ".sort.dup.dup_metrics.txt",
                                      "markdup")
    return itertools.chain(*results.values())
