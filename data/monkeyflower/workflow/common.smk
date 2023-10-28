import re
import itertools
import copy

VCFTOOLS_STATS = {
    "het": "het",
    "idepth": "depth",
    "ldepth": "site-depth",
    "lqual": "site-quality",
    "frq": "freq2",
    "imiss": "missing-indv",
    "lmiss": "missing-site",
}


def get_envmodules(envmodules):
    retmodules = []
    if not isinstance(envmodules, list):
        envmodules = [envmodules]
    try:
        retmodules = copy.deepcopy(config["envmodules"]["__site__"])
    except KeyError as e:
        retmodules = []
    for mod in envmodules:
        try:
            retmodules.extend(config["envmodules"][mod])
        except KeyError as e:
            pass
    return retmodules


def make_roi_bed_roi(wildcards):
    res = str()
    try:
        res = "\n".join(
            [
                "\t".join(re.split(":|-", x))
                for x in config["output"][wildcards.roi]["roi"]
            ]
        )
    except:
        pass
    return res


def get_read_group(wildcards):
    info = get_samplealias_dict(wildcards)
    run = info["srrun"]
    sample = info["samplename"]
    rg = f'-R "@RG\\tID:{run}\\tSM:{sample}\\tPL:ILLUMINA"'
    return rg


def get_roi_read_group(wildcards):
    info = get_samplealias_dict(wildcards)
    run = info["srrun"]
    sample = info["sample"]
    rg = f'-R "@RG\\tID:{run}\\tSM:{sample}\\tPL:ILLUMINA"'
    return rg


def get_samplealias_dict(wildcards):
    samplealias = wildcards.samplealias
    samplename = sampleinfo.set_index("SampleAlias").loc[samplealias].SampleName
    srrun = sampleinfo.set_index("SampleAlias").loc[samplealias].Run
    info = {"srrun": srrun, "samplename": samplename, "sample": samplealias}
    return info


def get_sequence_fai(wildcards):
    if wildcards.roi == "":
        path = config["reference"]
    else:
        path = os.path.join(wildcards.roi, config["reference"])
    return f"{path}.fai"


def get_sequence_dictionary(wildcards):
    if wildcards.roi == "":
        path = config["reference"]
    else:
        path = os.path.join(wildcards.roi, config["reference"])
    return re.sub(".fasta$", ".dict", path)


def subset_bam_for_roi_input(wildcards, *, ext=""):
    samplealias = wildcards.samplealias
    srrun = sampleinfo.set_index("SampleAlias").loc[samplealias].Run
    return f"{samplealias}/{srrun}.sort.md.bam{ext}"


def gatk_raw_or_bqsr_variant_filtration_options(wildcards):
    if wildcards.raw == ".raw":
        options = [
            "--filter-name",
            "FisherStrand",
            "--filter",
            "'FS > 50.0'",
            "--filter-name",
            "QualByDepth",
            "--filter",
            "'QD < 4.0'",
            "--filter-name",
            "MappingQuality",
            "--filter",
            "'MQ < 50.0'",
        ]
    else:
        options = [
            "--filter-name",
            "FisherStrand",
            "--filter",
            "'FS > 60.0'",
            "--filter-name",
            "QualByDepth",
            "--filter",
            "'QD < 2.0'",
            "--filter-name",
            "MappingQuality",
            "--filter",
            "'MQ < 40.0'",
        ]

    return " ".join(options)


def make_roi_save_script_input(wildcards):
    roi = config["output"][wildcards.roi]
    data = [os.path.join(wildcards.roi, "all.allsites.vcf.gz")]
    if roi["ubam"]:
        data.extend(
            expand(
                f"{wildcards.roi}/{{srrun}}/{{srrun}}.{{samplealias}}.unmapped.bam",
                zip,
                srrun=sampleinfo.Run.values,
                samplealias=sampleinfo.SampleAlias.values,
            )
        )
        data.extend(
            expand(
                f"{wildcards.roi}/{{srrun}}/{{srrun}}.{{samplealias}}.unmapped.bam.bai",
                zip,
                srrun=sampleinfo.Run.values,
                samplealias=sampleinfo.SampleAlias.values,
            )
        )
    if roi["fastq"]:
        data.extend(
            expand(
                f"{wildcards.roi}/{{srrun}}/{{srrun}}.{{samplealias}}_R1.fastq.gz",
                zip,
                srrun=sampleinfo.Run.values,
                samplealias=sampleinfo.SampleAlias.values,
            )
        )
        data.extend(
            expand(
                f"{wildcards.roi}/{{srrun}}/{{srrun}}.{{samplealias}}_R2.fastq.gz",
                zip,
                srrun=sampleinfo.Run.values,
                samplealias=sampleinfo.SampleAlias.values,
            )
        )
    if roi["reference"]:
        data.append(os.path.join(wildcards.roi, config["reference"]))
    return data


def gatk_combine_gvcfs_input(wildcards):
    fmt = (
        f"{wildcards.roi}{wildcards.sep}{{samplealias}}/{{samplealias}}"
        ".sort.dup.recal.hc.g.vcf.gz"
    )
    try:
        samples = config["output"][wildcards.roi]["callset"][wildcards.callset][
            "samples"
        ]
        info = sampleinfo[sampleinfo.SampleAlias.isin(samples)]
    except:
        info = sampleinfo
    vcf = expand(fmt, samplealias=info.SampleAlias.values)
    return vcf


def multiqc_roi_input(wildcards):
    def _expand_fmt(fmt, qc):
        fmt = fmt.format(qc=qc)
        return expand(
            fmt,
            samplealias=sampleinfo.SampleAlias.values,
        )

    fmt = f"{wildcards.roi}/{{qc}}/{{{{samplealias}}}}/{{{{samplealias}}}}"
    results = dict()
    results["fastqc"] = _expand_fmt(f"{fmt}_R1_fastqc/summary.txt", "fastqc")
    results["qualimap"] = _expand_fmt(
        f"{fmt}.sort_stats/genome_results.txt", "qualimap"
    )
    results["markdups"] = _expand_fmt(f"{fmt}.sort.dup.dup_metrics.txt", "markdup")
    fmt = f"{wildcards.roi}/vcftools/all.allsites.subset.{{stat}}"
    results["vcftools"] = expand(fmt, stat=VCFTOOLS_STATS.keys())
    if len(custom_all) != 0:
        results["vcftools.custom"] = custom_multiqc_roi_input(wildcards)
    return itertools.chain(*results.values())


def repeat_masker_dfam(wildcards):
    default = "Libraries/Dfam.h5"
    dfam = os.path.join(os.environ.get("REPEATMASKER_LIBDIR"), default)
    if os.path.exists(dfam):
        return dfam
    return default


def vcftools_stats(wildcards):
    return f"--{VCFTOOLS_STATS[wildcards.stat]}"


def vcftools_stats_options(wildcards):
    default = ""
    options = {
        "frq": "--max-alleles 2",
    }
    if wildcards.stat not in options.keys():
        return default
    return options[wildcards.stat]


def csvtk_compile_vcftools_stats_options(wildcards):
    d = {
        "idepth": "-f MEAN_DEPTH:mean",
        "ldepth": "-f 3:min,3:q1,3:median,3:mean,3:q3,3:max",
        "lqual": "-f QUAL:min,QUAL:q1,QUAL:median,QUAL:mean,QUAL:q3,QUAL:max",
        "lmiss": "-f 6:min,6:q1,6:median,6:mean,6:q3,6:max",
    }
    if wildcards.stat not in d.keys():
        return "-f 1"
    return d[wildcards.stat]


def csvtk_plot_vcftools_stats(wildcards):
    if wildcards.stat == "ldepth":
        cmd = [
            "csvtk summary -t -g SUM_DEPTH -f POS:count -w 0",
            "csvtk sort -t -k 1:n",
            'csvtk plot line -t - -x SUM_DEPTH -y POS:count --point-size 0.01 --xlab "Depth of coverage (X)" --ylab "Genome coverage (bp)" --width 9.0 --height 3.5',
        ]
    elif wildcards.stat == "lqual":
        cmd = [
            'csvtk filter -t -f "QUAL>0" -f "QUAL<1000"',
            "csvtk summary -t -g QUAL -f POS:count -w 0 -",
            "csvtk sort -t -k 1:n",
            'csvtk plot hist -t --bins 100 - --xlab "Quality value" --ylab "Count" --width 9.0 --height 3.5',
        ]
    elif wildcards.stat == "frq":
        cmd = [
            "csvtk fix -t -T",
            'csvtk filter -t -f "{FREQ}!=0" -f "{FREQ}!=1" -',
            'csvtk mutate2 -t -n maf -e \'${5} > ${6} ? "${6}" : "${5}" \' -',
            'csvtk plot hist -t --bins 20 -f maf - --xlab "Minor allele frequency" --ylab "Count" --width 9.0 --height 3.5',
        ]
    elif wildcards.stat == "imiss":
        cmd = ["csvtk plot hist -t --x-min 0 -f F_MISS"]
    elif wildcards.stat == "lmiss":
        cmd = ["csvtk plot hist --bins 20 -t -f F_MISS"]
    elif wildcards.stat == "het":
        cmd = ["csvtk plot hist --bins 20 -t -f F"]
    elif wildcards.stat == "idepth":
        cmd = ["csvtk plot hist --bins 20 -t -f MEAN_DEPTH"]

    return "|".join(cmd)


def csvtk_plot_vcftools_stats_options(wildcards):
    res = str()
    try:
        callset_sites = f"{wildcards.callset}.{wildcards.sites}"
        res = config["output"][wildcards.roi]["plot"][callset_sites][wildcards.stat]
    except KeyError as e:
        pass
    return res


def vcftools_filter_options(wildcards):
    options = []
    try:
        callset_sites = f"{wildcards.callset}.{wildcards.sites}"
        obj = config["output"][wildcards.roi]["filter"][callset_sites]
    except KeyError as e:
        return options
    if "miss" in obj.keys():
        options.append(f"--max-missing {obj['miss']}")
    if "qual" in obj.keys():
        options.append(f"--minQ {obj['qual']}")
    if "min_depth" in obj.keys():
        options.append(f"--min-meanDP {obj['min_depth']}")
        options.append(f"--minDP {obj['min_depth']}")
    if "max_depth" in obj.keys():
        options.append(f"--max-meanDP {obj['max_depth']}")
        options.append(f"--maxDP {obj['max_depth']}")
    if "options" in obj.keys():
        options.append(obj["options"])
    return " ".join(options)
