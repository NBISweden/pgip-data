# Custom rules for monkeyflower data
# Functions for custom.smk
include: "custom-common.smk"


try:
    sampleinfo = pd.read_csv(os.path.join(config["bioproject"], "sampleinfo.csv"))
except Exception as e:
    logging.error(
        "CUSTOM: No sampleinfo.csv available yet! Rerun workflow once it has been downloaded"
    )
    pass


results = ["sampleinfo.csv"]


rule custom_all:
    input:
        results,


rule make_sampleinfo:
    """Join SraRunInfo.csv and samples.csv (which must be manually generated)"""
    output:
        sampleinfo="sampleinfo.csv",
    input:
        runinfo="SraRunInfo.csv",
        samples="samples.csv",
    envmodules:
        get_envmodules("csvtk"),
    conda:
        "envs/environment.yml"
    benchmark:
        "benchmarks/sampleinfo.csv.benchmark.txt"
    log:
        "logs/sampleinfo.csv.log",
    threads: 1
    shell:
        """
        csvtk cut -f "Sample,Run,ScientificName,SampleName" SraRunInfo.csv |
        csvtk mutate - -n AuthorSample  -f SampleName -p "^([A-Z]+-.+[A-Z])[0-9]+-2016" |
        csvtk mutate - -n SampleAlias -f AuthorSample |
        csvtk replace -f SampleAlias,AuthorSample -p "(^[^A].+)_[0-9]+$" -r "\$1" |
        csvtk replace -f SampleAlias,AuthorSample -p "(^A[A-Z]+-T[0-9]+)_[0-9]+$" -r "\$1" |
        csvtk replace -f SampleAlias,AuthorSample -p "CLV" -r "CLV_" |
        csvtk replace -f SampleAlias,AuthorSample -p "-GH1" -r "GH" |
        csvtk replace -f AuthorSample -p "^[A-Z]+-" -r "" |
        csvtk replace -f SampleAlias -p "PUN-(BCRD|INJ|LO|PCT|POTR)" -r "PUN-Y-\$1" |
        csvtk replace -f SampleAlias -p "PUN-(ELF|JMC|LH|MT|UCSD)" -r "PUN-R-\$1" > sraruninfo.csv

        csvtk join sraruninfo.csv samples.csv -f "AuthorSample;Sample" > sampleinfo.csv
        """
