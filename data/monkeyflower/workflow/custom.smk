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
if "output" in config.keys():
    for roi, obj in config["output"].items():
        results.append(os.path.join(f"{roi}", "redyellow.allsites.vcf.gz.tbi"))


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


rule custom_gatk_combine_gvcfs:
    """Run GATK CombineGVCFs"""
    output:
        vcf="{roi}{sep}{label}combine.g.vcf.gz",
        tbi="{roi}{sep}{label}combine.g.vcf.gz.tbi",
    input:
        vcf=custom_gatk_combine_gvcfs_input,
        ref=os.path.join("{roi}", config["reference"]),
    wildcard_constraints:
        label="(redyellow.)",
    params:
        vcf=lambda wildcards, input: " ".join([f"-V {x}" for x in input.vcf]),
    conda:
        "envs/environment.yml"
    benchmark:
        "benchmarks/gatk_combine_gvcfs/{roi}{sep}{label}combined.g.vcf.gz.benchmark.txt"
    log:
        "logs/gatk_combine_gvcfs/{roi}{sep}{label}combined.g.vcf.gz.log",
    threads: 1
    shell:
        """
        gatk CombineGVCFs -OVI true --output {output.vcf} --reference {input.ref} {params.vcf} > {log} 2>&1
        """


rule custom_gatk_genotype_gvcfs:
    """GATK GenotypeGVCFs"""
    output:
        vcf="{roi}{sep}{label}allsites.vcf.gz",
        tbi="{roi}{sep}{label}allsites.vcf.gz.tbi",
    input:
        vcf="{roi}{sep}{label}combine.g.vcf.gz",
        ref=os.path.join("{roi}", config["reference"]),
    wildcard_constraints:
        label="(redyellow.)",
    conda:
        "envs/environment.yml"
    benchmark:
        "benchmarks/gatk_genotype_gvcfs/{roi}{sep}{label}allsites.vcf.gz.benchmark.txt"
    log:
        "logs/gatk_genotype_gvcfs/{roi}{sep}{label}allsites.log",
    threads: 1
    shell:
        """
        gatk GenotypeGVCFs -OVI true -R {input.ref} -V {input.vcf} -O {output.vcf} --all-sites > {log} 2>&1
        """
