# Custom rules for monkeyflower data


rule custom_all:
    input:
        ["sampleinfo.csv"],


rule make_sampleinfo:
    """Join SraRunInfo.csv and samples.csv (which must be manually generated)"""
    output:
        sampleinfo="sampleinfo.csv",
    input:
        runinfo="SraRunInfo.csv",
        samples="samples.csv",
    conda: "envs/environment.yml",
    benchmark:
        "benchmarks/sampleinfo.csv.benchmark.txt"
    log:
        "logs/sampleinfo.csv.log",
    threads: 1
    shell:
        """
        cat SraRunInfo.csv | cut -d "," -f 1,29,30 |
            csvtk mutate - -n Sample -f SampleName -p "^[A-Z]+-(.+)" |
            csvtk mutate - -R -n Sample -f Sample -p "^([A-Z].+)_[0-9]+" |
            csvtk mutate - -R -n Sample -f Sample -p "^([^TK1][A-Z]+)[0-9]*-[0-9]+" |
            csvtk replace -f Sample -p "CLV" -r "CLV_" |
            csvtk replace -f Sample -p "GH1" -r "CLV_GH" > sraruninfo.csv
        csvtk join samples.csv sraruninfo.csv -f "Sample" > sampleinfo.csv
        """
