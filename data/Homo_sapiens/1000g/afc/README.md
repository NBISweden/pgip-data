# 1000 genomes allele frequency data

<!-- markdownlint-disable MD013 -->

Allele frequency data was downloaded using [Ensembls Allele Frequency
Calculator](https://grch37.ensembl.org/Homo_sapiens/Tools/AlleleFrequency?db=core)
for region 1:10000001-12500000 for populations CEU, CHB, and YRI. The
output tsv files were manually cleaned to keep only biallelic SNPs:

<!-- markdownlint-enable MD013 -->

    for f in *.tsv; do
        awk '{if(length($4) == 1 && length($5) == 1 || $2 == "POS") print}' $f |\
        gzip -v > $f.gz;
    done
