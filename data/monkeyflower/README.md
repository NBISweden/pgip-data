# Monkeyflower data set

Reference sequence and unmapped bam files from study "Widespread
selection and gene flow shape the genomic landscape during a radiation
of monkeyflowers" by Stankowski et al, PLoS Biology, 2019.

samples.csv was manually created from
https://doi.org/10.1371/journal.pbio.3000391.s003

## Workflow

Bash script contains commands to

- download reference sequence from mimubase
- download SRR files from ncbi read archive
- map reads to reference

Since exercises concentrate on a region of interest consisting of 3Mbp
around gene MaMyb2, the script also contains commands to

- extract reads in ROI in ubam format
- remap reads to ROI
- perform variant calling in ROI following Stankowski et al

A small subset of the entries in the ubam files have been committed to
this repository to serve as test data set to generate the web site.
The full ROI files are tarred and uploaded to a server accessible to
students.
