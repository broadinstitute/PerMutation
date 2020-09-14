# PerMutation
Code for Tn-Seq analysis of _E. faecalis_ MMH594. Check out the wiki more details:

1. [Processing, Mapping, Replication Bias Correction, and Normalization Across Replicates](https://github.com/broadinstitute/PerMutation/wiki/1.-Processing,-Mapping,-Replication-Bias-Correction,-and-Normalization-Across-Replicates)
2. [Assessing the Fitness Cost of Genes For Growth on Nutrient Rich Media](https://github.com/broadinstitute/PerMutation/wiki/2.-Assessing-the-Fitness-Cost-of-Genes-For-Growth-on-Nutrient-Rich-Media)
3. [Predicting Whether Genes Contribute to Antibiotic Resistance](https://github.com/broadinstitute/PerMutation/wiki/3.-Predicting-Whether-Genes-Contribute-to-Antibiotic-Resistance)

All programs/scripts written in Python-3, no real installation required. Python-3 libraries used include:

* biopython 
* numpy
* scipy
* pysam

For processing and mapping of sequencing data to the reference genome the following third-party tools are necessary:

* trimmomatic
* samtools
* bowtie2
* fastx_toolkit

Progams should work if the user creates and uses a conda environment from the provided yml file.
