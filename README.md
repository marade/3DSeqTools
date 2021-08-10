# Software tools for <i>Genome-wide protein–DNA interaction site mapping using a double strand DNA-specific cytosine deaminase</i>
![logo](/title.png)
## Guide for Use
### Prerequisites
This software has been tested on the Linux operating system. It may be possible to adapt it for other operating systems. The following dependencies were used for the paper, though other versions may work:
* Python 3.8.5
* Python libraries:
  * Biopython 1.78
  * Pandas 1.3.0
  * colorama 0.4.3
  * PySAM 0.16.0.1
* HTStream 1.3.0 
* Minimap2 2.17-r974-dirty
* SAMTools 1.10
* Picard Tools 2.18.25

### Do Sequencing and Generate Fastq Files
We assume you have generated your sequencing data in roughly the manner described in the paper, using Illumina paired-end sequencing. We provide some example files for testing below.
### Prepare Fastq Files
This pipeline assumes your paired-end Fastq files are named like so:

    sampleX_1.fastq.gz sampleX_2.fastq.gz
    sampleY_1.fastq.gz sampleY_2.fastq.gz

### Run the Pipeline
The following pipeline works on a vanilla Ubuntu 20.04 installation with the prerequisites installed. We were able to process these data in the cloud with an AWS EC2 C5.2xlarge instance (8 compute cores; 16GB RAM).

    $ git clone https://github.com/marade/3DSeqTools.git
    $ cd 3DSeqTools && mkdir fastq && cd fastq
    # Note you are downloading nearly 4GB of fastq data here:
    $ for x in {1..4}; do for y in {1..2}; do wget https://3d-seq-01.s3.us-west-2.amazonaws.com/fastq/gcsR-dddA-Para-dddI-delta-ung-0Ara-P3-rep${x}_${y}.fastq.gz ; done ; done
    $ cd ../
    $ python3 AlignReads -n 8 fastq NCBI/Pseudomonas_aeruginosa_PAO1_GCF_000006765.1 tmf
    $ python3 Process3DSeq -n 3 tmf/align NCBI/Pseudomonas_aeruginosa_PAO1_GCF_000006765.1/Pseudomonas_aeruginosa_PAO1_GCF_000006765.1.fna mutation-types.tab
    $ Rscript allele_count_filtering.R tmf/align output.tab
    
Other stuff.
