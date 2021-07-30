# Software tools for <i>Genome-wide protein–DNA interaction site mapping using a double strand DNA-specific cytosine deaminase</i>
![logo](/title.jpg)
## Guide for Use
### Prerequisites
This software has been tested on the Linux operating system. It may be possible to adapt it for other operating systems. The following software dependencies are required:
* Python 3.8.5
* Python libraries:
  * Biopython 1.78
  * Pandas 1.3.0
  * colorama 0.4.3
  * PySAM 0.16.0.1
* HTStream 1.3.0 
* Minimap2 2.17-r974-dirty
* SAMTools 1.10

### Do Sequencing and Generate Fastq Files
We assume you have generated your sequencing data in roughly the manner described in the paper, using Illumina paired-end sequencing. We provide some example files for testing below.
### Prepare Fastq Files
This pipeline assumes your paired-end Fastq files are named like so:

    sampleX_1.fastq.gz sampleX_2.fastq.gz
    sampleY_1.fastq.gz sampleY_2.fastq.gz

### Run the Pipeline
Code blah blah...

    $ git clone https://github.com/marade/3DSeqTools.git
    $ cd 3DSeqtools
    $ python3 AlignReads -n 15 fq NCBI/Pseudomonas_aeruginosa_PAO1_GCF_000006765.1 tmf
    $ Process3DSeq -n 3 tmf/align NCBI/Pseudomonas_aeruginosa_PAO1_GCF_000006765.1/Pseudomonas_aeruginosa_PAO1_GCF_000006765.1.fna mutation-types.tab
    
Other stuff.
