# Software tools for <i>Genome-wide protein–DNA interaction site mapping using a double strand DNA-specific cytosine deaminase</i>
![logo](/title.png)
## Guide for Use
These tools can be used to analyze data produced as described in the pre-print:

  https://www.biorxiv.org/content/10.1101/2021.08.01.454665v1

### Prerequisites
This software has been tested on the Linux operating system. It may be possible to adapt it for other operating systems. The following dependencies were used for the paper, though other versions may work:
* Python 3.8.5
* Python libraries:
  * Biopython 1.78
  * Pandas 1.3.0
  * colorama 0.4.3
  * PySAM 0.16.0.1
  * gffutils 0.10.1
* R 3.6.3
* R libraries:
  * evobiR 1.1
* HTStream 1.3.0 
* Minimap2 2.17-r974-dirty
* SAMTools 1.10
* BCFTools 1.10
* Picard Tools 2.18.25

The prerequisites can be installed on a vanilla Ubuntu 20.04 machine using the installation script:

    $ git clone https://github.com/marade/3DSeqTools.git
    $ cd 3DSeqTools
    $ bash ubuntu-install.sh

Depending on your environment, installing the prerequisites with a package sytem like Bioconda or using Python virtual environments may be advisable.

### Do Sequencing and Generate Fastq Files
We assume you have generated your sequencing data in roughly the manner described in the paper, using Illumina paired-end sequencing. We provide some example files for testing below.
### Prepare Fastq Files
This pipeline assumes your paired-end Fastq files are named like so:

    sampleX_1.fastq.gz sampleX_2.fastq.gz
    sampleY_1.fastq.gz sampleY_2.fastq.gz

### Run the Pipeline
The following pipeline works on a vanilla Ubuntu 20.04 installation with the prerequisites installed. We were able to process these data in the cloud with an AWS EC2 C5.2xlarge instance (8 compute cores; 16GB RAM).

    # If you haven't cloned the repo already...
    $ git clone https://github.com/marade/3DSeqTools.git
    $ cd 3DSeqTools
    #
    $ mkdir fastq && cd fastq
    # Note you are downloading nearly 4GB of fastq data here:
    $ for x in {1..4}; do for y in {1..2}; do wget https://3d-seq-01.s3.us-west-2.amazonaws.com/fastq/gcsR-dddA-Para-dddI-delta-ung-0Ara-P3-rep${x}_${y}.fastq.gz ; done ; done
    $ cd ../
    $ python3 AlignReads -n 8 fastq NCBI/Pseudomonas_aeruginosa_PAO1_GCF_000006765.1 tmf
    $ python3 Process3DSeq -n 3 tmf/align NCBI/Pseudomonas_aeruginosa_PAO1_GCF_000006765.1/Pseudomonas_aeruginosa_PAO1_GCF_000006765.1.fna mutation-types.tab
    $ Rscript allele_count_filtering.R tmf/align output.tab
    $ awk '$3 != "0"' output.tab
    
The output shows a peak corresponding to the GcsR binding site:

        pos     maaf    mov_avg_75bp
  334409  2745370 0.0290809928898512      0.00193873285932341
  334424  2745437 0.0253027635979581      0.00253027635979581
  334555  2746633 0.0326857649188563      0.00662365117290723
  334556  2746636 0.033550746810216       0.0060215010662793
  334567  2746697 0.0372025356372053      0.00372025356372053
  334573  2746763 0.0493712494744486      0.0146799269917299
  334574  2746770 0.0565141345722283      0.0157231091251582
  334576  2746780 0.0409138858706218      0.0172954200376741
  334580  2746805 0.0261549304594418      0.0159523385938013
  334583  2746829 0.0518927736295225      0.0122110012825046
  334588  2746865 0.0562733100185859      0.0191399425312752
  334590  2746879 0.0832333416646441      0.0155007390759144
  334594  2746919 0.0247753203440407      0.00715313293050461
  334597  2746946 0.0252966101694915      0.0174156425821467
  334599  2746970 0.0299255477207476      0.0212261626242594
  334601  2746982 0.0593276624228934      0.0238794329522918
  334602  2746999 0.0247107555304515      0.0385466622605808
  334603  2747005 0.0517748877747507      0.0346919960345227
  334606  2747024 0.0909572221491726      0.0378420292181693
  334607  2747027 0.0902238847472113      0.0378420292181693
  334609  2747048 0.0450694835437686      0.0275057076486432
  334610  2747070 0.0213007783976358      0.0113953496494308
  334616  2747105 0.0361878849034729      0.024964450963095
  334617  2747109 0.127178383375039       0.0253715256925905
  334618  2747114 0.0649774629548026      0.0253715256925905
  334622  2747205 0.162003790782861       0.0270006317971434
  334626  2747245 0.138547006745569       0.0226411780409983
  334630  2747275 0.0425824175824176      0.0254992302772907
  334631  2747290 0.0484707392544127      0.0180916016332232
  334633  2747301 0.0253921391905075      0.0221119575517173
  334636  2747313 0.0274651102348634      0.0194229805045071
  334638  2747326 0.0550972117032544      0.0174806824540564
  334639  2747340 0.0183816241575261      0.0144205637279491
  334641  2747381 0.0420410024757851      0.018476190445428
  334642  2747383 0.0688161401967828      0.018476190445428
  334648  2747431 0.0338919781322152      0.00962179760514335
  334653  2747464 0.0623259979192183      0.00962179760514335
  334657  2747506 0.0503100640655996      0.0176961971013872
  334658  2747518 0.0735633156441108      0.0247746759419421
  334661  2747569 0.0187859794638032      0.0105536629295164
  334662  2747575 0.0445359981132955      0.0126643955154197
  334668  2747636 0.037156732875911       0.0037156732875911
  334675  2747680 0.0385734436509683      0.00624079596273641
  334680  2747711 0.0300753119391322      0.0114414592650167
  334688  2747837 0.0285755310827374      0.00680606088393682
  334696  2747873 0.0599032604084412      0.0111463514050466
  334698  2747879 0.0279746410137343      0.0121410940491452
  334700  2747895 0.0507418385707859      0.0149428849835634
  334703  2747911 0.0556377647933625      0.0196078438833174
  334705  2747935 0.0617241944552916      0.0130402176942949
  334739  2748231 0.0423671628246994      0.0058758771398685
  334741  2748251 0.0222674857138541      0.00646346485385535
  334774  2748549 0.0365225415989314      0.00365225415989314
  334779  2748595 0.0247510495596969      0.00247510495596969

