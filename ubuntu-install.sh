#!/bin/bash

sudo apt -y update

sudo apt -y install python3-pip pigz samtools bcftools minimap2 picard-tools

wget https://github.com/s4hts/HTStream/releases/download/v1.3.0-release/HTStream_v1.3.0-release.tar.gz && tar -xzvf HTStream_v1.3.0-release.tar.gz && sudo mv hts_* /usr/local/bin/ && rm HTStream_v1.3.0-release.tar.gz

sudo -H pip install pandas biopython colorama pysam
