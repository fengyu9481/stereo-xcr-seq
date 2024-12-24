#!/bin/bash
set -e
sampleID="$1"
chipID="$2"
sequence="GTCTTAGGAAGACAA"

# Step1: Fixed Sequence Split #
echo "====  Step1 start: Fixed Sequence Split ===="

python seq_search.py ${sampleID}_read_2.fq.gz ${sampleID}_read_2.fixed.fq.gz ${sequence}

echo "====  Step1 done  ===="

# Step2: BarcodeMap with h5 file #
echo "====  Step2 start: BarcodeMap with h5 file ===="

ST_BarcodeMap --in ${chipID}.barcodeToPos.h5 --in1 ${sampleID}_read_2.fixed.fq.gz --in2 ${sampleID}_read_2.fixed.fq.gz --out ${sampleID}_read_2.fixed.barcode.fq.gz --mismatch 1 --thread 2 --umiStart 25

echo "====  Step2 done  ===="





