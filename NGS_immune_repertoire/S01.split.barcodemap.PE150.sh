#!/bin/bash
set -e
sampleID="$1"
rawdataDir="$2"
sequence="GTCTTAGGAAGACAA"
chipID="$3"
maskDir="$4"

# software
pythonDir="/jdfssz1/ST_SUPERCELLS/P22Z10200N0739/liuyi/00.software/conda/bin"
scriptDir="/jdfssz2/ST_SUPERCELLS/P22Z10200N0739/yanyixin/Cyclone/analysis/PE150/script/split"
ST_BarcodeMapDir="/ldfssz1/ST_BIGDATA/USER/liuxing2/shareData/software"

# Step1: Fixed Sequence Split #
echo "====  Step1 start: Fixed Sequence Split ===="
date

${pythonDir}/python ${scriptDir}/seq_search.py ${rawdataDir}/${sampleID}_read_2.fq.gz ${sampleID}_read_2.fixed.fq.gz ${sequence}

date
echo "====  Step1 done  ===="

# Step2: BarcodeMap with h5 file #
echo "====  Step2 start: BarcodeMap with h5 file ===="
date

LD_LIBRARY_PATH=/ldfssz1/ST_BIGDATA/USER/liuxing2/lib/boost_1_73_0/lib/:/ldfssz1/ST_BIGDATA/USER/liuxing2/lib/hdf5-1.10.7/lib:/ldfssz1/ST_BIGDATA/USER/liuxing2/software/gcc-8.2.0/lib64

${ST_BarcodeMapDir}/ST_BarcodeMap --in ${maskDir}/${chipID}.barcodeToPos.h5 --in1 ${sampleID}_read_2.fixed.fq.gz --in2 ${sampleID}_read_2.fixed.fq.gz --out ${sampleID}_read_2.fixed.barcode.fq.gz --mismatch 1 --thread 2 --umiStart 25

date
echo "====  Step2 done  ===="





