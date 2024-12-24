#!/bin/bash
sampleID="$1"
rawdataDir="$2"
species="hsa"
# species="mmu"

# software
mixcrDir="/home/fengyu/liuzhong/software/mixcr-4.6.0"
fastpDir="/home/fengyu/miniconda3/envs/bioconda/bin"

# mixcr analyze
${mixcrDir}/mixcr align \
    -p rna-seq \
    -Xmx50g \
    --threads 20 \
    --species ${species} \
    -f \
    -OsaveOriginalReads=true \
    -OallowPartialAlignments=true \
    -OvParameters.geneFeatureToAlign="VTranscriptWithout5UTRWithP" \
    --report ./${sampleID}.align.report.txt \
    --json-report ./${sampleID}.align.report.json \
    ${rawdataDir}/${sampleID}.fq.gz \
    ./${sampleID}.vdjca

${mixcrDir}/mixcr exportReads ${sampleID}.vdjca > ${sampleID}.mixcr_align.fastq

${mixcrDir}/mixcr assemble \
    -Xmx50g \
    -f \
    -OassemblingFeatures='CDR3' \
    -OseparateByJ=true\
    -OseparateByV=true \
    -a \
    --report ./${sampleID}.assemble.report.txt \
    --json-report ./${sampleID}.assemble.report.json \
    ./${sampleID}.vdjca ./${sampleID}.clna


${mixcrDir}/mixcr exportClones\
    -Xmx50g \
    -f\
    --dont-split-files \
    --prepend-columns \
    -topChains \
    -isotype primary\
    ./${sampleID}.clna ./${sampleID}.contigs.tsv


${mixcrDir}/mixcr exportAlignments \
    -Xmx50g \
    -readIds \
    -descrsR1\
    -cloneId\
    -f \
    ./${sampleID}.clna ./${sampleID}.align.tsv
	
	
date

