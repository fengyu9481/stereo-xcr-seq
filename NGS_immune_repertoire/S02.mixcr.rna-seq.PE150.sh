#!/bin/bash
sampleID="$1"
species="hsa"
# species="mmu"

# mixcr analyze
mixcr align \
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
    ./${sampleID}.fq.gz \
    ./${sampleID}.vdjca

mixcr exportReads ${sampleID}.vdjca > ${sampleID}.mixcr_align.fastq

mixcr assemble \
    -Xmx50g \
    -f \
    -OassemblingFeatures='CDR3' \
    -OseparateByJ=true\
    -OseparateByV=true \
    -a \
    --report ./${sampleID}.assemble.report.txt \
    --json-report ./${sampleID}.assemble.report.json \
    ./${sampleID}.vdjca \
    ./${sampleID}.clna

mixcr exportClones \
    -Xmx50g \
    -f \
    --dont-split-files \
    --prepend-columns \
    -topChains \
    -isotype primary \
    ./${sampleID}.clna \
    ./${sampleID}.contigs.tsv


${mixcrDir}/mixcr exportAlignments \
    -Xmx50g \
    -readIds \
    -descrsR1 \
    -cloneId \
    -f \
    ./${sampleID}.clna \
    ./${sampleID}.align.tsv



