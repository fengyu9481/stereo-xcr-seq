version 1.0

workflow cyclone_immune_repertoire_v1{
    input{
        String Outdir
        String SampleID
        String Cyclone_fastq
        String species
        String BarcodeWhiteList
        String BarcodeIndex
        String chain
        Int OminimalQuality
        Int ANN
        Int Kmer
        Int Chunk_size
        Int BarcodeMapCPU
        Boolean DryRun
    }

    String Root="/jdfssz2/ST_SUPERCELLS/P22Z10200N0739/liuyi/scirpt/cyclone_immune_repertoire"

    call makedir{
        input:
        outdir=Outdir
    }
    if (!DryRun){
        call splitFix{
            input:
            root=Root,
            cyclone_fastq=Cyclone_fastq,
            outdir=Outdir,
            sampleId=SampleID
        }
        call MixcrAlign{
            input:
            root=Root,
            inputfastq=splitFix.SplitR2Fastq,
            outdir=Outdir,
            chain=chain,
            sampleId=SampleID,
            species=species,
            score=OminimalQuality
        }
        call SubsetAlign{
            input:
            root=Root,
            sampleId=SampleID,
            mixcr_align=MixcrAlign.mixcrAlign,
            fastq_need_subset=splitFix.SplitR2Fastq,
            clone_info=MixcrAlign.CloneInfo,
            outdir=Outdir,
            barcode_file=splitFix.SplitBarcode
        }
        call BarocdeMap{
            input:
            root=Root,
            sampleId=SampleID,
            alignbarcode=SubsetAlign.AlignBarcode,
            white_list=BarcodeWhiteList,
            chip_index=BarcodeIndex,
            outdir=Outdir,
            ann=ANN,
            k=Kmer,
            chunk_size=Chunk_size,
            cpu=BarcodeMapCPU
        }
        call ModifyCDR3{
            input:
            root=Root,
            sampleId=SampleID,
            fastq_file=splitFix.SplitR2Fastq,
            match_file=BarocdeMap.MappedR1,
            cdr3_file=SubsetAlign.ReadsCDR3,
            align_file=MixcrAlign.mixcrAlign,
            outdir=Outdir,
            chain=chain
        }
    }
    if (DryRun){
        call wdl_dryrun{
            input:
            root=Root,
            cyclone_fastq=Cyclone_fastq,
            outdir=Outdir,
            sampleId=SampleID,
            species=species,
            chain=chain,
            score=OminimalQuality,
            white_list=BarcodeWhiteList,
            chip_index=BarcodeIndex,
            ann=ANN,
            k=Kmer,
            chunk_size=Chunk_size,
            cpu=BarcodeMapCPU
        }
    }
}

## task.1 Create directory ##
task makedir{
    input{
        String outdir
    }
    command{
        mkdir -p ${outdir}
        mkdir -p ${outdir}/01.splitFix
        mkdir -p ${outdir}/02.mixcr
        mkdir -p ${outdir}/02.mixcr/report
        mkdir -p ${outdir}/02.mixcr/clone_fastq
        mkdir -p ${outdir}/03.SubsetAlign
        mkdir -p ${outdir}/04.BarcodeMap
    }
}

## task.2 Split Fixed Sequence ##
task splitFix{
    input{
        String root
        String cyclone_fastq
        String outdir
        String sampleId
    }
    command{
        ${root}/bin/python ${root}/scripts/ReadSplit.py ${cyclone_fastq} ${outdir}/01.splitFix ${sampleId} 10
    }

    output{
        File SplitBarcode="${outdir}/01.splitFix/${sampleId}.barcode.csv.gz"
        File SplitR2Fastq="${outdir}/01.splitFix/${sampleId}.r2.fastq.gz"
        File FixMatch="${outdir}/01.splitFix/fix_search.pdf"
        File R1ToFixLenth="${outdir}/01.splitFix/r1_to_fix.lenth.pdf"
        File R1ToR2Lenth="${outdir}/01.splitFix/r1_to_r2.lenth.pdf"
        File BarcodeLenth="${outdir}/01.splitFix/barcode.lenth.pdf"
        File R2Lenth="${outdir}/01.splitFix/r2.lenth.pdf"
        File PattenPie="${outdir}/01.splitFix/patten.pie.pdf"
        File ValidRead="${outdir}/01.splitFix/valid_read.info"
        File SplitReport="${outdir}/01.splitFix/step1.Split.report"
    }
}

## task.3 Mixcr Align R2 ##
task MixcrAlign{
    input{
        String root
        String inputfastq
        String outdir
        String sampleId
        String species
        String chain
        Int score
    }
    command{
        ${root}/bin/mixcr align -Xmx100g -f --report ${outdir}/02.mixcr/report/${sampleId}.align.report.txt --json-report ${outdir}/02.mixcr/report/${sampleId}.align.report.json --threads 10 --use-local-temp --preset generic-ont --species ${species} -M saveOriginalReads=true ${inputfastq} ${outdir}/02.mixcr/${sampleId}.vdjca
        ${root}/bin/mixcr assemble -f -OassemblingFeatures='CDR3' -OseparateByJ=true -OseparateByC=true -OseparateByV=true -OminimalQuality=${score} -a --report ${outdir}/02.mixcr/report/${sampleId}.assemble.report.txt --json-report ${outdir}/02.mixcr/report/${sampleId}.assemble.report.json ${outdir}/02.mixcr/${sampleId}.vdjca ${outdir}/02.mixcr/${sampleId}.clna 
        ${root}/bin/mixcr exportClones --dont-split-files --chains ${chain} ${outdir}/02.mixcr/${sampleId}.clna ${outdir}/02.mixcr/${sampleId}.clone.tsv
        ${root}/bin/mixcr exportAlignments --drop-default-fields -readIds -descrsR1 -cloneId -vHitsWithScore -dHitsWithScore -jHitsWithScore -cHitsWithScore -topChains -allNFeatures -f ${outdir}/02.mixcr/${sampleId}.clna ${outdir}/02.mixcr/${sampleId}.align.tsv
    }

    output{
        File vdjca="${outdir}/02.mixcr/${sampleId}.vdjca"
        File assemble_cla="${outdir}/02.mixcr/${sampleId}.clna"
        File mixcrAlign="${outdir}/02.mixcr/${sampleId}.align.tsv"
        File CloneInfo="${outdir}/02.mixcr/${sampleId}.clone.tsv"
    }
}

## task.4 Subset MIXCR result From R2 and Barcode ##
task SubsetAlign{
    input{
        String root
        String sampleId
        String mixcr_align
        String fastq_need_subset
        String clone_info
        String outdir
        String barcode_file
    }
    command{
        ${root}/bin/python ${root}/scripts/SubsetAlign.py ${sampleId} ${mixcr_align} ${fastq_need_subset} ${clone_info} ${outdir}/03.SubsetAlign ${barcode_file}
    }

    output{
        File AlignR2="${outdir}/03.SubsetAlign/${sampleId}.align.r2.fastq.gz"
        File AlignBarcode="${outdir}/03.SubsetAlign/${sampleId}.align.barcode.csv.gz"
        File ReadsCDR3="${outdir}/03.SubsetAlign/${sampleId}.ReadsCDR3.csv"
        File AlignReport="${outdir}/03.SubsetAlign/step2.Mixcr.report"
    }
}


## task.5 Mapping Barcode ##
task BarocdeMap{
    input{
        String root
        String sampleId
        String alignbarcode
        String white_list
        String chip_index
        String outdir
        Int ann=10000
        Int k=5
        Int chunk_size=5000
        Int cpu=2
    }
    command{
        ${root}/bin/python ${root}/scripts/BarcodeMap.py ${sampleId} ${alignbarcode} ${white_list} ${chip_index} ${outdir}/04.BarcodeMap ${cpu} ${ann} ${k} ${chunk_size}
    }

    output{
        File MappedR1="${outdir}/04.BarcodeMap/${sampleId}.match.csv"
        File MapReport="${outdir}/04.BarcodeMap/step3.BarcodeMap.report"
    }
}

## task.6 ModifyCDR3 ##
task ModifyCDR3{
    input{
        String root
        String sampleId
        String fastq_file
        String match_file
        String cdr3_file
        String align_file
        String outdir
	String chain
    }
    command{
        ${root}/bin/python ${root}/scripts/ModifyCDR3.py ${match_file} ${cdr3_file} ${outdir}/04.BarcodeMap/${sampleId}
    }
}


task wdl_dryrun{
    input{
        String root
        String cyclone_fastq
        String outdir
        String sampleId
        String species
        String chain
        Int score
        String white_list
        String chip_index
        Int ann=10000
        Int k=5
        Int chunk_size=5000
        Int cpu=2
    }
    command{
        echo "## Split" >> "${outdir}/${sampleId}.dryrun.sh"
        echo "${root}/bin/python ${root}/scripts/ReadSplit.py ${cyclone_fastq} ${outdir}/01.splitFix ${sampleId} 10" > "${outdir}/${sampleId}.dryrun.sh"
        echo "## MIXCR" >> "${outdir}/${sampleId}.dryrun.sh"
        echo "${root}/bin/mixcr align -Xmx100g -f --report ${outdir}/02.mixcr/report/${sampleId}.align.report.txt --json-report ${outdir}/02.mixcr/report/${sampleId}.align.report.json --threads 10 --use-local-temp --preset generic-ont --species ${species} -M saveOriginalReads=true ${outdir}/01.splitFix/${sampleId}.r2.fastq.gz ${outdir}/02.mixcr/${sampleId}.vdjca" >> "${outdir}/${sampleId}.dryrun.sh"
        echo "${root}/bin/mixcr assemble -f -OassemblingFeatures='CDR3' -OseparateByJ=true -OseparateByC=true -OseparateByV=true -OminimalQuality=${score} -a --report ${outdir}/02.mixcr/report/${sampleId}.assemble.report.txt --json-report ${outdir}/02.mixcr/report/${sampleId}.assemble.report.json ${outdir}/02.mixcr/${sampleId}.vdjca ${outdir}/02.mixcr/${sampleId}.clna">> "${outdir}/${sampleId}.dryrun.sh" 
        echo "${root}/bin/mixcr exportClones --dont-split-files --chains ${chain} ${outdir}/02.mixcr/${sampleId}.clna ${outdir}/02.mixcr/${sampleId}.clone.tsv" >> "${outdir}/${sampleId}.dryrun.sh"
        echo "${root}/bin/mixcr exportAlignments --drop-default-fields -readIds -descrsR1 -cloneId -vHitsWithScore -dHitsWithScore -jHitsWithScore -cHitsWithScore -topChains -allNFeatures -f ${outdir}/02.mixcr/${sampleId}.clna ${outdir}/02.mixcr/${sampleId}.align.tsv" >> "${outdir}/${sampleId}.dryrun.sh"
        echo "## SubsetAlign" >> "${outdir}/${sampleId}.dryrun.sh"
        echo "${root}/bin/python ${root}/scripts/SubsetAlign.py ${sampleId} ${outdir}/02.mixcr/${sampleId}.align.tsv ${outdir}/01.splitFix/${sampleId}.r2.fastq.gz ${outdir}/02.mixcr/${sampleId}.clone.tsv ${outdir}/03.SubsetAlign ${outdir}/01.splitFix/${sampleId}.barcode.csv.gz" >> "${outdir}/${sampleId}.dryrun.sh"
        echo "## Barcode Mapping" >> "${outdir}/${sampleId}.dryrun.sh"
        echo "${root}/bin/python ${root}/scripts/BarcodeMap.py ${sampleId} ${outdir}/03.SubsetAlign/${sampleId}.align.barcode.csv.gz ${white_list} ${chip_index} ${outdir}/04.BarcodeMap ${cpu} ${ann} ${k} ${chunk_size}" >> "${outdir}/${sampleId}.dryrun.sh"
        echo "## ModifyCDR3" >> "${outdir}/${sampleId}.dryrun.sh"
        echo "${root}/bin/python ${root}/scripts/ModifyCDR3.py ${outdir}/04.BarcodeMap/${sampleId}.match.csv ${outdir}/03.SubsetAlign/${sampleId}.ReadsCDR3.csv ${outdir}/04.BarcodeMap/${sampleId} " >> "${outdir}/${sampleId}.dryrun.sh"
    }
}
