
sample=$1
mask=$2
fastq1=$3
fastq2=$4


out="${sample}.fastq.gz"
out_csv="${sample}.location.csv"

LD_LIBRARY_PATH=/ldfssz1/ST_BIGDATA/USER/liuxing2/lib/boost_1_73_0/lib/:/ldfssz1/ST_BIGDATA/USER/liuxing2/lib/hdf5-1.10.7/lib:/ldfssz1/ST_BIGDATA/USER/liuxing2/software/gcc-8.2.0/lib64
/ldfssz1/ST_BIGDATA/USER/liuxing2/shareData/software/ST_BarcodeMap --in $mask\
        --in1 $fastq1\
        --in2 $fastq2\
        --out $out\
        --mismatch 1\
        --thread 2

zcat $out|grep "|||" | awk -v FS="|" '{print $4}'|cut -b 6- | sort | uniq >> $out_csv

