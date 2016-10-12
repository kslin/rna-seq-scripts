# !/bin/sh
#
# Pipeline for analyzing RNA-Seq data from Bioo NEXTflex library prep
#

# check for the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo -e "-h print this screen and exit"
    printf "Usage: `basename $0` <fastq file> <STAR index> <annotation file> <output directory>";
    exit
fi

fastq_file="$1"
star_index="$2"
gtf_file="$3"
out_dir="$4"

filename=$(basename $1 "$fullfile")
filename="${filename%.*}"

# check existence of input files and directories
if [ ! -f "$fastq_file" ] ;
then
   echo -e "Cannot find fastq file $fastq_file";
   exit
fi
if [ ! -d "$star_index" ] ;
then
   echo -e "Cannot find directory $star_index";
   exit
fi
if [ ! -f "$gtf_file" ] ;
then
   echo -e "Cannot find annotation file $gtf_file";
   exit
fi
if [ ! -d "$out_dir" ] ;
then
   echo -e "Cannot find directory $out_dir";
   exit
fi


## reverse complement sequences - required for Bioo NEXTflex library prep
## -Q 33: base 33 quality scores
# bsub -q bartel -J "revcomp" "fastx_reverse_complement -i $fastq_file -o $out_dir/${filename}_revcomp.fastq -Q 33"

## align to STAR
## flags: borrowed from ENCODE settings
# bsub -q bartel -w "ended(revcomp)" -J "STAR" "STAR --genomeDir $star_index --runThreadN 30 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSJfilterReads Unique --outReadsUnmapped Fastx --readFilesIn $out_dir/${filename}_revcomp.fastq --outFileNamePrefix $out_dir/${filename}_star_ > $out_dir/${filename}_stdOut_logFile.txt"
# bsub -q bartel -J "STAR" "STAR --genomeDir $star_index --runThreadN 30 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSJfilterReads Unique --outReadsUnmapped Fastx --readFilesIn $out_dir/${filename}_revcomp.fastq --outFileNamePrefix $out_dir/${filename}_star_ > $out_dir/${filename}_stdOut_logFile.txt"

## sort SAM file
## flags: SAM outformat, outfile, infile sorted by name, temp filename, number of threads
# bsub -q bartel -w "ended(STAR)" -J "sort" "samtools sort -O sam -o $out_dir/${filename}_revcomp_star_sorted.sam -n -T $out_dir/${filename}_TEMP -@ 24"
# bsub -q bartel -J "sort" "samtools sort -O sam -o $out_dir/${filename}_revcomp_star_sorted.sam -n -T $out_dir/${filename}_TEMP -@ 24 $out_dir/${filename}_star_Aligned.out.sam"

## use htseq-count to count reads, only use uniquely mapped reads
## flags: SAM infile format, infile sorted by name, annotate exons
# bsub -q bartel -w "ended(STAR)" -J "htseq-count" "htseq-count --format sam --order name --type exon $out_dir/${filename}_revcomp_star_sorted.sam $gtf_file"
bsub -q bartel -J "htseq-count" "htseq-count --format sam --order name --type exon $out_dir/${filename}_revcomp_star_sorted.sam $gtf_file"

## remove extraneous files
# TODO
