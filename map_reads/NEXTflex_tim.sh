#!/bin/bash

#  TAILseq_pipeline_v1_TJE.sh
#
#
#  Created by Timothy Eisen on 6/19/15.
#  Updated by Timothy Eisen on 10/03/15.
#  Updated by Timothy Eisen on 1/30/16 to create TAIL-seq modifications
#  Updated by Timothy Eisen on 2/22/16
#  Version 3 makes the following modifications:
#   Performs two parallel mappings on the adapters_removed fastq1 file
#   one mapping uses star with a transcriptome that does not contain the standard indices (mouse mm10)
#   the second mapping uses bowtie with an index of the STD13 standards only (chrS_...) and maps these two a sam file
#   awk, grep, cat to combine these this second sam file with the intersect bed output converted sam file (gene_assignments.sam)
#   the output of adapters removed isn't gzipped until after it is used in the alignments (gzip...)

#Requires that at least one argument is present
if (! getopts ":1:s:b:d:" opt); then
    echo -e "-h print this screen and exit"
    printf "Usage: `basename $0` \nRequired arguments: \n-1 Fastq1 \n-s STAR index \n-b Bed annotation file \n-d Directory\n";
    exit $E_OPTERROR;
fi

#requires that the correct number of arguments is present
if [ "$#" -ne 8 ]; then
    echo -e "\nIllegal number of arguments (Error #2)"
    echo -e "-h print this screen and exit"
    printf "Usage: `basename $0` \nRequired arguments: \n-1 Fastq1 \n-s STAR index \n-b Bed annotation file\n-d Directory\n";
    exit
fi

##parses argument flags
while getopts ":1:s:b:i:d:h" opt; do
    case $opt in
        1)
            echo "$OPTARG"
            fastq1=$OPTARG
            ;;
        s)
            star_index=$OPTARG
            ;;
        b)
            bed_file=$OPTARG
            ;;
   
        d)
            echo "$OPTARG"
            directory=$OPTARG
            ;;
        h)
            echo -e "-h print this screen and exit"
            exit $E_OPTERROR
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

echo "$directory used as the working directory"

##filename extraction
filename=$(basename ${fastq1%%.*})

##creates log file with input commands
# echo "Log file $(date +%y_%m_%d)" > $directory/log_file.txt;
# echo "FASTQ1 file: $fastq1" >>$directory/log_file.txt;
# echo "STAR index: $star_index" >>$directory/log_file.txt;
# echo "BED file: $bed_file" >>$directory/log_file.txt;
# echo "Output file directory: $directory" >>$directory/log_file.txt;


##unzip data
# bsub -q bartel -J "unzip_fastq1" "tar xzfO $fastq1>$directory/${filename}_fastq1.txt"

##rev comp the fastq file
# bsub -q bartel -J "revcomp" "fastx_reverse_complement -i $fastq1 -o $directory/${filename}_fastq_revcomp.txt"

##sort bed file
# bsub -q bartel -J "sort" "sort -k1,1V -k2,2n $bed_file > $directory/${filename}_bed_sorted.txt"

##STAR commands, arguments
# bsub -q bartel -w "ended(revcomp)" -J "STAR" "STAR --genomeDir $star_index --runThreadN 30 --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSJfilterReads Unique --outReadsUnmapped Fastx --readFilesIn $directory/${filename}_fastq_revcomp.txt --outFileNamePrefix $directory/${filename}_map1_ > $directory/${filename}_stdOut_logFile.txt"

##collapse intron spanning reads and removes reads with soft clipping. 
# bsub -q bartel -w "ended(STAR)" -J "remove_clipped" "bash /lab/solexa_bartel/teisen/RNAseq/Scripts/intron_spanning_read_collapse_TAILseq_edit.sh $directory/${filename}_map1_Aligned.out.sam $directory/${filename}_map1_Aligned.out_intron_spanning_read_collapse.sam"

##sam to bam conversion-This line was edited by TJE on 20160916 to add a proper sorting prefix
bsub -q bartel -J "sam_to_bam" "samtools view -bS $directory/${filename}_map1_Aligned.out.sam | samtools sort -T $directory/${filename}_temp -@ 8 -O bam > $directory/${filename}_map1.bam"

##bam to interesect bed output. This script takes a significant portion of the pipeline in terms of runtime.
bsub -q bartel -w "ended(sam_to_bam)" -J "bam_to_intersectBed" "intersectBed -bed -abam $directory/${filename}_map1.bam -b $directory/${filename}_bed_sorted.txt -wo -s > $directory/${filename}_wo_output.txt"

bsub -q bartel -w "ended(bam_to_intersectBed)" -J "compile" "cut -f 16 $directory/${filename}_wo_output.txt | sort | uniq -c | sort -rn > $directory/${filename}_gene_assignment_compiled.txt"


###zip everything and remove intermediate files
#bsub -q bartel w "ended(plots)" -J "gzip_txt" "gzip $directory/*.txt"
#bsub -q bartel w "ended(plots)" -J "gzip_sam" "gzip $directory/*.sam"




