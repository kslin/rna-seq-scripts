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

# ex: bash NEXTflex_pipeline.sh /lab/solexa_bartel/klin/RNA_Seq/Nam_HEK_transfections/mini_fastq/mir124_rep1_mini/mir_124_rep1_mini.fastq /nfs/genomes/human_hg38_dec13_no_random/STAR/GRCh38.81.canonical_overhang_100/ /nfs/genomes/human_hg38_dec13_no_random/gtf/Homo_sapiens.GRCh38.81.canonical.gtf /lab/solexa_bartel/klin/RNA_Seq/Nam_HEK_transfections/mini_fastq/mir124_rep1_mini
# bash NEXTflex_pipeline.sh /lab/solexa_bartel/klin/RNA_Seq/Nam_HEK_transfections/mir_124_rep1/mir_124_rep1.fastq /nfs/genomes/human_hg38_dec13_no_random/STAR/GRCh38.81.canonical_overhang_100/ /nfs/genomes/human_hg38_dec13_no_random/gtf/Homo_sapiens.GRCh38.81.canonical.gtf /lab/solexa_bartel/klin/RNA_Seq/Nam_HEK_transfections/mir_124_rep1

fastq_file="$1"
star_index="$2"
annot_file="$3"
# intron_file="$4"
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
if [ ! -f "$annot_file" ] ;
then
   echo -e "Cannot find exon annotation file $annot_file";
   exit
fi
# if [ ! -f "$intron_file" ] ;
# then
#    echo -e "Cannot find intron annotation file $intron_file";
#    exit
# fi
if [ ! -d "$out_dir" ] ;
then
   echo -e "Cannot find directory $out_dir";
   exit
fi

## reverse complement sequences - required for Bioo NEXTflex library prep
## -Q 33: base 33 quality scores
# fastx_reverse_complement -i "$fastq_file" -o "$out_dir/${filename}_revcomp.fastq" -Q "33" &&
fastx_reverse_complement -i "$fastq_file" -o "$out_dir/${filename}_revcomp.fastq" &&

## align to STAR
## flags: borrowed from ENCODE settings
STAR --genomeDir "$star_index" --runThreadN "30" --outFilterMultimapNmax "1" --outFilterMismatchNoverLmax "0.04" --outFilterIntronMotifs "RemoveNoncanonicalUnannotated" --outSJfilterReads "Unique" --outReadsUnmapped "Fastx" --readFilesIn "$out_dir/${filename}_revcomp.fastq" --outFileNamePrefix "$out_dir/${filename}_revcomp_star_" > "$out_dir/${filename}_revcomp_star_logFile.txt" &&

## convert SAM to BAM, then sort BAM file
samtools view -bS "$out_dir/${filename}_revcomp_star_Aligned.out.sam" | samtools sort -T "$out_dir/${filename}_TEMP2" -@ "24" -o "$out_dir/${filename}_revcomp_star_sorted.bam" &&

## annotate
htseq-count -f "bam" -s "yes" -t "transcript" -i "ID" -o "$out_dir/${filename}_revcomp_star_sorted_counts.sam" "$out_dir/${filename}_revcomp_star_sorted.bam" "$annot_file" &&

## compile counts
python "map_reads/count_reads.py" "$out_dir/${filename}_revcomp_star_sorted_counts.sam" "$out_dir/${filename}_revcomp_star_sorted_counts_compiled.txt" &&

## intersect with annotation files
## flags: write output as bed, use features from both files, stranded, first input bam, second input BED
# intersectBed -bed -wo -s -abam "$out_dir/${filename}_revcomp_star_sorted.bam" -b "$annot_file" > "$out_dir/${filename}_revcomp_star_sorted_annotated.bed" &&

# compile counts
# cut -f "16" "$out_dir/${filename}_revcomp_star_sorted_annotated.bed" | "sort" | "uniq" -c > "$out_dir/${filename}_revcomp_star_sorted_annotated_counts.txt"
# cut -f "16" "$out_dir/${filename}_revcomp_star_sorted_introns.bed" | "sort" | "uniq" -c > "$out_dir/${filename}_revcomp_star_sorted_intron_counts.txt" &&

## remove extraneous files
# rm -r "$out_dir/${filename}_revcomp_star__STARtmp"
# rm "$out_dir/${filename}_revcomp.fastq"
rm "$out_dir"/*.sam
# rm "$out_dir/${filename}_revcomp_star_logFile.txt"
# rm "$out_dir/${filename}_revcomp_star_SJ.out.tab"
# rm "$out_dir"/*.bed
# rm "$out_dir"/*.out
