# !/bin/sh
#
# Pipeline for analyzing RNA-Seq data from Bioo NEXTflex library prep
#

# check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo -e "-h print this screen and exit"
    printf "Usage: `basename $0` <bam file> <rRNA annotation file> <output directory>";
    exit
fi

# ex: bash NEXTflex_pipeline.sh /lab/solexa_bartel/klin/RNA_Seq/Nam_HEK_transfections/mini_fastq/mir124_rep1_mini/mir_124_rep1_mini.fastq /nfs/genomes/human_hg38_dec13_no_random/STAR/GRCh38.81.canonical_overhang_100/ /nfs/genomes/human_hg38_dec13_no_random/gtf/Homo_sapiens.GRCh38.81.canonical.gtf /lab/solexa_bartel/klin/RNA_Seq/Nam_HEK_transfections/mini_fastq/mir124_rep1_mini
# bash NEXTflex_pipeline.sh /lab/solexa_bartel/klin/RNA_Seq/Nam_HEK_transfections/mir_124_rep1/mir_124_rep1.fastq /nfs/genomes/human_hg38_dec13_no_random/STAR/GRCh38.81.canonical_overhang_100/ /nfs/genomes/human_hg38_dec13_no_random/gtf/Homo_sapiens.GRCh38.81.canonical.gtf /lab/solexa_bartel/klin/RNA_Seq/Nam_HEK_transfections/mir_124_rep1

bam_file="$1"
rRNA_file="$2"
out_dir="$3"

filename=$(basename $1 "$fullfile")
filename="${filename%.*}"

# check existence of input files and directories
if [ ! -f "$bam_file" ] ;
then
   echo -e "Cannot find bam file $bam_file";
   exit
fi
if [ ! -f "$rRNA_file" ] ;
then
   echo -e "Cannot find rRNA annotation file $rRNA_file";
   exit
fi
if [ ! -d "$out_dir" ] ;
then
   echo -e "Cannot find directory $out_dir";
   exit
fi

## intersect with annotation files
## flags: write output as bed, use features from both files, stranded, first input bam, second input BED
intersectBed -bed -wo -s -abam "$bam_file" -b "$rRNA_file" > "$out_dir/${filename}_sorted_rRNAs.bed" &&

# compile counts
cut -f "16" "$out_dir/${filename}_sorted_rRNAs.bed" | "sort" | "uniq" -c > "$out_dir/${filename}_sorted_rRNA_counts.txt"

## remove extraneous files
# rm -r "$out_dir/${filename}_revcomp_star__STARtmp"
# rm "$out_dir/${filename}_revcomp.fastq"
# rm "$out_dir"/*.sam
# rm "$out_dir/${filename}_revcomp_star_logFile.txt"
# rm "$out_dir/${filename}_revcomp_star_SJ.out.tab"
# rm "$out_dir"/*.bed
# rm "$out_dir"/*.out
