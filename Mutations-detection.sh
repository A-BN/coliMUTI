#################################################################################
							Quality-trim_Illumina-reads
#################################################################################

module load fastqc/0.11.9
module load multiqc/1.13
module load trim-galore/0.6.5

cd /shared/projects/tolerance_mutations_ecoli/breseq_all/data-for-breseq/illumina_reads/

# Create the files list_sample.txt, list_R1.txt and list_R2.txt in the directory
ls /shared/projects/tolerance_mutations_ecoli/breseq_all/data-for-breseq/illumina_reads/fastq/*.fastq* > list.samples.txt
find /shared/projects/tolerance_mutations_ecoli/breseq_all/data-for-breseq/illumina_reads/fastq/ -maxdepth 1 -name "*_R1.fastq*" > list_R1.txt
sed 's/_R1.fastq/_R2.fastq/g' list_R1.txt > list_R2.txt

## (1) Quality
# Run fastqc
mkdir fastqc_results

while read line
do
	fastqc -o fastqc_results $line
done < list.samples.txt

# Run multiqc
mkdir multiqc_results
cd multiqc_results/

multiqc -ip ../fastqc_results/

## (2) Trimming by pairs of reads
cd ..

mkdir trim_results

count=1
while read lineA
    do 
        lineB=`sed -n "$count"p list_R2.txt`
        count=`expr $count + 1`
        trim_galore -q 30 --illumina --paired -o trim_results --length 50 $lineA $lineB

done < list_R1.txt

echo '######################'
echo 'Job finished' $(date --iso-8601=seconds)

#################################################################################
									BreSeq
#################################################################################

module load breseq/0.35.0

cd /shared/projects/tolerance_mutations_ecoli/breseq_all/

mkdir output_breseq
cd output_breseq

#Create the files, list_val_R1.txt and list_val_R2.txt in the directory
find ../data-for-breseq/illumina_reads/trim_results/ -maxdepth 1 -name "*R1_val_1.fq" > list_val_R1.txt
sed 's/_R1_val_1.fq/_R2_val_2.fq/g' list_val_R1.txt > list_val_R2.txt

count=1
while read lineA
    do 
     	lineB=`sed -n "$count"p list_val_R2.txt`
        count=`expr $count + 1`
        
        d=`basename $lineA`
        d="${d%%_*}"
        
		breseq -j 8 -p -o $d -r ../data_reference/PAS_annot/PROKKA_12052022.gbk $lineA $lineB

done < list_val_R1.txt


echo '######################'
echo 'Job finished' $(date --iso-8601=seconds)

#################################################################################
									gdtool
#################################################################################

module load breseq/0.35.0

gd_folder="/shared/projects/tolerance_mutations_ecoli/breseq_all/output_breseq/gd_output"
cd ${gd_folder}
 
list_gd=$(ls -1 *.gd)
 
for curr_gd in ${list_gd}
do
	echo $curr_gd
	gdtools ANNOTATE \
		-o ${curr_gd}.tsv \
		-f TSV \
		-r ../PROKKA_12052022.gbk \
		${curr_gd}
done

echo '######################'
echo 'Job finished' $(date --iso-8601=seconds)