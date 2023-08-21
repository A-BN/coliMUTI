#############################################################################
						Filtlong-Nanopore-reads
#############################################################################

module load filtlong/0.2.0

cd /shared/projects/tolerance_mutations_ecoli/breseq_all/data_reference/nanopore_reads/

## unzip fastq files

gzip -d *.fastq.gz

## Concatenation des fastq

cat *.fastq > PAS_nanopore.fastq

filtlong --min_length 1000 --keep_percent 95 ./PAS_nanopore.fastq | gzip > ./PAS_long.fastq.gz

echo '######################'
echo 'Job finished' $(date --iso-8601=seconds)

#############################################################################
						QC_bbmap-Illumina-reads
#############################################################################
module load fastqc/0.11.9
module load multiqc/1.13
module load bbmap/39.00

cd /shared/projects/tolerance_mutations_ecoli/breseq_all/data_reference/illumina_reads/

## (1) Quality
# Run fastqc
mkdir fastqc_results

fastqc -o fastqc_results ./*.fastq.gz

# Run multiqc
mkdir multiqc_results
cd multiqc_results/

multiqc -ip ../fastqc_results/

## (2) Cleaning

cd /shared/projects/tolerance_mutations_ecoli/breseq_all/data_reference/illumina_reads/

mkdir fastq_clean
cd fastq_clean

ls /shared/projects/tolerance_mutations_ecoli/breseq_all/data_reference/illumina_reads/*.fastq* > list.samples.txt
find /shared/projects/tolerance_mutations_ecoli/breseq_all/data_reference/illumina_reads/ -maxdepth 1 -name "*R1_001.fastq*" > list_R1.txt
sed 's/_R1_001.fastq/_R2_001.fastq/g' list_R1.txt > list_R2.txt

count=1
while read lineA
    do 
        lineB=`sed -n "$count"p list_R2.txt`
        count=`expr $count + 1`
        
        clean_lineA="${lineA/illumina_reads/fastq_clean}"
        clean_lineB="${lineB/illumina_reads/fastq_clean}"
		
		bbduk.sh \
			in1=$lineA \
			in2=$lineB \
			out1=$clean_lineA \
			out2=$clean_lineB \
			maq=25 \
			minlen=50 \
			qtrim=lr \
			trimq=25

	done < list_R1.txt		
		
echo '######################'
echo 'Job finished' $(date --iso-8601=seconds)

##############################################################################
								Unicycler_PAS
##############################################################################

module load unicycler/0.4.8

cd /shared/projects/tolerance_mutations_ecoli/breseq_all/data_reference/

unicycler \
			-1 /shared/projects/tolerance_mutations_ecoli/breseq_all/data_reference/illumina_reads/fastq_clean/PAS133_S26_R1_001.fastq.gz \
			-2 /shared/projects/tolerance_mutations_ecoli/breseq_all/data_reference/illumina_reads/fastq_clean/PAS133_S26_R2_001.fastq.gz \
			-l /shared/projects/tolerance_mutations_ecoli/breseq_all/data_reference/nanopore_reads/PAS_long.fastq.gz \
			-o PAS \
			--threads 10 \
			--keep 0

echo '######################'
echo 'Job finished' $(date --iso-8601=seconds)

##############################################################################
							Prokka-annotation
##############################################################################
module load prokka/1.14.6

cd /shared/projects/tolerance_mutations_ecoli/breseq_all/data_reference/

mkdir PAS_annot

prokka --outdir PAS_annot --force --genus Escherichia --species coli ./PAS/assembly.fasta 


echo '######################'
echo 'Job finished' $(date --iso-8601=seconds)