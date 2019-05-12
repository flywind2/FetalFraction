#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=priitpaluoja@gmail.com
#SBATCH --mem=14000
#SBATCH -J FetalFractionCleaner4060
# The job requires 1 compute node
#SBATCH -N 1
# The job requires 1 task per node
#SBATCH --ntasks-per-node=1
# Number of CPU cores per task
#SBATCH --cpus-per-task=10

# Initial flow

# module load bowtie-2.3.2
module load samtools-1.6
module load python-3.6.0

# Build index
# bowtie2-build GCA_000001405.15_GRCh38_genomic.fna bwtie38

for dir in /gpfs/rocket/samba/CCHT/BelgiaNIPT/fastq/*
do
  test -d "$dir" || continue
  echo "$dir - remove old files"
  rm "$dir/filtered.sam"
  rm "$dir/results.sam"
  rm "$dir/converted.csv"
  # Concatenate gz to single gz
  echo "$dir - concatenate"
  cat $dir/*L001* $dir/*L002* $dir/*L003* $dir/*L004* > $dir/concatenated.fastq.gz
  # Use the concatenated for mapping
  echo "$dir - bowtie"
  /gpfs/hpchome/ppaluoja/software/bowtie2-2.3.3.1/bowtie2 --very-sensitive --norc -x /gpfs/rocket/samba/CCHT/BelgiaNIPT/fastq/bwtie38 -q "$dir/concatenated.fastq.gz" -S "$dir/results.sam" --no-unal -p 10
  # Filter by quality
  samtools view -q 35 "$dir/results.sam" > "$dir/filtered35.sam"
  # No need to hold the concatenated file
  echo "$dir - remove concatenated file"
  rm "$dir/concatenated.fastq.gz"
  # Convert positions to chromosomes
  echo "$dir - convert chromosomes"
  python3 convert_chromosomes.py "$dir/filtered35.sam" "$dir/converted35.csv"  
  # Remove file with old positions
  echo "$dir - rm old positions file"
  rm "$dir/results.sam"
  rm "$dir/filtered.sam"
  # Separate to bins
  echo "$dir - separate to bins"
  python3 chromosome_y.py "$dir/converted35.csv" "bins_cleaner_50000_q35_filter_uniq_y.csv" "hglft_genome_f72_b036e0.bed" "50000"
  echo "$dir - rm converted file"
  rm "$dir/converted.csv"
  echo "$dir - DONE!"
  echo " "
done

echo "FILE operations DONE!"

