#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu=60G
#SBATCH --account=genomic-struct #Namae of your account
#SBATCH  --time=03:00:00

#WHAT: run purecn to create coverage files of normal and pon 
#Make coverage file for normal

set -e

#Activate conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
source $CONDA_BASE/bin/activate gatk-purecn

#export path to the folder with PCN scripts
export PURECN=$(conda env list | grep "*" | rev | cut -f1 -d' ' | rev)/lib/R/library/PureCN/extdata

#Paths to the downloaded
OUTDIR="/home/mnesic/PureCN/test_pipeline_paper/PON_test/retest" #Path/to/Output_folder
fasta_file="/home/mnesic/genomic-struct/Workspaces/mnesic/WEStest/PureCN/extrafiles/hg38.fa" #Path/to/downloaded/fasta_file
map_file="/home/mnesic/genomic-struct/Workspaces/mnesic/WEStest/PureCN/extrafiles/GCA_000001405.15_GRCh38_no_alt_analysis_set_76.bw" #Path/to/downloaded/mapability_file
capture_bed_file="/faststorage/project/genomic-struct/Workspaces/mnesic/WEStest/PureCN/extrafiles/exome_captures_improve/exome_v1_primary_targets.interval.bed" #Path/to/your/capture/bed_file

#Create interval file
##Recommendation for aech used kit create separate interval_file

Rscript $PURECN/IntervalFile.R \
    --in-file $capture_bed_file \
    --fasta $fasta_file\
    --out-file $OUTDIR/baits_hg38_intervals.txt \
    --off-target --genome hg38 \
    --export $OUTDIR/baits_optimized_hg38.bed \
    --mappability $map_file \
    --force   
	
	
#Iterate through folders with bam files 
for folder in /home/mnesic/WES_data/231001*/; do #Path/to/directorium with bam files
  # Find the BAM file in the current folder
    bam_file=$(find "$folder" -type f -name "*_buffycoat.normal.alignment.bam")
  
  if [ -f "$bam_file" ]; then
    echo "Processing BAM file: $bam_file"
    
    # Run the Rscript with the specified parameters
    Rscript $PURECN/Coverage.R \
      --out-dir $OUTDIR \
      --bam $bam_file \
      --intervals $OUTDIR/baits_hg38_intervals.txt \
      --force
  else
    echo "No BAM file found in folder: $folder"
  fi
done

#Create pon
##Create pon for each kit used, do not mix different kits  
ls -a /home/mnesic/PureCN/test_pipeline_paper/PON_test/retest/*_loess.txt.gz | cat > normal_coverages.list #Path/to/created/coverage_files/*

Rscript $PURECN/NormalDB.R --out-dir "$OUTDIR" \
--coverage-files normal_coverages.list \
--genome hg38 \
--assay twist
