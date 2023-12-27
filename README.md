# Feature-Extraction-Workflow-for-Target-Mutation-Selection
Bioinformatics pipeline for the extraction of all tumor features necessary for the selection of clonal mutations from whole-exome sequencing (WES) data in tumor-normal settings
Employs two publicly accessible tools, Mutect2 for the detection of single nucleotide variants (SNV) and insertion and deletion (indel) variants (available at https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2), and PureCN for the extraction of essential tumor features such as cancer cell fraction (CCF), multiplicity of mutations, tumor sample purity, and ploidy. (available at https://bioconductor.org/packages/devel/bioc/vignettes/PureCN/inst/doc/Quick.html).
***paper link will be available***
### Requirements
Installed conda (available instructions: https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

### Input files
Tumor and normal bam files of queried sample

### Provided here
* A bash script create_interval_and_pon.sh, which is for a one-time run to create necessary files (NormalDB and interval.file) for running sample workflow.py
* workflow.py for creating output for the queried sample

## Quick start

**1. Install conda environment from provided gatk-purecn.yml**
```{bash}
conda env create -f gatk-purecn.yml
```

**2. In the project folder create a sub-folder for generating NormalDB and the interval file needed to run the pipeline**

* Copy provided one-time run bash script
* Download GRCh38 reference file (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz)
* Download mappability file for GRCh38 (GCA_000001405.15_GRCh38_no_alt_analysis_set_76.bw from https://bioconductor.org/packages/devel/bioc/vignettes/PureCN/inst/doc/Quick.html)
* Change paths and run the script
```{bash}
./create_interval_and_pon.sh
```
**3.To set up Slurm back-end, run the following command (described at https://gwf.app/guide/tutorial)**

```{bash}
gwf config set backend slurm
```

**4. After successfully generating of necessary files, run gwf workflow.py for the queried sample(s) after adjusting paths within workflow.py file**

```{bash}
gwf run
``` 
