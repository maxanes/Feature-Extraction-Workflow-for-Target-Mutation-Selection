import os
from gwf import Workflow
from gwf import AnonymousTarget

### TEMPLATES ###

def run_mutect(normal_bam, tumor_bam, config, genome_fasta, db_snps, output_folder):
    """
    Run Mutect2 in tumor with matched normal settings according to GATK best practices,
    filter VCF, and annotate with DB SNP annotation, which is an input for PureCN tool
    """

    inputs = [normal_bam, tumor_bam, genome_fasta, db_snps]
    normal_sample = config['sample-names']['normal']
    tumor_sample = config['sample-names']['tumor']

    outputs = {
        'vcf_file': os.path.join(output_folder, f'{tumor_sample}.mutect.vcf.gz'),
        'vcf_file_index': os.path.join(output_folder, f'{tumor_sample}.mutect.vcf.gz.tbi'),
        'stats_file': os.path.join(output_folder, f'{tumor_sample}.mutect.vcf.gz.stats'),
        'f1r2_stats_file': os.path.join(output_folder, f'{tumor_sample}.f1r2.tar.gz'),
        'model': os.path.join(output_folder, f'{tumor_sample}.read-orientation-model.tar.gz'),
        'filtered_vcf': os.path.join(output_folder, f'{tumor_sample}.mutect.filtered.vcf.gz'),
        'filtered_vcf_index': os.path.join(output_folder, f'{tumor_sample}.mutect.filtered.vcf.gz.tbi'),
        'filtered_annot_vcf': os.path.join(output_folder, f'{tumor_sample}.mutectDB.filtered.vcf.gz'),
        'filtered_annot_vcf_index': os.path.join(output_folder, f'{tumor_sample}.mutectDB.filtered.vcf.gz.tbi')
    }

    options = dict(cores='20', memory='32g', walltime='65:00:00')

    spec = f'''
    set -e

    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    source $CONDA_BASE/bin/activate gatk-purecn

    # Run Mutect2
    gatk Mutect2 \
        -R {genome_fasta} \
        -I {tumor_bam} \
        -I {normal_bam} \
        -normal {normal_sample} \
        --genotype-germline-sites true \
        --genotype-pon-sites true \
        --f1r2-tar-gz {outputs['f1r2_stats_file']} \
        -O {outputs['vcf_file']}

    # Make statistics and model about fw/rv reads
    gatk LearnReadOrientationModel \
        -I {outputs['f1r2_stats_file']} \
        -O {outputs['model']}

    # Insert in Filter of mutations PASS/NO
    gatk FilterMutectCalls -V {outputs['vcf_file']} \
        -R {genome_fasta} \
        --ob-priors {outputs['model']} \
        -O {outputs['filtered_vcf']}

    ## Annotate with DB flag
    #CONDA_BASE=$(conda info --base)
    #source $CONDA_BASE/etc/profile.d/conda.sh
    #source $CONDA_BASE/bin/activate bcftools

    bcftools annotate --annotation {db_snps} --columns ID,INFO/DB --output {outputs['filtered_annot_vcf']} {outputs['filtered_vcf']}
    tabix {outputs['filtered_annot_vcf']}
    '''

    return AnonymousTarget(
        inputs=inputs,
        outputs=outputs,
        options=options,
        spec=spec)

def run_purecn(tumor_bam, pon, purecn_vcf, intervals, config, output_folder):
    """
    This runs PureCN tool according to PureCN best practices
    for the estimation of tumor purity, ploidy, LOH, and cancer cell fraction of mutations
    in tumor with matched normal settings
    """

    inputs = [tumor_bam, pon, purecn_vcf, intervals]
    normal_sample = config['sample-names']['normal']
    tumor_sample = config['sample-names']['tumor']
    
    outputs = {
        'tumor_coverage': os.path.join(output_folder, f'{tumor_sample}.tumor.alignment_coverage_loess.txt.gz')
    }

    options = dict(cores='1', memory='50g', walltime='05:00:00')

    spec = f'''
    set -e

    export PURECN=$(conda env list | grep "*" | rev | cut -f1 -d' ' | rev)/lib/R/library/PureCN/extdata

    # Calculate and GC-normalize coverage from a tumor BAM file
    Rscript ${{PURECN}}/Coverage.R \
        --out-dir {output_folder} \
        --bam {tumor_bam} \
        --intervals {intervals}

    # With a matched normal tumor sample
    Rscript ${{PURECN}}/PureCN.R \
        --out {output_folder} \
        --tumor {outputs['tumor_coverage']} \
        --normaldb {pon} \
        --sampleid {tumor_sample} \
        --vcf {purecn_vcf} \
        --intervals {intervals} \
        --fun-segmentation PSCBS  \
        --genome hg38 \
        --model betabin \
        --force \
        --post-optimize \
        --seed 42
    '''

    return AnonymousTarget(
        inputs=inputs,
        outputs=outputs,
        options=options,
        spec=spec)

### WORKFLOW ###

gwf = Workflow(defaults={'account': 'genomic-struct'})

config = {
    "normal-bam": "/faststorage/project/genomic-struct/Workspaces/mnesic/WEStest/PureCN/test_pipeline_paper/C45A10589D_buffycoat.normal.alignment.bam", #Path/to/normal_bam_file
    "tumor-bam": "/faststorage/project/genomic-struct/Workspaces/mnesic/WEStest/PureCN/test_pipeline_paper/C23T10589D_ffpe.tumor.alignment.bam",#Path/to/tumor_bam_file
    "genome-fasta": "/faststorage/project/MomaReference/BACKUP/hg38/reference_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.GRC_exclusions_masked.fna", #Path/to/fasta_file
    "db-snps": "/home/mnesic/PureCN/extrafiles/1000G_phase1.snps.high_confidence.hg38.vcf.gz", #Path/to/DB_file
    "pon": "/faststorage/project/MomaReference/BACKUP/hg38/research/purecn/pon/normalDB_twist_hg38.rds", #Path/to/created_NormalDB_file
    "target-intervals": "/faststorage/project/MomaReference/BACKUP/hg38/research/purecn/intervals/baits_hg38_intervals.txt", #Path/to/created_interval_file
    "sample-names": {
        "tumor": "C23T10589D_ffpe",
        "normal": "C45A10589D_buffycoat"
    }
}

output_folder = os.path.join('output')
os.makedirs(output_folder, exist_ok=True)

# Define the targets using corrected function calls
gwf.target_from_template(
    "mutect2",
    run_mutect(normal_bam=config["normal-bam"],
               tumor_bam=config["tumor-bam"],
               genome_fasta=config["genome-fasta"],
               db_snps=config["db-snps"],
               output_folder=output_folder,
               config=config)
)


tumor_sample = config['sample-names']['tumor']

gwf.target_from_template(
    "purecn",
    run_purecn(tumor_bam=config["tumor-bam"],
               intervals=config["target-intervals"],
               pon=config["pon"],
               purecn_vcf = os.path.join(output_folder, f'{tumor_sample}.mutectDB.filtered.vcf.gz'),
               config=config,
               output_folder=output_folder)
)
