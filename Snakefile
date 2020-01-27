import os
from glob import glob
import yaml
import itertools

configfile: "/projectnb/bradham/indrops/2019-10-24/2019-10-24-config.yaml"
shell.prefix("source activate indrops; PYTHONWARNINGS=ignore::yaml.YAMLLoadWarning; ")

with open(config['yaml'], 'r') as f:
    project_config = yaml.load(f, yaml.FullLoader)

def symlink_output(input_dir, output_dir):
    files = glob(os.path.join(input_dir, '*.fastq.*'))
    output_files = [os.path.join(output_dir, os.path.split(x)[-1])\
                    for x in files]
    return output_files

PROJECT_DIR = project_config['project_dir']
BASECALLS = os.path.join(config['base_dir'], 'Data', 'Intensities', 'BaseCalls')
RUN_LIBRARIES = []
LIBRARIES = {}
for each in project_config['sequencing_runs']:
    run = each['name']
    for lib in each['libraries']:
        RUN_LIBRARIES.append((run, lib['library_name']))
        LIBRARIES[lib['library_name']] = lib['library_index']
WORKERS = range(config['cores'])
SPLITS = ['L001', 'L002', 'L003', 'L004']
READS = ['R1', 'R2', 'R3', 'R4']
# expected output from calling bcl2fastq
MULTIPLEXED = [os.path.join(BASECALLS,
                            "Undetermined_S0_{split}_{R}_001.fastq.gz".format(
                            split=split, R=R))\
               for split, R in itertools.product(SPLITS, READS)]
# new symlinked locations for fastq
SYMLINKS = [os.path.join(config['base_dir'], 'fastq', os.path.split(x)[-1])\
            for x in MULTIPLEXED]

def sorted_output(library):
    out = os.path.join(PROJECT_DIR, library, 'filtered_parts',
                       "{library}_{run}_{index}_{split}.fastq.sorted.fastq.gz")
    for run, lib in RUN_LIBRARIES:
        if lib == library:
            return([out.format(library=lib, run=run, index=LIBRARIES[lib],
                               split=split) for split in SPLITS])

def aggregate_input(wildcards):
    library_quant = [os.path.join(PROJECT_DIR, wildcards.library, "quant_dir",
                     "worker{i}_".format(i=i) + str(config['cores']) + ".counts.tsv") \
                    for i in WORKERS]
    return library_quant

rule all:
    input:
        '{}.transcripts.fa'.format(config['bowtie_index']),
        [os.path.join(PROJECT_DIR, "{library}/{library}.filtering_stats.csv".\
                      format(library=library)) for library in LIBRARIES.keys()],
        [os.path.join(config['base_dir'], config['run_id'],
                     '{library}_{worker}_quantify_reads.out'.format(library=library,
                                                                    worker=worker))\
         for library, worker in itertools.product(LIBRARIES.keys(), WORKERS)],
        [os.path.join(config['base_dir'], config['run_id'],
                     '{library}_extract_barcodes.out'.format(library=library))\
         for library in LIBRARIES.keys()],
        [os.path.join(PROJECT_DIR, "{library}/{library}.meta.csv".\
                       format(library=library))\
            for library in LIBRARIES.keys()]      

rule extract_fastqs:
    output:
        MULTIPLEXED
    params:
        base_dir=config['base_dir']
    shell:
        """
        cd {params.base_dir} 
        module load bcl2fastq
        bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0 --minimum-trimmed-read-length 0
        """

rule symlink_fastq_files:
    input:
        MULTIPLEXED,
    output:
        SYMLINKS
    params:
        link_dir=os.path.join(config['base_dir'], 'fastq')
    shell:
        "ln -s {input} {params.link_dir}"
    
rule build_bowtie_index:
    input:
        fasta=config['genome_fa'],
        gtf=config['genome_gtf'],
        yaml=ancient(config['yaml'])
    params:
        indrops=os.path.join(config['indrops_dir'], 'indrops.py')
    output:
        out='{}.transcripts.fa'.format(config['bowtie_index'])
    shell:
        "python {params.indrops} {input.yaml} build_index "
        "--genome-fasta-gz {input.fasta} --ensembl-gtf-gz {input.gtf}"

rule filter_reads:
    input:
        index='{}.transcripts.fa'.format(config['bowtie_index']),
        symlinks=SYMLINKS,
        yaml=ancient(config['yaml'])
    output:
        os.path.join(config['base_dir'], config['run_id'],
                     "{library}_{run}_{worker}_filter.out")
    #     [os.path.join(PROJECT_DIR, '{library}',
    #     'filtered_parts', "{library}_{run}_{index}_" + "{split}.fastq".format(
    #         split=split)) for split in SPLITS]
    params:
        workers=config['cores']
    log:
        "logs/{run}_{library}_{worker}_filter.log"
    shell:
        """
        if (python indrops.py {input.yaml} filter --runs {wildcards.run} --libraries {wildcards.library} --total-workers {params.workers} --worker-index {wildcards.worker}) > {log}; then
            echo "Filtered successful!" > {output}
        else
            echo "Filtered Failed!"
        fi
        """

rule abundant_barcodes:
    input:
        [os.path.join(config['base_dir'], config['run_id'],
                     "{library}_{run}_{worker}_filter.out".format(library=run_lib[1],
                                                           run=run_lib[0],
                                                           worker=worker))\
        for run_lib, worker in itertools.product(RUN_LIBRARIES, WORKERS)],
        yaml=ancient(config['yaml'])
    output:
        os.path.join(PROJECT_DIR, "{library}", "{library}.filtering_stats.csv")
    shell:
        """
        python indrops.py {input.yaml} identify_abundant_barcodes --libraries {wildcards.library}
        """

rule sort_reads:
    input:
        os.path.join(PROJECT_DIR, "{library}", "{library}.filtering_stats.csv"),
        yaml=ancient(config['yaml'])
    output:
        os.path.join(config['base_dir'], config['run_id'],
                     '{library}_sort_reads.out')
    log:
        "logs/{library}_sort.log"
    shell:
        """
        if (python indrops.py {input.yaml} sort --libraries {wildcards.library}) 2> {log}; then
            echo "Reads sorted!" > {output}
        else
            echo "Fail!"
        fi
        """

rule extract_barcodes:
    input:
        reads=os.path.join(config['base_dir'], config['run_id'],
                           '{library}_sort_reads.out'),
        yaml=ancient(config['yaml'])
    output:
        # flag=os.path.join(config['base_dir'], config['run_id'],
        #                   '{library}_extract_barcodes.out'),
        os.path.join(PROJECT_DIR, "{library}", "barcode_fastq",
                     "{library}.{barcode}.fastq")
    # log:
    #     os.path.join(config['base_dir'], "logs",
    #                  "{libary}_extract_barcodes.log")
    shell:
        """
        python indrops.py {input.yaml} output_barcode_fastq --libraries {wildcards.library};       fi
        """

rule star_align_barcodes:
    input:
        fastq=os.path.join(PROJECT_DIR, "{library}", "barcode_fastq",
                           "{library}.{barcode}.fastq"),
        fasta=config['genome_fa'],
        gtf=config['genome_gtf'],
    params:
        prefix=os.path.join(PROJECT_DIR, "{library}", "barcode_fastq",
                            "{library}", "bams", "{barcode}",
                            "{library}.{barcode}")
    output:
        os.path.join(PROJECT_DIR, "{library}", "barcode_fastq",
                     "{library}", "bams", "{barcode}",
                     "{library}.{barcode}Aligned.bam")
    params:
        index="/projectnb/bradham/data/ReferenceSequences/GenomeAnnotations/indrops_star",
    shell:
        "STAR --runMode alignReads --outSAMtype BAM Unsorted --readFilesCommand "
        "cat --genomeDir {params.index} --readFilesIn {fastq} "
        "--outFileNamePrefix {output}"


# def aggregate_barcodes(wildcards):
#     '''
#     aggregate the file names of the random number of files
#     generated at the scatter step
#     '''
#     checkpoint_output = checkpoints.extract_barcodes.get(**wildcards).output[0]
#     return expand('scatter_copy_head/{i}_head.txt',
#            i=glob_wildcards(os.path.join(checkpoint_output, '{i}.txt')).i)

# rule aggreate_alignments:
#     input:



rule quantify_barcodes:
    input:
        os.path.join(config['base_dir'], config['run_id'], '{library}_sort_reads.out'),
        yaml=ancient(config['yaml']),
        # This only necessary to allow worker
        # filtered=os.path.join(config['base_dir'], config['run_id'],
        #                       "{library}_{run}_{worker}_filter.out")
    output:
        os.path.join(PROJECT_DIR, "{library}", "quant_dir", "worker{worker}_"\
                     + str(config['cores']) + ".counts.tsv")
    params:
        workers=config['cores']
    shell:
        """
        python indrops.py {input.yaml} quantify --libraries {wildcards.library} --total-workers {params.workers} --worker-index {wildcards.worker}
        """

rule aggregate_umis:
    input:
        lambda wildcards: aggregate_input(wildcards),
        # [os.path.join(PROJECT_DIR, "{library}/quant_dir/" "worker{i}_".format(i=i) \
        #                             + str(config['cores']) + ".counts.tsv") \
        #             for i in WORKERS]
        yaml=ancient(config['yaml'])
    params:
        workers=config['cores']
    output:
        os.path.join(PROJECT_DIR, "{library}/{library}.counts.tsv.gz")
    shell:
        """
        python indrops.py {input.yaml} aggregate --total-workers {params.workers} --libraries {wildcards.library}
        """

rule unzip_counts:
    input:
        os.path.join(PROJECT_DIR, "{library}/{library}.counts.tsv.gz")
    output:
        os.path.join(PROJECT_DIR, "{library}/{library}.counts.tsv")
    shell:
        "gunzip {input}"

rule annotate_cells:
    input:
        counts=os.path.join(PROJECT_DIR, "{library}/{library}.counts.tsv"),
        json = os.path.join(config['run_id'], 'library_meta.json')
    output:
        csv=os.path.join(PROJECT_DIR, "{library}/{library}.meta.csv")
    script:
        "scripts/annotate_libraries.py"


    

