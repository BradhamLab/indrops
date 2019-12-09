import os
from glob import glob
import yaml
import itertools

configfile: "/projectnb/bradham/indrops/2019-10-24/2019-10-24-config.yaml"
shell.prefix("source activate indrops; ")

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

rule all:
    input:
        [os.path.join(PROJECT_DIR, '{library}'.format(library=run_lib[1]),
        'filtered_parts', "{library}_{run}_{index}_{split}.fastq".format(
            library=run_lib[1], run=run_lib[0], index=LIBRARIES[run_lib[1]], split=split))\
        for run_lib, split in itertools.product(RUN_LIBRARIES, SPLITS)],
        ["logs/{run}_{library}_abundant.log".format(run=run, library=library)\
        for run, library in RUN_LIBRARIES],
        [os.path.join(config['base_dir'], config['run_id'],
                      '{library}_sort_reads.out'.format(library=library))\
        for __, library in RUN_LIBRARIES],
        [os.path.join(config['base_dir'], config['run_id'],
                     '{library}_quantify_reads.out'.format(library=library))\
        for __, library in RUN_LIBRARIES]
        

rule extract_fastqs:
    output:
        MULTIPLEXED
    params:
        base_dir=config['base_dir']
    shell:
        """
        cd {params.base_dir}; 
        module load bcl2fastq;
        bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0
        --minimum-trimmed-read-length 0
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
        [os.path.join(PROJECT_DIR, '{library}',
        'filtered_parts', "{library}_{run}_{index}_" + "{split}.fastq".format(
            split=split)) for split in SPLITS]
    log:
        "logs/{run}_{library}_{index}_filter.log"
    shell:
        """
        (python indrops.py {input.yaml} filter --runs {wildcards.run} --libraries {wildcards.library}) > {log}
        """

rule abundant_barcodes:
    input:
        [os.path.join(PROJECT_DIR, '{library}'.format(library=run_lib[1]),
        'filtered_parts', "{library}_{run}_{index}_{split}.fastq".format(
            library=run_lib[1], run=run_lib[0], index=LIBRARIES[run_lib[1]], split=split))\
        for run_lib, split in itertools.product(RUN_LIBRARIES, SPLITS)],
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
        os.path.join(config['base_dir'], config['run_id'], '{library}_sort_reads.out')
    log:
        "logs/{library}_sort.log"
    shell:
        """
        if (python indrops.py {input.yaml} sort --libraries {wildcards.library}) > {log}; then
            echo "Reads sorted!" > {output}
        else
            echo "Fail!"
        fi
        """

rule quantify_barcodes:
    input:
        os.path.join(config['base_dir'], config['run_id'], '{library}_sort_reads.out'),
        yaml=ancient(config['yaml'])
    output:
        os.path.join(config['base_dir'], config['run_id'],
                     '{library}_quantify_reads.out')
    log:
        "logs/{library}_quantify.log"
    shell:
        """
        if python indrops.py {input.yaml} quantify --libraries {wildcards.library}; then
            echo "Reads quantified!" > {output}
        else
            echo "Fail!"
        fi
        """