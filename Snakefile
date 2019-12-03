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

BASECALLS = os.path.join(config['base_dir'], 'Data', 'Intensities', 'BaseCalls')
RUN_LIBRARIES = []
for each in project_config['sequencing_runs']:
    run = each['name']
    for lib in each['libraries']:
        RUN_LIBRARIES.append((run, lib['library_name']))

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

rule all:
    input:
        [os.path.join(config['indrops_dir'], config['run_id'],
                      "{run}_{library}_filtered.out".format(run=run,
                                                            library=library))\
         for run, library in RUN_LIBRARIES]

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
        os.path.join(config['indrops_dir'], config['run_id'],
                            "{run}_{library}_filtered.out")
    shell:
        """
        if python indrops.py {input.yaml} filter --runs {wildcards.run} --libraries {wildcards.library}; then 
            echo "Filtering complete!" > {output}
        else
            echo "Exit code $?, failure"
        fi
        """
    
