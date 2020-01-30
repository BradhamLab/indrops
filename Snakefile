import os
from glob import glob
import yaml
import itertools
import subprocess as sbp
import numpy as np

configfile: "/projectnb/bradham/indrops/2019-10-24/2019-10-24-config.yaml"
shell.prefix("source activate indrops; PYTHONWARNINGS=ignore::yaml.YAMLLoadWarning; ")

with open(config['yaml'], 'r') as f:
    project_config = yaml.load(f, yaml.FullLoader)

def symlink_output(input_dir, output_dir):
    files = glob(os.path.join(input_dir, '*.fastq.*'))
    output_files = [os.path.join(output_dir, os.path.split(x)[-1])\
                    for x in files]
    return output_files

# function to get genomeChrBinNBits parameter for STAR alignment.
def estimate_STAR_ChrBinNbits(genome_file, read_length):
    """
    Estimate the `ChrBinNBits` parameter for genome indexing in STAR
    Estimate the `ChrBinNBits` parameter for genome indexing in STAR. Value
    must be estimated due to memory constraints caused by the large number
    of scaffolds present in some genomes (i.e. the LV genome). If estimation
    is unnecessary, flag `star_est_ChrBinNbits: False` in configuration file.
    Args:
        genome_file (string): path to fasta file containing genome reference
            sequences.
        read_length (int): length of reads from RNAseq experiment.
    Return:
        (int) new value for scaling RAM consumption
    References:
    https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf (p. 7)
    https://github.com/alexdobin/STAR/issues/103
    """
    len_call = 'grep -v ">" {} | wc | awk '.format(genome_file)\
               + "'{print $3-$1}'"
    n_ref_call = 'grep "^>" {} | wc -l'.format(genome_file)

    return_values = [None, None]
    for i, call in enumerate([len_call, n_ref_call]):
        p = sbp.Popen(call, stdin=sbp.PIPE, stdout=sbp.PIPE, stderr=sbp.PIPE,
                      shell=True)
        output, err = p.communicate()
        if p.returncode == 0:
            return_values[i] = int(output.strip())
        else:
            raise OSError(err)
    estimate = max([int(np.log2(return_values[0] / return_values[1])),
                    int(np.log2(read_length))])
    return min(18, estimate)

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
        # '{}.transcripts.fa'.format(config['bowtie_index']),
        os.path.join(config['star_index'], 'Genome'),
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

rule fastqc_biological_reads:
    input:
        [x in SYMLINKS if 'R1' in x]
    params:
        outdir=os.path.join(config['base_dir'], config['run_id'], 'fastqc')
    output:
        [os.path.join(config['base_dir'], config['run_id'],
                      'fastqc', x.replace('.fastq.gz', '_fastqc.html'))\
         for x in SYMLINKS if 'R1' in x]
    shell:
        "fastqc {input} -o {params.outdir}"
        

if not os.path.exists(config['star_index']):
    os.makedirs(config['star_index'])
    

rule build_star_index:
    input:
        fasta=config['genome_fa'],
        gff3=config['genome_gtf'],
        yaml=ancient(config['yaml'])
    params:
        index_dir=config['star_index'],
        chr_n_bits=estimate_STAR_ChrBinNbits(config['genome_fa'], 56),
    output:
        os.path.join(config['star_index'], 'Genome')
    shell:
        "STAR --runMode genomeGenerate --genomeDir {params.index_dir} "
        "--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gff3} "
        "--sjdbOverhang 56 --genomeChrBinNbits {params.chr_n_bits} "
        "--sjdbGTFtagExonParentTranscript Parent"


rule filter_reads:
    input:
        index=os.path.join(config['star_index'], 'Genome'),
        symlinks=SYMLINKS,
        yaml=ancient(config['yaml'])
    output:
        os.path.join(config['base_dir'], config['run_id'],
                     "{library}_{run}_{worker}_filter.out")
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
        flag=os.path.join(config['base_dir'], config['run_id'],
                          '{library}_extract_barcodes.out')
    shell:
        """
        if python indrops.py {input.yaml} output_barcode_fastq --libraries {wildcards.library}; then
            echo "Barcode FASTQs extracted" > {output}
        else
            echo "Fail!"
        fi
        """

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


    

