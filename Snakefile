import os
configfile: "/projectnb/bradham/indrops/2019-10-24/2019-10-24-config.yaml"
shell.prefix("source activate indrops; ")

rule all:
    input:
        index='{}.transcripts.fa'.format(config['bowtie_index'])

rule demultiplex_bcl:
    output:
        os.path.join(config['indrops_dir'],
                    "{}/demultiplexed.out".format(config['run_id']))
    params:
        run_dir=config['run_dir']
    shell:
        """
        cd {params.run_dir}; 
        module load bcl2fastq;
        bcl2fastq --use-bases-mask y*,y*,y*,y* --mask-short-adapter-reads 0
        --minimum-trimmed-read-length 0;
        echo "Demultiplexing complete!" > {output};
        """


rule build_bowtie_index:
    input:
        fasta=config['genome_fa'],
        gtf=config['genome_gtf'],
        yaml=config['yaml']
    params:
        indrops=os.path.join(config['indrops_dir'], 'indrops.py')
    output:
        out='{}.transcripts.fa'.format(config['bowtie_index'])
    shell:
        "python {params.indrops} {input.yaml} build_index "
        "--genome-fasta-gz {input.fasta} --ensembl-gtf-gz {input.gtf}"
    
