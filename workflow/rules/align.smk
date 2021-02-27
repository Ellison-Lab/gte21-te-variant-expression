rule get_trna_miscrna:
    output:
        "results/bait/bait.fasta"
    params:
        trna = config.get("TRNA_FASTA"),
        misc = config.get("MISC_RNA_FASTA"),
    shell:
        "wget -q -O - {params.trna} | zcat > {output} && "
        "wget -q -O - {params.misc} | zcat >> {output}"

rule star_bait_index:
    input:
        fasta = rules.get_trna_miscrna.output
    output:
        directory("results/idx/bait")
    threads:
        12
    params:
        extra = "--genomeSAindexNbases 6"
    log:
        "results/logs/star_bait_index.log"
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/star/index"

rule star_genome_index:
    input:
        fasta = custom_genome('results/custom-genome/combined.fasta'),
        gtf = custom_genome('results/custom-genome/combined.fixed.gtf')
    output:
        directory("results/idx/genome")
    threads:
        12
    params:
        extra = "--genomeSAindexNbases 12"
    resources:
        time=60,
        mem=20000,
        cpus=12
    log:
        "results/logs/star_genome_index.log"
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/star/index"

def get_fqr1(wc):
    tmp = SUBSAMPLE_TABLE[SUBSAMPLE_TABLE['sample_name'] == wc.sample]
    tmp2 = tmp[tmp['subsample_name'] == wc.subsample]
    return tmp2.get('fastq_r1')

def get_fqr2(wc):
    tmp = SUBSAMPLE_TABLE[SUBSAMPLE_TABLE['sample_name'] == wc.sample]
    tmp2 = tmp[tmp['subsample_name'] == wc.subsample]
    return tmp2.get('fastq_r2')

rule star_align_bait:
    input:
        fq1 = get_fqr1,
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = get_fqr2,
        idx = rules.star_bait_index.output,
    output:
        # see STAR manual for additional output files
        "results/aln/bait/{sample}/{subsample}/Aligned.out.bam",
        r1 = "results/aln/bait/{sample}/{subsample}/Unmapped.out.mate1",
        r2 = "results/aln/bait/{sample}/{subsample}/Unmapped.out.mate2",
    log:
        "results/logs/star/bait/{sample}/{subsample}.log"
    params:
        # path to STAR reference genome index
        index="results/idx/bait",
        # optional parameters
        extra="--outReadsUnmapped Fastx --outSAMtype BAM Unsorted"
    resources:
        time=360,
        mem=128000,
        cpus=24
    threads:
        4
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/star/align"


rule star_align_genome:
    input:
        fq1 = rules.star_align_bait.output.r1,
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = rules.star_align_bait.output.r2,
        idx = rules.star_genome_index.output,
        vcf = te_variants('results/snps/snps.vcf'),
        gtf = custom_genome('results/custom-genome/combined.fixed.gtf')
    output:
        # see STAR manual for additional output files
        "results/aln/genome/{sample}/{subsample}/Aligned.sortedByCoord.out.bam"
    log:
        "results/logs/star/genome/{sample}/{subsample}.log"
    resources:
        time=360,
        mem=128000,
        cpus=24
    params:
        # path to STAR reference genome index
        index="results/idx/genome",
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1 --outSAMattributes NH HI AS nM vA vG --varVCFfile {v} --sjdbGTFfile {g}".format(g=custom_genome('results/custom-genome/combined.fixed.gtf'), v=te_variants('results/snps/w1118_male-snps.vcf'))
    threads:
        24
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/star/align"

rule filter_reads:
    input:
        rules.star_align_genome.output
    output:
        "results/aln/filt/{sample}/{subsample}.bam"
    resources:
        time=20,
        mem=10000,
    params:
        "-O BAM -F 256"
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/samtools/view"
