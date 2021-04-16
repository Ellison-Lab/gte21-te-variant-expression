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


rule trim_pe:
    input:
        fq1 = lambda wc: get_fqr1(wc),
        fq2 = lambda wc: get_fqr2(wc),
    output:
        r1 = "results/fastq/{sample}_{subsample}_r1.trimmed.fq.gz",
        r2 = "results/fastq/{sample}_{subsample}_r2.trimmed.fq.gz",
        html = "results/fastq/{sample}_{subsample}_fastp.html",
        json = "results/fastq/{sample}_{subsample}_fastp.json"
    threads:
        12
    resources:
        time=60,
        mem=20000,
        cpus=12
    conda:
        "../envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp --in1 {input.fq1} --in2 {input.fq2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.sample}_fastp"


rule star_align_bait:
    input:
        fq1 = rules.trim_pe.output.r1,
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = rules.trim_pe.output.r2,
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
        time=640,
        mem=128000,
        cpus=24
    threads:
        24
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
        "results/aln/genome/{sample}/{subsample}/Aligned.out.bam"
    log:
        "results/logs/star/genome/{sample}/{subsample}.log"
    resources:
        time=720,
        mem=64000,
        cpus=24
    params:
        # path to STAR reference genome index
        index="results/idx/genome",
        # optional parameters
        extra="--outSAMtype BAM Unsorted --outSAMmultNmax 1 --outSAMattributes NH HI AS nM vA vG --varVCFfile {v} --sjdbGTFfile {g}".format(g=custom_genome('results/custom-genome/combined.fixed.gtf'), v=te_variants('results/snps/snps.vcf'))
    threads:
        24
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/star/align"


rule sort_star:
    input:
        rules.star_align_genome.output
    output:
        "results/aln/sort/{sample}/{subsample}.bam"
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    params:
        extra=""
    resources:
        time=240,
        mem=24000,
        cpus=8
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/samtools/sort"


rule filter_reads:
    input:
        rules.sort_star.output
    output:
        "results/aln/filt/{sample}/{subsample}.bam"
    resources:
        time=20,
        mem=10000,
    params:
        "-O BAM -F 256"
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/samtools/view"

rule samtools_index:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/samtools/index"
