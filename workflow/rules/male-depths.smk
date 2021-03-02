rule get_te_reads:
    input:
        rules.filter_reads.output
    output:
        tmp = temp("results/aln/tes/{sample}/{subsample}.tmp.bam"),
        tmpbai = temp("results/aln/tes/{sample}/{subsample}.tmp.bam.bai"),
        bam = "results/aln/tes/{sample}/{subsample}.tes.bam",
        bai = "results/aln/tes/{sample}/{subsample}.tes.bam.bai",
    params:
        tes = TES
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        cp {input} {output.tmp} && samtools index {output.tmp}
        samtools view -b {output.tmp} {params.tes} > {output.bam} && \
        samtools index {output.bam}
        """

rule get_male_specific_te_reads:
    input:
        rules.get_te_reads.output.bam
    output:
        sam = temp("results/aln/sex-te/{sample}/{subsample}.male.sam"),
        bam = "results/aln/sex-te/{sample}/{subsample}.male.bam",
        bai = "results/aln/sex-te/{sample}/{subsample}.male.bam.bai"
    params:
        tes = TES
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -H {input} > {output.sam} && \
        samtools view {input} | grep "vA:B" | grep ",2[,\\s]" >> {output.sam} && \
        samtools view -b {output.sam} > {output.bam} && \
        samtools index {output.bam}
        """

rule get_non_male_te_reads:
    input:
        rules.get_te_reads.output.bam
    output:
        sam = temp("results/aln/sex-te/{sample}/{subsample}.unknown.sam"),
        bam = "results/aln/sex-te/{sample}/{subsample}.unknown.bam",
        bai = "results/aln/sex-te/{sample}/{subsample}.unknown.bam.bai"
    params:
        tes = TES
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -H {input} > {output.sam} && \
        samtools view {input} | grep -v "vA:B:\\w+,2[,\\s]" | grep -v "c,2" >> {output.sam} && \
        samtools view -b {output.sam} > {output.bam} && \
        samtools index {output.bam}
        """

rule get_male_depth_per_snp:
    input:
        bam = rules.get_te_reads.output.bam,
        vcf = te_variants('results/snps/snps.vcf'),
    output:
        csv = "results/depths/{sample}-{subsample}-depth-at-male-snps.csv.gz",
    resources:
        time=20,
        mem=10000,
        cpus=2
    conda:
        "../envs/bioc-general.yaml"
    script:
        '../scripts/depth-at-snps.R'

rule get_male_reads_per_snp:
    input:
        male = rules.get_male_specific_te_reads.output.bam,
        unk = rules.get_non_male_te_reads.output.bam,
        vcf = te_variants('results/snps/snps.vcf'),
    output:
        csv = "results/depths/{sample}-{subsample}-reads-at-male-snps.csv.gz",
    resources:
        time=20,
        mem=10000,
        cpus=2
    conda:
        "../envs/bioc-general.yaml"
    script:
        '../scripts/reads-at-snps.R'
