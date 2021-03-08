rule bigwigs:
    input:
        bam = "results/aln/filt/{sample}/{subsample}.bam",
        bai = "results/aln/filt/{sample}/{subsample}.bam.bai"
    output:
        "results/bigwigs/tes/{sample}.{subsample}.tes.strand-{strand}.rpkm.bw"
    threads:
        24
    conda:
        "../envs/deeptools.yaml"
    resources:
        time=60,
        mem=20000,
        cpus=24
    shell:
        """
        bamCoverage -b {input.bam} \
            -o {output} -p {threads} -v --normalizeUsing RPKM \
            --smoothLength 150 --filterRNAstrand {wildcards.strand}
        """
