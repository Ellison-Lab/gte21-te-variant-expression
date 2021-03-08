rule bigwigs:
    input:
        bam = "results/aln/filt/{sample}/{subsample}.bam",
        bai = "results/aln/filt/{sample}/{subsample}.bam.bai"
    output:
        "results/bigwigs/tes/{sample}.{subsample}.tes.strand-{strand}.rpkm.bw"
    threads:
        4
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        bamCoverage -b {input.bam} \
            -o {output} -p {threads} -v --normalizeUsing RPKM \
            --smoothLength 150 --filterRNAstrand {wildcards.strand}
        """
