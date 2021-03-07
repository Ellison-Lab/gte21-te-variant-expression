rule bigwigs:
    input:
        bam = rules.get_te_reads.output.bam,
        bai = rules.get_te_reads.output.bai
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
