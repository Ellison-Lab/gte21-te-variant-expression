library(VariantAnnotation)
library(Rsamtools)
library(tidyverse)

male <- snakemake@input[['male']] #male<- 'results/aln/sex-te/w1118_testes/rep1.male.bam'
unk <- snakemake@input[['unk']] # unk <- 'results/aln/sex-te/w1118_testes/rep1.unknown.bam'

vcf <- snakemake@input[['vcf']] #vcf <- "~/work/transposon-variants-hts/results/snps/snps.vcf"
sample_name <- snakemake@wildcards[['sample']]# sample_name <- 'aa'
subsample_name <- snakemake@wildcards[['subsample']]# subsample_name <- 'bb'

gr <-VariantAnnotation::readVcfAsVRanges(vcf) %>% GRanges()

pup <- PileupParam(max_depth = 10e6,
                   distinguish_nucleotides = F,
                   distinguish_strands = F,
                   min_minor_allele_depth = 0,
                   min_nucleotide_depth = 1)

scbp <- ScanBamParam(which = gr)

pileups.df <- list(male=male, unknown =unk) %>%
  map(~pileup(., scanBamParam = scbp, pileupParam = pup)) %>%
  map_df(as_tibble, .id = 'sex') %>%
  mutate(sample = sample_name) %>%
  mutate(subsample = subsample_name) %>%
  dplyr::select(sample, subsample, sex, seqnames,pos,total.depth=count) %>%
  arrange(seqnames, pos)


write_csv(pileups.df, snakemake@output[['csv']])
