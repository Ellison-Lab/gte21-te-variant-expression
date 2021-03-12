library(VariantAnnotation)
library(Rsamtools)
library(tidyverse)

bam <- snakemake@input[['bam']]
#bam <- "results/aln/tes/w1118_testes/rep1.tes.bam"

vcf <- snakemake@input[['vcf']]
#vcf <- "~/work/transposon-variants-hts/results/snps/snps.vcf"

sample_name <- snakemake@wildcards[['sample']]
# sample_name <- 'aa'
subsample_name <- snakemake@wildcards[['subsample']]
# subsample_name <- 'bb'

vr <- VariantAnnotation::readVcfAsVRanges(vcf)

gr <- vr %>% GRanges()

allele.lookup <- vr %>% as.data.frame() %>% as_tibble() %>%
  dplyr::select(seqnames, pos=start, ref, alt, specificity) %>%
  filter(specificity == 'w1118_male')


pup <- PileupParam(max_depth = 10e6,
                   distinguish_nucleotides = T,
                   distinguish_strands = F,
                   min_minor_allele_depth = 0,
                   min_nucleotide_depth = 1)

scbp <- ScanBamParam(which = gr)

pileups.df <- pileup(bam, scanBamParam = scbp, pileupParam = pup) %>%
  as_tibble(.) %>%
  mutate(sample = sample_name) %>%
  mutate(subsample = subsample_name)

pileups.df <- pileups.df %>%
  left_join(allele.lookup, by=c('seqnames','pos')) %>%
  mutate(sex = ifelse(nucleotide == alt, 'male','unknown')) %>%
  group_by(seqnames, pos, sex, sample, subsample) %>%
  #summarise(depth = sum(`count`),.groups = 'drop') %>%
  dplyr::select(sample, subsample, seqnames,pos,sex, count) %>%
  arrange(seqnames, pos)


write_csv(pileups.df, snakemake@output[['csv']])
