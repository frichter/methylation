# Felix Richter
# Overlap Alex's 850K DMRs with the location of DNVs
##############################################################################


setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

p = c("limma", "IRanges", "GenomicRanges", "edgeR", "variancePartition", "qvalue",
      "stringi", "VennDiagram")
lapply(p, require, character.only = TRUE)
p = c("dplyr", "ggplot2", "tidyr", "readr")
lapply(p, require, character.only = TRUE)

dn = readRDS(file = "whole_genome/pcgc350/dn_all.RDS")

meth = read.table("whole_genome/methylation/Analysis12.List.DMRs.coordinates.custom.1.txt",
                  header = T)

# goal find methylation - de novo overlaps in same individuals
meth = meth %>% mutate(Chr_hg19 = gsub("^", "chr", Chr_hg19)) %>%
  select(-Size, -c(enh.FetalHeart:prom.ahorta))

meth_dna_id_blinded_id_map = read.csv("whole_genome/methylation/meth_dn_id_blinded_id.csv") %>%
  rename(DNA_ID_meth = DNA.Sample.ID) %>%
  select(Blinded.ID, DNA_ID_meth)

# how many surrounding base pairs until you find
# at least 1 de novo mutation with
meth = meth %>% mutate(Start_hg19 = Start_hg19 - 10^6, End_hg19 = End_hg19 + 10^6) %>%
  mutate(Start_hg19 = ifelse(Start_hg19 < 0, 0, Start_hg19))

dn_cases = dn %>% filter(case_ctrl == "Case")
# turn into GRanges object for easier intersections
meth.ir = IRanges(start = meth$Start_hg19, end = meth$End_hg19)
meth.gr = GRanges(seqnames = meth$Chr_hg19, ranges = meth.ir)
mcols(meth.gr) = meth %>% select(-Chr_hg19, -Start_hg19, -End_hg19)

dn.ir = IRanges(start = dn_cases$Start, end = dn_cases$End)
dn.gr = GRanges(seqnames = dn_cases$Chrom, ranges = dn.ir)  #, strand = "+"
mcols(dn.gr) = dn_cases %>% select(-Chrom, -Start, -End)

# intersect dn with near TSS pos
# ranges = subsetByOverlaps(dn.gr, meth.gr)
ranges = subsetByOverlaps(meth.gr, dn.gr)
# transfer metadata for hits
# hits = findOverlaps(dn.gr, meth.gr)
hits = findOverlaps(meth.gr, dn.gr)
snp.id = CharacterList(split(dn.gr$snp.id[subjectHits(hits)], queryHits(hits)))
mcols(ranges) = DataFrame(mcols(ranges), snp.id)

# confirm max sample count
# as.data.frame(ranges) %>%
#   rowwise() %>%
#   mutate(sample_count = length(strsplit(SAMPLE, ",")[[1]]) ) %>%
#   arrange(desc(sample_count)) %>%
#   as.data.frame %>% head

# expand by individual, then see if de novo DNA_ID matches methylation ID
meth_hits = as.data.frame(ranges) %>%
  mutate(SAMPLE = gsub(".AVG_Beta", "", SAMPLE)) %>%
  mutate(SAMPLE = gsub("\\.", "-", SAMPLE)) %>%
  select(-strand) %>%
  rename(Chrom_meth = seqnames, Start_meth = start, End_meth = end) %>%
  separate(SAMPLE, paste("DNA_ID_", c(1:77), sep = ""), sep = ",") %>%
  gather(key = dna_id_key, value = DNA_ID_meth, -c(Chrom_meth:Number_samples), -c(Closest_Feature:snp.id)) %>%
  filter(!is.na(DNA_ID_meth)) %>%
  mutate(snp.id = as.character(snp.id)) %>% #Chrom_meth = gsub("^", "chr", Chrom_meth),
  left_join(meth_dna_id_blinded_id_map) %>%
  select(-dna_id_key)

# expand SNP ID hits
meth_hits %>% rowwise() %>%
    mutate(snp.id_count = length(strsplit(snp.id, ",")[[1]]) ) %>%
    arrange(desc(snp.id_count)) %>% as.data.frame %>% head

meth_hits_snp_expanded = meth_hits %>%
  mutate(snp.id = as.character(snp.id) %>%
           gsub("c\\(", "", .) %>%
           gsub("\\)", "", .) %>%
           gsub("\"", "", .)) %>%
  separate(snp.id, paste("snp.id.", c(1:28), sep = ""), sep = ", ") %>%
  gather(snp_id_count, snp.id, -c(Chrom_meth:Class), -c(DNA_ID_meth:Blinded.ID)) %>%
  filter(!is.na(snp.id)) %>%
  select(-snp_id_count)

# how many of the blinded IDs with methylation data also have WGS?
sum(unique(meth_hits_snp_expanded$Blinded.ID) %in% dn$Blinded.ID)
dn_meth = meth_hits_snp_expanded %>% inner_join(dn_cases)
dim(dn_meth)
# dn_meth[1:5, 1:30]
