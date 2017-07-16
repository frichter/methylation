# Felix Richter
# felix.richter@icahn.mssm.edu
# 4/15/2017
# Purpose: perform gene set enrichment using either mouse CHD genes or all genes
#          as the background
################################################################################

setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

p = c("limma", "IRanges", "GenomicRanges", "edgeR", "variancePartition", "purrr",
      "org.Hs.eg.db", "annotate")
lapply(p, require, character.only = TRUE)
p = c("dplyr", "ggplot2", "tidyr", "stringi", "readr")
lapply(p, require, character.only = TRUE)

source("whole_genome/current_pipeline/baylor_pipeline_funcs.R")

heart_lists = ImportHeartLists()
heart_list_lens = map_int(heart_lists, length)
mouse_all = read.table("gene_lists.chd/mouse_all.txt")$V1
mouse_bg = length(mouse_all)
human_bg = 25000

rna_info = readRDS("rna_seq_other/expression_data/info.RDS")
rna_exp = readRDS("rna_seq_other/expression_data/lib_tissue_plat/residuals_z.RDS")

ctcf_dmrs = read_tsv("methylation/dmrs/alex_2017_04_03_dmrs_ctcf.txt")
ctcf_dmrs %>% filter(Blind.ID.DMR %in% rna_info$Blinded.ID) %>% as.data.frame
ctcf_dmrs

ctcf_dmrs = ctcf_dmrs %>% filter(Gene_refGene != ".") #!grepl(",", Gene_refGene),
rna_exp %>% filter(grepl("1-01984", Sample.ID), gene == "CELSR1")

ctcf_dmr_genes = ctcf_dmrs$Gene_refGene %>%
  strsplit(ctcf_dmr_genes, split = ",") %>% unlist %>% unique
ctcf_obs_count = length(ctcf_dmr_genes)

map_df(heart_lists, ~ sum(ctcf_dmr_genes %in% .)) %>%
  gather("heart_list", "heart_obs_count") %>%
  mutate(heart_count = heart_list_lens,
         ctcf_obs_count = ctcf_obs_count,
         bg_length = ifelse(heart_list != "mm_Cardio_KO", human_bg, mouse_bg)) %>%
  rowwise() %>%
  mutate(or = (heart_obs_count/heart_count)/(ctcf_obs_count/bg_length),
         binom_p = binom.test(heart_obs_count, heart_count, ctcf_obs_count/bg_length)$p.value) %>%
  ungroup()

hhe_mouse_chd = c(heart_lists$hhe, heart_lists$mm_Cardio_KO, heart_lists$chromatin) %>% unique %>% sort
hhe_mouse_chd_chromatin_all = heart_lists$hhe[heart_lists$hhe %in% heart_lists$mm_Cardio_KO &
                                                heart_lists$hhe %in% heart_lists$chromatin]
hhe_mouse_chd = heart_lists$hhe[heart_lists$hhe %in% heart_lists$mm_Cardio_KO]

ctcf_dmr_genes[ctcf_dmr_genes %in% heart_lists$chromatin]

heart_obs_count = 2
heart_count = length(hhe_mouse_chd)
bg_length = human_bg
(heart_obs_count/heart_count)/(ctcf_obs_count/bg_length)
binom.test(heart_obs_count, heart_count, ctcf_obs_count/bg_length)
