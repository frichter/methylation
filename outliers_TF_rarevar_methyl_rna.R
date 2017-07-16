# Felix Richter
# 7/2/2017
# Check expression of all genes with rare variants in TFBS and methylation outliers
############################


setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

p = c("limma", "IRanges", "GenomicRanges", "edgeR", "variancePartition", "purrr",
      "org.Hs.eg.db", "annotate")
lapply(p, require, character.only = TRUE)
p = c("dplyr", "ggplot2", "tidyr", "stringi", "readr")
lapply(p, require, character.only = TRUE)

source("methylation/outliers_TF_rarevar_methyl_rna_functions.R")

#######################################
# load TF and methylation data 
#######################################

# 2017_06_23_rarevar_TFs
tf_meth_path = "methylation/dnm_meth_results/2017_07_05_rarevar_TFs/"
tf_meth_df = LoadTFmeth(tf_meth_path)
meth_ids = tf_meth_df$Samples %>% unique

sig_tfs = read_tsv("methylation/dnm_meth_results/2017_07_metadata/significant_TFs.txt")$Factor
## hyper TFs: CTCF, PU1, c("CTCF", "PU1", "SP1")
## hypo TFs: AP1, CTCF, YY1 c("AP1", "CTCF", "YY1")

tf_meth_df = tf_meth_df %>%
  filter(TF %in% sig_tfs)

# #################################################
# in how many ids was each rare variant observed
# #################################################

p = tf_meth_df %>% select(Samples, Rare.variant.MapInfo) %>% unique %>% 
  group_by(Rare.variant.MapInfo) %>% tally %>% ungroup %>% arrange(desc(n)) %>% 
  filter(!grepl(";", Rare.variant.MapInfo)) %>%
  # filter(n > 2) %>%
  # as.data.frame %>% 
  ggplot(., aes(x = n)) + 
  geom_histogram(bins = 9) +
  # xlim(0,12) + 
  theme_classic()

p

# ggsave("methylation/dnm_meth_results/2017_07_05_rare_vars_per_id.png", width = 3, height = 3)

### how many probes per TFBS with a rare variant?

meth_rank_summary = tf_meth_df %>% 
  mutate(rank_category = ifelse(Rank < mean(Rank), "hypOmeth/low", "hypERmeth/high")) %>% 
  group_by(Samples, TF, Rare.variant.MapInfo) %>%  # gene, Blinded.ID
  arrange(rank_category) %>% 
  summarise(rank_categories = paste(unique(rank_category), collapse = ", "), probes_per_locus = n()) %>% 
  ungroup()

p = meth_rank_summary %>% 
  ggplot(., aes(x = probes_per_locus)) +
  geom_histogram(bins = 20) +
  facet_wrap(~rank_categories) + 
  theme_classic()
p
# ggsave("methylation/dnm_meth_results/2017_07_metadata/probes_per_id_per_TFBS_per_variant.png", 
#        width = 7, height = 2.5)

#######################################
# create a location dataframe for genes
# and annotate with gene lists
#######################################

# use same name and position annotations as RNA-Seq
source("whole_genome/current_pipeline/baylor_pipeline_funcs.R")
fc_loc = "rna_seq_other/expression_data/counts.gene.RDS"
gene_info_df = ExtractGeneDFfromFeatureCountsRefSeq(fc_loc) 
gene_df = CreateGenePosDf(gene_info_df, 10^6)

# annotate with causative gene lists
heart_lists = ImportHeartLists()

# return df of TRUE and FALSE corresponding to if a gene is in an interesting gene list
gene_df_anno = sapply(heart_lists, function(hl_i) gene_df$gene %in% hl_i) %>%
  data.frame %>%
  bind_cols(gene_df, .)

# remove genes with duplicate TSS, keep the ones that are heart relevant or longer
gene_df_anno = RemoveGenesWithDupTSS(gene_df_anno)

gene_df_anno$mm_Cardio_KO %>% sum

#################################################################
# intersection between gene TSS +/- 10^6 and methylation outliers
#################################################################

gene.meth.df.expanded = IntersectTFBSandGeneTSS(gene_df_anno, tf_meth_df)

## finding closest gene to each TFBS
gene.meth.df.closest = gene.meth.df.expanded %>%
  group_by(Chrom, Start, End) %>%
  filter(abs(tss.dist) == min(abs(tss.dist))) %>% 
  # filter(tss.dist <= 0) %>% # | abs(tss.dist) < 5000
  ungroup() %>%
  unique

# gene.meth.df.closest %>% group_by(TF) %>% tally

gene.meth.df.closest %>% 
  filter(gene == "SMAD5")

#####################
# gene set enrichment
#####################

## Need correct background list! All genes with a TSS within 1kb of TFBS, where the TFBS
## also has a probe within 1kb

# k_all_genes = gene_df$gene %>% unique
# gene.meth.df.closest %>% head %>% as.data.frame
# 
# ## calculate median rank per gene (should this be calculated separately for separate TFs?)
# med_rank_per_gene = gene.meth.df.closest %>% 
#   filter(abs(tss.dist) < 10^3) %>% 
#   rename(Blinded.ID = Samples) %>% 
#   mutate(rank_category = ifelse(Rank < mean(Rank), "hypOmeth/low", "hypERmeth/high")) %>%
#   group_by(Blinded.ID, gene, 
#            hhe, hbe, omim_sarah, cilia, wes_multihit, chea_TFs, mm_Cardio_KO) %>%  
#   arrange(rank_category) %>% 
#   summarise(rank_categories = paste(unique(rank_category), collapse = ", "), 
#             median_rank = median(Rank), 
#             probes_per_locus = n()) %>% 
#   # filter(probes_per_locus > 1) %>% 
#   # filter(grepl(",", rank_categories))
#   ungroup %>% 
#   # ggplot(., aes(x = median_rank)) + geom_histogram()
#   mutate(median_rank_category = ifelse(median_rank < 50, "hypOmeth/low", "none")) %>% 
#   mutate(median_rank_category = ifelse(median_rank > 200, "hypERmeth/high", 
#                                        median_rank_category)) %>% 
#   gather(key = "gene_set", "gene_set_member", hhe:mm_Cardio_KO) 
# 
# ## dataframe of total number of genes in each set
# heart_list_ct_df = map(heart_lists, ~ sum(. %in% all_genes)) %>% stack
# names(heart_list_ct_df) = c("gene_set_total", "gene_set")
# heart_list_ct_df %<>% 
#   select(gene_set, gene_set_total) %>% 
#   mutate(gene_set = as.character(gene_set),
#          all_gene_n = length(k_all_genes)) %>% 
#   mutate(bg_fraction = gene_set_total/all_gene_n)
# 
# ### calculate gene set enrichments
# med_rank_per_gene %>% 
#   group_by(median_rank_category, gene_set) %>% 
#   summarise(gene_set_member_obs = sum(gene_set_member), total_genes_obs = n()) %>% 
#   ungroup %>% 
#   left_join(heart_list_ct_df, by = "gene_set") %>% 
#   rowwise() %>% 
#   mutate(p.binom = binom.test(gene_set_member_obs, total_genes_obs, p = bg_fraction)$p.value,
#          or = (gene_set_member_obs/total_genes_obs)/bg_fraction) %>% 
#   as.data.frame
# # p.fet = fisher.test(cbind("Case" = c(Case, n_cases_muts - Case), 
# #                           "Ctrl" = c(Ctrl, n_ctrl_muts - Ctrl)))$p.value)

################################
# join with gene expression data
################################

# 0.01*247
# 0.025 tails instead of 0.05

## use library and tissue or PEER corrected?
# r_long_file = "rna_seq_other/expression_data/lib_tissue_plat/residuals_z_avg.RDS"
# r_long_file = "rna_seq_other/expression_data/eqtl_prep/residual_expr_long.RDS"
# r_long_file = "whole_genome/rare_variants_eqtl/r_long_full_cohort_pcgc.RDS"
r_long_file = "whole_genome/rare_variants_eqtl/r_long_full_cohort_z_pcgc.RDS"
# r_long_file = "rna_seq_other/expression_data/R_long_atr_vent_art_library_plat_sep_calc.RDS"
R_long = readRDS(r_long_file)

## calculate z expression (if not already done)
head(R_long)
# R_mean_sd = R_long %>% group_by(gene) %>%
#   summarise(mean_expr = mean(expr), sd_expr = sd(expr)) %>% ungroup
# R_long = R_long %>%
#   left_join(R_mean_sd) %>%
#   mutate(z_expr = (expr - mean_expr)/sd_expr)

## confirm these are the same length:
R_long$Sample.ID %>% unique %>% length
R_long$Blinded.ID %>% unique %>% length

## if not then remove the extra Samples:
R_long = R_long %>% filter(Blinded.ID %in% meth_ids) %>% 
  ## For Blinded IDs with multiple samples from different tissues, keep 1 (most relevant) sample
  ## 1-01984 is HLHS so keeping Ventricle
  filter(!(Sample.ID %in% c("1-01984--LA1", "1-04056--DuctArt1")))

# #### calculating ranks
# R_long = R_long %>% 
#   group_by(gene) %>% 
#   mutate(expr_rank = rank(-z_expr)) %>% 
#   ungroup 
# 
### join RNA and methylation data

rna_ids = R_long$Blinded.ID %>% unique
gene_meth_tf_expr_df = gene.meth.df.closest %>%
  mutate(Blinded.ID = Samples) %>%
  filter(Blinded.ID %in% rna_ids) %>%
  left_join(R_long)

genes_wo_expr = gene_meth_tf_expr_df %>% filter(is.na(expr)) %>% select(gene) %>% unique
gene_meth_tf_expr_df = gene_meth_tf_expr_df %>% filter(!is.na(expr))

gene_meth_tf_expr_df %>% 
  ggplot(., aes(x = z_expr, fill = mm_Cardio_KO)) +
  geom_histogram()

gene_meth_tf_expr_df %>% 
  # filter(tss.dist < 0) %>% 
  filter(abs(tss.dist) < 10^3) %>% 
  # filter(mm_Cardio_KO) %>%
  mutate(z_expr_abs = abs(z_expr)) %>% 
  summarise(mean_expr = mean(z_expr_abs), median_expr = median(z_expr_abs))

meth_genes = gene_meth_tf_expr_df$gene %>% unique

##############################################
# plot results. are there any obvious trends?
##############################################

R_long %>% 
  filter(gene %in% meth_genes) %>% 
  mutate(z_expr_abs = abs(z_expr)) %>% 
  summarise(mean_expr = mean(z_expr_abs), median_expr = median(z_expr_abs))

gene_meth_tf_expr_df %>% 
  filter(abs(tss.dist) < 10^3) %>% 
  ggplot(., aes(x = tss.dist, y = z_expr, col = mm_Cardio_KO)) +
  geom_point()

gene_meth_tf_expr_df %>% 
  filter(abs(tss.dist) < 10^3) %>%
  mutate(rank_category = ifelse(Rank < mean(Rank), "low", "high")) %>% 
  ggplot(., aes(x = rank_category, y = z_expr, col = mm_Cardio_KO)) +
  geom_boxplot()

################################
# calculate enrichments
################################
gene_meth_tf_expr_df$TF %>% unique %>% length

gene_meth_tf_expr_df %>% 
  # filter(mm_Cardio_KO) %>%
  # filter(Rank < 7 | Rank > 240) %>%
  filter(abs(tss.dist) < 1*10^3) %>% 
  # ### RANK-based:
  # mutate(rank_category = ifelse(Rank < mean(Rank), "hypOmeth/low", "hypERmeth/high"),
  #        outlier_category = ifelse(expr_rank >= 199, "negative expr outlier", "none")) %>%
  # mutate(outlier_category = ifelse(expr_rank <= 5, "positive expr outlier", outlier_category)) %>%
  #### z-score based:
  mutate(rank_category = ifelse(Rank < mean(Rank), "hypOmeth/low", "hypERmeth/high")) %>%
  mutate(outlier_category = ifelse(z_expr < -2, "negative expr outlier", "none")) %>%
  # mutate(outlier_category = ifelse(z_expr > 2, "positive expr outlier", "none")) %>%
  mutate(outlier_category = ifelse(abs(z_expr) > 2, "expr outlier", "none")) %>%
  group_by(Blinded.ID, outlier_category, gene) %>%  
  arrange(rank_category) %>% 
  summarise(rank_categories = paste(unique(rank_category), collapse = ", "), 
            median_rank = median(Rank), 
            probes_per_locus = n()) %>% 
  # filter(probes_per_locus > 1) %>% 
  # filter(grepl(",", rank_categories))
  ungroup %>% 
  # ggplot(., aes(x = median_rank)) + geom_histogram()
  mutate(median_rank_category = ifelse(median_rank < 50, "hypOmeth/low", "none")) %>% 
  mutate(median_rank_category = ifelse(median_rank > 200, "hypERmeth/high", median_rank_category)) %>% 
  # group_by(median_rank_category) %>% tally
  # expr_outlier_neg expr_outlier_status
  group_by(outlier_category, median_rank_category) %>% #  outlier_category rank_categories median_rank_category
  tally

R_long %>% 
  # filter(gene %in% heart_lists$mm_Cardio_KO) %>%
  filter(Blinded.ID %in% meth_ids) %>% 
  ## z-score based
  # mutate(outlier_category = ifelse(abs(z_expr) > 2, "outlier", "none")) %>%
  # mutate(outlier_category = ifelse(z_expr < -2, "negative expr outlier", "none")) %>%
  mutate(outlier_category = ifelse(z_expr > 2, "positive expr outlier", "none")) %>%
  ## rank based
  # # mutate(outlier_category = ifelse(expr_rank >= 203, "negative expr outlier", "none")) %>%
  # mutate(outlier_category = ifelse(expr_rank <= 2, "positive expr outlier", "none")) %>%
  # # mutate(outlier_category = ifelse(expr_rank <= 2 | expr_rank >= 67, "outlier", "none")) %>%
  group_by(outlier_category) %>%  ## expr_outlier_status, expr_outlier_neg
  tally


## all genes:
outs_near_meth_site = 18
non_outs_near_meth_site = 682


## using residuals from all 204, known covariates, and 15 PEER SVs
# expr_outs = 12688 # 21907 12688 9219
# expr_non_outs = 428176 # 418957 428176 431645

## using residuals from all 204, known covariates, and 15 PEER SVs (z-score ALSO from full)
expr_outs = 10756 # 14363 # 
expr_non_outs = 449276 # 445669 # 
## using residuals from each tissue after correcting for known covariates
# expr_outs = 9721
# expr_non_outs = 381385
## other
# expr_outs = 4222
# expr_non_outs = 384823

fet_mtx = as.matrix(cbind("outlier" = c(outs_near_meth_site, expr_outs-outs_near_meth_site), 
                          "non_outlier" = c(non_outs_near_meth_site, expr_non_outs - non_outs_near_meth_site)))
fisher.test(fet_mtx)

# outlier assocation with methylation level
hyper_outs = 28
hypo_outs = 25
hyper_non_outs = 1363
hypo_non_outs = 848
fet_mtx = as.matrix(cbind("outlier" = c(hyper_outs, hypo_outs), 
                          "non_outlier" = c(hyper_non_outs, hypo_non_outs)))
fisher.test(fet_mtx)


outs_fullz = 10756 # 10756 14363
non_outs_fullz = 449276 # 449276 445669
outs_tissueZ = 9721 # 14019
non_outs_tissueZ = 381385 # 377087
fet_mtx = as.matrix(cbind("fullZ" = c(outs_fullz, non_outs_fullz), 
                          "tissueZ" = c(outs_tissueZ, non_outs_tissueZ)))
fisher.test(fet_mtx) #$p.value


#################################################
# Continuous data comparison
#################################################

expr_rank_tf = gene_meth_tf_expr_df %>% 
  filter(mm_Cardio_KO | omim_sarah) %>%
  filter(TF == "CTCF") %>% 
  # filter(hhe) %>% 
  ## with TFBS within 1kb of TSS
  filter(abs(tss.dist) < 10^3) %>% 
  ## annotate for high/hypERmethylation outliers and gene expression outliers
  mutate(rank_category = ifelse(Rank < mean(Rank), "hypOmeth/low", "hypERmeth/high")) %>% 
  group_by(Blinded.ID, gene, z_expr, mm_Cardio_KO, omim_sarah, hhe) %>%  
  arrange(rank_category) %>% 
  summarise(rank_categories = paste(unique(rank_category), collapse = " + "), 
            median_rank = median(Rank), 
            probes_per_locus = n()) %>% 
  ungroup %>% 
  mutate(median_rank_category = ifelse(median_rank < 50, "hypOmeth/low", "none")) %>% 
  mutate(median_rank_category = ifelse(median_rank > 200, "hypERmeth/high", median_rank_category))

expr_rank_tf %>%
  # filter(probes_per_locus > 1) %>% 
  group_by(median_rank_category) %>% 
  summarise(z_mean = mean(z_expr), z_sd = sd(z_expr), n = n()) %>% 
  mutate(z_summary = (z_mean + 0.077869090)/(z_sd/sqrt(n)))
  # mutate(z_summary = z_mean/(z_sd/sqrt(n)))

expr_rank_tf %>% 
  ggplot(., aes(x = median_rank_category, y = z_expr, col = omim_sarah)) +
  geom_boxplot()

#################################################
# Plotting methylation + gene expression outliers
#################################################
tf_meth_df %>% filter(Rare.variant.MapInfo == 10971417, TF == "CTCF") %>% slice(1) %>% 
  as.data.frame %>% unlist %>% as.character %>% paste(collapse = ", ")

gene_meth_tf_expr_df %>% 
  # filter(TF == "CTCF") %>% 
  filter(abs(tss.dist) < 10^3) %>% 
  filter(gene == "CIITA") %>% as.data.frame
## label the expression and methylation outliers like you do above
meth_and_expr_outliers = gene_meth_tf_expr_df %>% 
  filter(chea_TFs | omim_sarah | mm_Cardio_KO) %>%
  filter(TF == "CTCF") %>% 
  # filter(hhe) %>% 
  ## with TFBS within 1kb of TSS
  filter(abs(tss.dist) < 10^3) %>% 
  ## annotate for high/hypERmethylation outliers and gene expression outliers
  mutate(rank_category = ifelse(Rank < mean(Rank), "hypOmeth/low", "hypERmeth/high")) %>%
  mutate(outlier_category = ifelse(z_expr < -2, "negative expr outlier", "none")) %>%
  mutate(outlier_category = ifelse(z_expr > 2, "positive expr outlier", outlier_category)) %>%
  group_by(Blinded.ID, outlier_category, gene) %>%  
  arrange(rank_category) %>% 
  summarise(rank_categories = paste(unique(rank_category), collapse = " + "), 
            median_rank = median(Rank), 
            probes_per_locus = n()) %>% 
  ungroup %>% 
  mutate(median_rank_category = ifelse(median_rank < 50, "hypOmeth/low", "none")) %>% 
  mutate(median_rank_category = ifelse(median_rank > 200, "hypERmeth/high", median_rank_category))

## only keep methylation adn expression outliers
outlier_genes = meth_and_expr_outliers %>% 
  filter(outlier_category != "none", median_rank_category != "none") %>%
  select(gene) %>% unlist %>% unique %>% as.character

r_long_file = "whole_genome/rare_variants_eqtl/r_long_full_cohort_z_pcgc.RDS"
# r_long_file = "rna_seq_other/expression_data/R_long_atr_vent_art_library_plat_sep_calc.RDS"
R_long = readRDS(r_long_file)
## in case using the old Blinded and Sample IDs..
R_long = R_long %>% 
  mutate(Sample.ID = ifelse(grepl("02520", Sample.ID), "1-06596--RA1", Sample.ID),
         Blinded.ID = ifelse(grepl("02520", Blinded.ID), "1-06596", Blinded.ID))

info = readRDS(file = "rna_seq_other/expression_data/info.RDS")
info = info %>% filter(Sample.ID %in% unique(R_long$Sample.ID))
tissue_df = info %>% select(Sample.ID, TissueCorrected)

## use either full cohort or a subset
Expr_anno = R_long %>% 
  # filter(Blinded.ID %in% meth_ids) %>% 
  filter(gene %in% outlier_genes) %>% 
  left_join(meth_and_expr_outliers) %>% 
  left_join(tissue_df)

hyper_outlier_genes = meth_and_expr_outliers %>% 
  filter(outlier_category != "none", median_rank_category == "hypERmeth/high") %>%
  select(gene) %>% unlist %>% unique %>% as.character

hypOmeth_genes = meth_and_expr_outliers %>% 
  filter(outlier_category != "none", median_rank_category == "hypOmeth/low") %>%
  select(gene) %>% unlist %>% unique %>% as.character

out_cat_order = c("none", "negative expr outlier", "positive expr outlier")
Expr_anno = Expr_anno %>% unique
p = Expr_anno %>% 
  filter(gene %in% outlier_genes) %>% 
  mutate(outlier_category = ifelse(is.na(outlier_category), "none", outlier_category)) %>%
  mutate(outlier_category = factor(outlier_category, levels = out_cat_order)) %>%
  # filter(median_rank_category == "hypOmeth/low" | is.na(median_rank_category)) %>%
  filter(gene %in% hypOmeth_genes) %>%
  # filter(median_rank_category == "hypERmeth/high" | is.na(median_rank_category)) %>%
  # filter(gene %in% hyper_outlier_genes) %>%
  ggplot(., aes(x = gene, y = z_expr)) +
  scale_colour_manual(values = c("grey60", "blue")) + # "red", 
  geom_boxplot(aes(col = outlier_category), position = "identity") +
  # facet_wrap( ~ all_tissues_per_id) +
  # facet_wrap( ~ rank_categories, ncol = 1) +
  theme_classic()
p

ggsave("methylation/rna_seq_meth_results/2017_07_05_mouseCHD_omim_genes/hypOmeth_CITED2.png", 
       p, width = 3, height = 4)

###################################################
# Printing out Methylation outliers for comparisons
###################################################

# write_tsv(gene_meth_tf_expr_df, "methylation/rna_seq_meth_results/all_TFs_expr_TissueZ.txt")
write_tsv(gene_meth_tf_expr_df, "methylation/rna_seq_meth_results/all_TFs_expr_totalZ.txt")


gene_meth_tf_expr_df %>% 
  filter(abs(tss.dist) < 1*10^3) %>% 
  filter(omim_sarah) %>% 
  #### z-score based:
  mutate(rank_category = ifelse(Rank < mean(Rank), "hypOmeth/low", "hypERmeth/high")) %>% 
  filter(z_expr < -2 | z_expr > 2) %>% 
  mutate(outlier_category = ifelse(z_expr < -2, "negative expr outlier", "none")) %>%
  mutate(outlier_category = ifelse(z_expr > 2, "positive expr outlier", outlier_category)) %>%
  # mutate(outlier_category = ifelse(abs(z_expr) > 2, "expr outlier", "none")) %>%
  group_by(Blinded.ID, outlier_category, gene, z_expr) %>%  
  arrange(rank_category) %>% 
  summarise(rank_categories = paste(unique(rank_category), collapse = ", "), probes_per_locus = n()) %>% 
  ungroup %>% select(-outlier_category) %>% as.data.frame

mean(gene_meth_tf_expr_df$Rank)
median(gene_meth_tf_expr_df$Rank)

gene_meth_tf_expr_df %>% 
  filter(gene == "SMC3", Blinded.ID == "1-00508") %>% 
  summarise(median_rank = median(Rank), mean_rank = mean(Rank))


## using residuals from all 150 and known covariates
# expr_outs = 15247 # 432442 15247 12343
# expr_non_outs = 444785 # 432442 444785 447689
## using residuals from only 68 IDs and PEER covariates
# expr_outs = 12669 # 21135 12669 8466
# expr_non_outs = 376376 # 367910 376376 380579


# ### use the GTEx definition of outliers:
# R_long = R_long %>% 
#   group_by(gene) %>%
#   mutate(gtex_expr_outlier = abs(z_expr) == max(abs(z_expr)) & abs(z_expr) > 2) %>%
#   mutate(gtex_expr_outlier_neg = gtex_expr_outlier & z_expr < 0) %>%
#   ungroup()
# 
# 
# ## option play: keep only those with Atrial/Ventricle
# rna_info = readRDS("rna_seq_other/expression_data/eqtl_prep/info.RDS")
# # rna_info = readRDS("rna_seq_other/expression_data/lib_tissue_plat/info.RDS")
# atrial_vent_ids = rna_info %>% 
#   filter(TissueCorrected %in% c("Atrial", "Ventricle")) %>% 
#   select(Blinded.ID) %>% unique %>% unlist %>% as.character
# # R_long = R_long %>% 
# #   filter(Blinded.ID %in% atrial_vent_ids)
# 