# Felix Richter
# felix.richter@icahn.mssm.edu
# 7/3/2017
# Purpose: script to compare the post-PEER RNAseq results when using only the
#          IDs with WGS (68) vs the full dataset (204) to calculate the SVs
#############################################################################



setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

p = c("limma", "IRanges", "GenomicRanges", "edgeR", "variancePartition", "purrr",
      "org.Hs.eg.db", "annotate")
lapply(p, require, character.only = TRUE)
p = c("dplyr", "ggplot2", "tidyr", "stringi", "readr")
lapply(p, require, character.only = TRUE)


## use library and tissue or PEER corrected?
r204_long_file = "rna_seq_other/expression_data/lib_tissue_plat/residuals_z_avg.RDS"
r68_long_file = "rna_seq_other/expression_data/eqtl_prep/residual_expr_long.RDS"
R204_long = readRDS(r204_long_file)
R68_long = readRDS(r68_long_file)

wgs_rna_ids = R68_long$Blinded.ID %>% unique %>% as.character

CalculateZexpr = function(R_long) {
  R_mean_sd = R_long %>% group_by(gene) %>%
    summarise(mean_expr = mean(expr), sd_expr = sd(expr)) %>% ungroup
  R_long = R_long %>%
    left_join(R_mean_sd) %>%
    mutate(z_expr = (expr - mean_expr)/sd_expr)
  return(R_long)
}

R68_long = CalculateZexpr(R68_long)

R204_long = R204_long %>% filter(Blinded.ID %in% wgs_rna_ids)

R68_long = R68_long %>% select(gene, Blinded.ID, z_expr) %>% rename(z_expr_r68 = z_expr)
R204_long = R204_long %>% select(gene, Blinded.ID, z_expr) %>% rename(z_expr_r204 = z_expr)

R_long = full_join(R68_long, R204_long, by = c("gene", "Blinded.ID"))

# R_long %>% unique %>% dim

R_long %>%
  mutate(r68_na = !is.na(z_expr_r68), r204_na = !is.na(z_expr_r204)) %>%
  group_by(r68_na, r204_na) %>% tally

p = R_long %>%
  filter(!is.na(z_expr_r68), !is.na(z_expr_r204)) %>%
  # filter(z_expr_r68 < -2 & z_expr_r204 < -2) %>%
  filter(z_expr_r68 > 2 & z_expr_r204 > 2) %>%
  ggplot(., aes(x = z_expr_r204, y = z_expr_r68)) +
  geom_point(size = 0.3) +
  theme_classic()

p

R_long %>%
  filter(!is.na(z_expr_r68), !is.na(z_expr_r204)) %>%
  # filter(z_expr_r68 < -2 & z_expr_r204 < -2) %>%
  filter(z_expr_r68 > 2 & z_expr_r204 > 2) %>%
  summarise(cor(z_expr_r68, z_expr_r204))

## outliers in both or outliers in

z_cut_off = -2

R_long %>%
  mutate(r68_na = !is.na(z_expr_r68), r204_na = !is.na(z_expr_r204),
         r68_out = z_expr_r68 < z_cut_off, r204_out = z_expr_r204 < z_cut_off) %>%
  group_by(r68_na, r204_na, r68_out, r204_out) %>% tally %>%
  ungroup %>% mutate(pct = 100 * (n/1309544) %>% signif(1))

R_long %>%
  mutate(r68_na = !is.na(z_expr_r68), r204_na = !is.na(z_expr_r204)) %>%
  filter(!r68_na) %>% select(gene) %>% unique %>%
  dim

## 90 genes not in r204, 2343 not in r68
## use r204 because get increased power to detect decreased expression, especially
## negative outliers
## Perform PEER corrections for 204 IDs

## 2% outliers unique to each, 0.8% common to both

## Is there a bias in the IDs or genes that have outliers in either dataset?
## Yes but this seems normally distributed..

z_cut_off = -2

p = R_long %>%
  mutate(r68_no_na = !is.na(z_expr_r68), r204_no_na = !is.na(z_expr_r204),
         r68_out = z_expr_r68 < z_cut_off, r204_out = z_expr_r204 < z_cut_off) %>%
  filter(r68_no_na, r204_no_na) %>%
  filter(r68_out | r204_out) %>%
  group_by(Blinded.ID) %>% # Blinded.ID gene
  summarise(n = n(), r68_out_pct = (sum(r68_out & !r204_out)/n) %>% signif(4),
            r204_out_pct = (sum(!r68_out & r204_out)/n) %>% signif(4),
            both_out_pct = (sum(r68_out & r204_out)/n) %>% signif(4)) %>%
  ungroup %>%
  ggplot(., aes(x = both_out_pct, y = r68_out_pct)) + # r68_out_pct r204_out_pct
  geom_jitter() +
  theme_classic()

p

#############################################################################
# # eqtl_outliers
# full_dataset_outliers = gene_meth_tf_expr_df %>%
#   # filter(mm_Cardio_KO) %>%
#   filter(abs(tss.dist) < 1*10^3) %>%
#   # ### RANK-based:
#   # mutate(rank_category = ifelse(Rank < mean(Rank), "hypOmeth/low", "hypERmeth/high"),
#   #        outlier_category = ifelse(expr_rank >= 199, "negative expr outlier", "none")) %>%
#   # mutate(outlier_category = ifelse(expr_rank <= 5, "positive expr outlier", outlier_category)) %>%
#   #### z-score based:
#   mutate(rank_category = ifelse(Rank < mean(Rank), "hypOmeth/low", "hypERmeth/high"),
#          outlier_category = ifelse(z_expr < -2, "negative expr outlier", "none")) %>%
#   mutate(outlier_category = ifelse(z_expr > 2, "positive expr outlier", outlier_category)) %>%
#   # select(Blinded.ID, outlier_category, gene, rank_category) %>% unique %>% # rank_category
#   group_by(Blinded.ID, outlier_category, gene, z_expr) %>%
#   arrange(rank_category) %>%
#   summarise(rank_categories = paste(unique(rank_category), collapse = ", "), probes_per_locus = n()) %>%
#   # filter(probes_per_locus > 1) %>%
#   # filter(grepl(",", rank_categories))
#   ungroup %>%
#   filter(outlier_category != "none")
#
# wgs_id_expr = R_long %>%
#   filter(gene %in% eqtl_outliers$gene | gene %in% full_dataset_outliers$gene)
# wgs_id_expr = wgs_id_expr %>%
#   select(gene, Blinded.ID, z_expr) %>%
#   mutate(cohort = "wgs_ids")
#
#
# ##########
# full_data_expr = R_long %>%
#   filter(gene %in% eqtl_outliers$gene | gene %in% full_dataset_outliers$gene)
# full_data_expr = full_data_expr %>%
#   select(gene, Blinded.ID, z_expr) %>%
#   mutate(cohort = "all_ids")
#
# ### comparisons
# eqtl_outliers %>% dim
# full_dataset_outliers %>% dim
#
# outlier_df_list = list(eqtl_outliers, full_dataset_outliers)
# names(outlier_df_list) = c("wgs_ids", "all_ids")
#
# all_outliers = bind_rows(outlier_df_list, .id = "cohort")
#
# head(eqtl_outliers)
# head(full_dataset_outliers)
#
# eqtl_outlier_ids = eqtl_outliers %>% select(Blinded.ID) %>% unlist %>% as.character
# full_dataset_outlier_ids = full_dataset_outliers %>% select(Blinded.ID) %>% unlist %>% as.character
#
# wgs_outlier_rows_to_keep = eqtl_outliers %>% select(Blinded.ID, gene)
# full_dataset_outlier_rows_to_keep = full_dataset_outliers %>% select(Blinded.ID, gene)
#
# full_data_expr %>%
#   inner_join(wgs_outlier_rows_to_keep)
#
# all_outliers %>%
#   filter(cohort == "all_ids") %>%
#   select(gene, Blinded.ID, z_expr, cohort) %>%
#   # bind_rows(full_data_expr) %>%
#   # inner_join(wgs_outlier_rows_to_keep) %>%
#   bind_rows(wgs_id_expr) %>%
#   inner_join(full_dataset_outlier_rows_to_keep) %>%
#   # filter(rank_categories == "hypERmeth/high") %>%
#   # filter(rank_categories == "hypOmeth/low") %>%
#   mutate(shared_genes = gene %in% shared_genes) %>%
#   filter(shared_genes) %>%
#   ggplot(., aes(x = gene, y = z_expr, col = cohort)) +
#   geom_boxplot(position = "identity") +
#   # facet_wrap( ~ rank_categories, ncol = 1) +
#   theme_classic()
#
# shared_genes = eqtl_outliers %>%
#   filter(gene %in% full_dataset_outliers$gene) %>%
#   select(gene) %>% unlist %>% as.character
#
