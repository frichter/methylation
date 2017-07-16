# Felix Richter
# Look at expression associated with epimutations from Mafalda's 450K analysis
# 7/16/2017
################################################################################

## load dictionary for blinded ID to DNA methylation sample ID
blinded_to_dna_id = read_csv("methylation/metadata/dmr_id_blinded_id_mafalda.csv")
blinded_to_dna_id = blinded_to_dna_id %>%
  rename(Blinded.ID = `Blinded ID`,
         Sample_ID = `DNA Sample ID`) %>%
  select(Blinded.ID, Sample_ID)

## load DMRs
dmrs = read_tsv("methylation/dnm_meth_results/Felix_07_14_2017_DMRs.txt")
dmrs = dmrs %>% left_join(blinded_to_dna_id)
dmr_blinded_ids = dmrs$Blinded.ID %>% unique


# R_per_sample_loc = "rna_seq_other/expression_data/non_syndromic/lib_tissue_plat_gender/residuals.gene.RDS"
# R_per_sample = readRDS(R_per_sample_loc)

## non-syndromic
R_long_loc = "rna_seq_other/expression_data/non_syndromic/lib_tissue_plat_gender/r_long.gene.RDS"
info_loc = "rna_seq_other/expression_data/non_syndromic/info.07_2017.RDS"

R_long = readRDS(R_long_loc)
info = readRDS(info_loc)

rna_seq_ids = R_long$Blinded.ID %>% unique

sum(rna_seq_ids %in% dmr_blinded_ids)

R_long_dmr_ids = R_long %>%
  filter(Blinded.ID %in% dmr_blinded_ids) %>%
  rename(RNA.Sample.ID = Sample.ID) %>%
  select(Blinded.ID, gene, z_expr, expr_rank, TissueCorrected)

## only Ventricle and Atrial
R_long_dmr_ids$TissueCorrected %>% unique %>% paste(collapse = ", ")

dmrs_long = dmrs %>%
  select(Blinded.ID, Coordinates_hg19, Gene_Promoter_4kbTSS:Closest_Feature) %>%
  filter(Blinded.ID %in% rna_seq_ids) %>%
  gather(key = "gene_association", "gene", Gene_Promoter_4kbTSS:Closest_Feature) %>%
  filter(!is.na(gene))

dmrs_expr = dmrs_long %>%
  left_join(R_long_dmr_ids, by = c("Blinded.ID", "gene" = "gene")) %>%
  filter(!is.na(z_expr)) ## 4 genes without expression (2 lncrna, 1 AS1, one gene-doublet)

dmrs_expr %>%
  group_by(gene_association, TissueCorrected) %>%
  summarise(z_abs_mean = mean(abs(z_expr)), z_abs_median = median(abs(z_expr)), n) %>%
  ungroup

dmrs_expr %>%
  group_by(gene_association) %>%
  summarise(names = paste(collapse = ", ", gene), n()) %>% as.data.frame
## 5 individuals
##

all_gene_expr = R_long %>%
  ## remove the blinded IDs
  filter(!(Blinded.ID %in% dmr_blinded_ids)) %>%
  ## only keep the genes with DMRs
  filter(gene %in% unique(dmrs_expr$gene)) %>%
  filter(TissueCorrected %in% c("Ventricle", "Atrial")) %>%
  mutate(Coordinates_hg19 = "NA", gene_association = "Background") %>%
  select(Blinded.ID, Coordinates_hg19, gene_association, gene, z_expr, expr_rank, TissueCorrected) %>%
  bind_rows(dmrs_expr) %>%
  mutate(z_expr_absolute = abs(z_expr))

dmrs_expr %>% filter(gene == "ELANE")
all_gene_expr$gene %>% unique %>% length
dmrs_expr$gene %>% unique %>% length

p = all_gene_expr %>%
  filter(TissueCorrected == "Ventricle") %>%
  ggplot(., aes(x = z_expr_absolute, fill = gene_association)) +
  geom_density(alpha = 0.5) +
  theme_classic()

p
ggsave("methylation/rna_seq_meth_results/population_dmrs_facet.png", p, width = 6, height = 4)

## list of genes ordered by expression (to make the plots prettier)
prom_genes = dmrs_expr %>% filter(grepl("Promoter", gene_association)) %>%
  filter(TissueCorrected == "Ventricle") %>%
  arrange(z_expr) %>%
  select(gene) %>%
  unique %>% unlist %>% as.character
body_genes = dmrs_expr %>% filter(grepl("Body", gene_association)) %>%
  filter(TissueCorrected == "Ventricle") %>%
  arrange(z_expr) %>% select(gene) %>%
  unique %>% unlist %>% as.character
closest_genes = dmrs_expr %>% filter(grepl("Closest", gene_association)) %>%
  filter(TissueCorrected == "Ventricle") %>%
  arrange(z_expr) %>% select(gene) %>% unique %>% unlist %>% as.character

## bind each gene name to the blinded ID, so that you have a separate
## column for multiple blinded IDs with the same gene and similar expression
## To answer the question: How do I graph multiple HCP5 IDs?
prom_geneID_dict = dmrs_expr %>% filter(grepl("Promoter", gene_association)) %>%
  filter(TissueCorrected == "Ventricle") %>%
  ## pick the order in which teh genes are sorted:
  arrange(desc(z_expr)) %>% ## Blinded.ID,
  mutate(gene_blinded_id = paste(gene, Blinded.ID, sep = "\n")) %>%
  select(gene, gene_blinded_id)

## plot expression by gene, excel table of numbers
p = all_gene_expr %>%
  filter(TissueCorrected == "Ventricle") %>%
  ## Promoter Body Closest
  filter(grepl("Promoter|Background", gene_association)) %>%
  ## prom_genes body_genes closest_genes
  filter(gene %in% prom_genes) %>%
  full_join(prom_geneID_dict) %>%
  # filter(gene_association == "Gene_Promoter_4kbTSS", gene == "HCP5") %>% as.data.frame
  # group_by(gene_blinded_id, gene_association) %>% tally
  mutate(gene_blinded_id = factor(gene_blinded_id, levels = prom_geneID_dict$gene_blinded_id)) %>%
  ggplot(., aes(x = gene_blinded_id, y = z_expr, col = gene_association)) +
  geom_boxplot(size = 0.5, position = "identity") + ##
  # scale_colour_manual(values = c("grey60", "darkgreen", "blue", "red")) + #
  scale_colour_manual(values = c("black", "red")) + #
  theme_classic()

p
## dmr_promoter_gene_expr.png dmr_body_gene_expr.png dmr_closest_gene_expr.png
ggsave("methylation/case_ctrl_dmrs/dmr_promoter_gene_expr_2.png", p, width = 8, height = 5)


dmrs_expr %>% filter(grepl("Promoter", gene_association)) %>%
  filter(TissueCorrected == "Ventricle") %>%
  mutate(z_abs = z_expr %>% abs %>% round(digits = 2)) %>% as.data.frame

R_long %>%
  ## remove the blinded IDs
  filter(!(Blinded.ID %in% dmr_blinded_ids)) %>%
  mutate(z_round = z_expr %>% round(digits = 2)) %>%
  mutate(z_abs = z_round %>% abs) %>%
  mutate(outlier = z_round < -2) %>%
  group_by(outlier) %>% tally

## all outliers
binom.test(5, 9, 162768/3449967)
fisher.test(cbind("obs" = c(5, 4), "exp" = c(162768, 3449967-162768)))
## negative outliers
binom.test(4, 9, 100797/3511938)


########## Preparing an excel table of expression values

dim(dmrs)

R_long_dmr_ids %>% head
R_long_dmr_ids_prom = R_long_dmr_ids %>% select(-TissueCorrected) %>%
  rename(Gene_Promoter_4kbTSS = gene,
         Gene_Promoter_4kbTSS_z_expr = z_expr,
         Gene_Promoter_4kbTSS_expr_rank = expr_rank)
R_long_dmr_ids_body = R_long_dmr_ids %>% select(-TissueCorrected) %>%
  rename(Gene_Body = gene,
         Gene_Body_z_expr = z_expr,
         Gene_Body_expr_rank = expr_rank)
R_long_dmr_ids_closest = R_long_dmr_ids %>%
  rename(Closest_Feature = gene,
         Closest_Feature_z_expr = z_expr,
         Closest_Feature_expr_rank = expr_rank)

dmrs_expr = dmrs %>%
  left_join(R_long_dmr_ids_prom) %>%
  left_join(R_long_dmr_ids_body) %>%
  left_join(R_long_dmr_ids_closest)


dmrs_expr %>% select(Gene_Promoter_4kbTSS_z_expr:Gene_Promoter_4kbTSS_expr_rank) %>% as.data.frame

dmrs_expr %>%
  select(Coordinates_hg19, Blinded.ID, Gene_Promoter_4kbTSS_z_expr:TissueCorrected) %>%
  mutate_at(vars(Gene_Promoter_4kbTSS_z_expr:Closest_Feature_expr_rank), function(x) round(x, digits = 2)) %>%
  write_tsv(., "methylation/dnm_meth_results/Felix_07_14_2017_DMRs_expr.txt")
