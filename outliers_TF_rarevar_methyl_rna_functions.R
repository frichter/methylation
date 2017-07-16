
LoadTFmeth = function(tf_meth_path) {
  # tf_names = gsub(".*2017_07_05_rarevar_TFs/", "", tf_meth_files) %>% 
  #   gsub(".rare.variants.1000.*", "", .)
  tf_meth_files = list.files(tf_meth_path, full.names = T)
  tf_names = gsub(tf_meth_path, "", tf_meth_files) %>% 
    gsub("^/", "", .) %>% 
    gsub(".rare.variants.1000.*", "", .)
  names(tf_meth_files) = tf_names
  tf_meth_df = map_df(tf_meth_files, read_tsv, .id = "TF", col_types = c("ciiiicccdic"))
  tf_meth_df = tf_meth_df %>% 
    mutate(Chr.hg19 = paste0("chr", Chr.hg19))
  return(tf_meth_df)
}


IntersectTFBSandGeneTSS = function(gene_df_anno, tf_meth_df) {
  # intersection between gene TSS +/- 10^6 and methylation outliers
  # all data should be 1-based 
  gene.ir = IRanges(start = gene_df_anno$Start, end = gene_df_anno$End)
  gene.gr = GRanges(seqnames = gene_df_anno$Chr, ranges = gene.ir)
  mcols(gene.gr) = gene_df_anno %>% select(-Chr, -Start, -End)
  
  ## switching from methylation (Start.hg19.1000.bp/End.hg19.1000.bp) 
  ## to TFBS site
  tf_meth_df.ir = IRanges(start = tf_meth_df$Start.TFBS, end = tf_meth_df$End.TFBS)
  tf_meth_df.gr = GRanges(seqnames = tf_meth_df$Chr.hg19, ranges = tf_meth_df.ir)  #, strand = "+"
  # mcols(tf_meth_df.gr) = tf_meth_df %>% select(-Chr.hg19, -Start.hg19.1000.bp, -End.hg19.1000.bp)
  mcols(tf_meth_df.gr) = tf_meth_df %>% 
    select(Samples, Meth, Rank, Start.hg19.1000.bp, End.hg19.1000.bp, TF, Rare.variant.MapInfo)
  
  # intersect dn with near TSS pos
  ranges = subsetByOverlaps(tf_meth_df.gr, gene.gr)
  # ranges = subsetByOverlaps(gene.gr, tf_meth_df.gr)
  # transfer metadata for hits
  hits = findOverlaps(tf_meth_df.gr, gene.gr)
  # hits = findOverlaps(gene.gr, tf_meth_df.gr)
  gene = CharacterList(split(gene.gr$gene[subjectHits(hits)], queryHits(hits)))
  mcols(ranges) = DataFrame(mcols(ranges), gene)
  
  # gene.dn.df = as.data.frame(ranges) %>%
  #   mutate(gene = as.character(gene))
  
  # account for TFBS regions being near multiple genes
  gene.meth.df.expanded = as.data.frame(ranges) %>% 
    mutate(gene = as.character(gene) %>%
             gsub("c\\(", "", .) %>%
             gsub("\\)", "", .) %>%
             gsub("\"", "", .) %>% 
             gsub("\n", "", .)) %>% 
    # mutate(gene = gsub("c\\(", " ", gene)) %>% head
    separate(gene, paste("g", c(1:600), sep = ""), sep = ", ", fill = "right") %>%
    gather(key = gene_count, value = gene, c(g1:g600)) %>%
    filter(!is.na(gene))
  
  # find distance btw TFBS/DMR and TSS
  gene.meth.df.expanded = gene.meth.df.expanded %>%
    left_join(gene_df_anno, by = c("gene" = "gene")) %>%
    dplyr::rename(Start_gene = Start, End_gene = End, Chrom_gene = Chr, Strand_gene = Strand) %>%
    dplyr::rename(Start = start, End = end, Chrom = seqnames) %>%
    select(-strand, -width) %>%
    # negative distance if upstream of TSS, positive distance if downstream
    mutate(tss.dist = ifelse(Strand_gene == "+", Start - tss.pos, tss.pos - Start))
  return(gene.meth.df.expanded)
}
