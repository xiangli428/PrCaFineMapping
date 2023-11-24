options(stringsAsFactors = F, check.names = F)

library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggvenn)
library(ggpubr)
library(Hmisc)

setwd("PrCaFineMapping/gene_pathway_enrichment")

variant_info_CS95 = read.csv("../results/variant_info_CS95.csv", header = T,
                             stringsAsFactors = F, check.names = F)
variant_info_CS95$trans_genes_TCGA = ""

# CS95 genes
genes = list()
genes[["dbSNP"]] = foreach(j = 31, .combine = "union") %do%
{
  foreach(gs = unique(variant_info_CS95[,j]), .combine = "union") %do%
  {
    strsplit(gs, ",")[[1]]
  }
}

genes[["eQTL"]] = foreach(j = 32:34, .combine = "union") %do%
{
  foreach(gs = unique(variant_info_CS95[,j]), .combine = "union") %do%
  {
    strsplit(gs, ",")[[1]]
  }
}

genes[["Hi-C"]] = foreach(j = 35:37, .combine = "union") %do%
{
  foreach(gs = unique(variant_info_CS95[,j]), .combine = "union") %do%
  {
    strsplit(gs, ",")[[1]]
  }
}

genes[["H3K27ac HiChIP"]] = foreach(j = 38, .combine = "union") %do%
{
foreach(gs = unique(variant_info_CS95[,j]), .combine = "union") %do%
  {
    strsplit(gs, ",")[[1]]
  }
}

genes_union = foreach(j = names(genes), .combine = "union") %do%
{
  genes[[j]]
}

write.table(genes_union, "CS95_genes.txt", sep = '\t', row.names = F,
            col.names = F, quote = F)

p_venn = ggvenn(genes[c(2,1,4,3)], set_name_size = 4, text_size = 2.5)

# Enrichment analsysis
dbs = c("GO_BP", "WikiPathways")
gene_enrichment = foreach(db = dbs, .combine = "rbind") %dopar%
{
  tmp = read.delim(sprintf("enrich-input1-%s.tsv", db))
  tmp$Database = db
  tmp
}

gene_enrichment = gene_enrichment[gene_enrichment$term_genes >= 5,]
gene_enrichment$OR = gene_enrichment$genes_found *
  (gene_enrichment$universe - gene_enrichment$term_genes -
     gene_enrichment$input_size + gene_enrichment$genes_found) /
  (gene_enrichment$input_size - gene_enrichment$genes_found) /
  (gene_enrichment$term_genes - gene_enrichment$genes_found)
gene_enrichment = gene_enrichment[,c(12,1:9,13,11)]
gene_enrichment = arrange(gene_enrichment, pval)

# Pathway clustering and remove reduntant pathways
n = nrow(gene_enrichment)
gs = strsplit(gene_enrichment$genes, ", ")
simi_matrix = matrix(1, n, n)
for(i in 1:(n-1))
{
  li = length(gs[[i]])
  for(j in (i+1):n)
  {
    simi_matrix[i,j] = simi_matrix[j,i] = 
      2 * length(intersect(gs[[i]], gs[[j]])) / (li + length(gs[[j]]))
  }
}

keep_annos = c(1)
for(j in 2:n)
{
  if(all(simi_matrix[j, keep_annos] <= 0.3))
  {
    keep_annos = c(keep_annos,j)
  }
}

gene_enrichment_sub = gene_enrichment[keep_annos,]
simi_matrix_sub = simi_matrix[keep_annos,keep_annos]

gene_enrichment_sub$pval_adj = p.adjust(gene_enrichment_sub$pval, method = "BH")
write.table(gene_enrichment_sub, "gene_pathway_enrichment.txt", 
            sep = '\t', row.names = F, col.names = T, quote = F)

gene_enrichment_sig = gene_enrichment_sub[gene_enrichment_sub$pval_adj < 5e-3,]
gene_enrichment_sig$description = capitalize(gene_enrichment_sig$description)
gene_enrichment_sig$description = factor(
  gene_enrichment_sig$description, levels = rev(gene_enrichment_sig$description))

p_pathway = ggplot(gene_enrichment_sig, 
                   aes(x = -log10(pval_adj), y = description)) +
  geom_bar(stat = "identity", aes(color = log(OR), fill = log(OR))) +
  theme_classic() + theme(aspect.ratio = 0.5) +
  scale_color_gradient(low = "white", high = "red", limits = c(0,4)) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0,4)) +
  labs(x = TeX("$-log_{10}(P_{adj})$"), y = "", color = "log(OR)", fill = "log(OR)")

ggsave("gene_pathway_enrichment.pdf",
       ggarrange(p_venn, p_pathway,
                 nrow = 2, ncol = 1,
                 labels = c("A","B"),
                 heights = c(1,1),
                 font.label = list(size = 40/3, color = "black", face = "bold", 
                                   family = NULL)),
       device = "pdf", width = 8, height = 6, units = "in")
