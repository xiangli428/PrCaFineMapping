options(stringsAsFactors = F, check.names = F)

library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)
library(ggplot2)
library(scales)
library(Hmisc)

setwd("~/Documents/github/PrCaFineMapping")

registerDoParallel(22)

#1. Enrichment of cistrome annotation
variant_info = foreach(i = 1:22, .combine = "rbind") %dopar%
{
  read.delim(gzfile(sprintf(
    "results/variant_info_chr%d.txt.gz", i)), 
    stringsAsFactors = F, check.names = F)
}
variant_info$trans_genes_TCGA[is.na(variant_info$trans_genes_TCGA)] = ""

cistrome_meta = read.delim("cistrome_meta.txt")

variant_cistrome = separate_longer_delim(variant_info[,c(1,18,28)], features, ',')
variant_cistrome$value = 1
variant_cistrome = spread(variant_cistrome, features, value, fill = 0)[,-3]

m = sum(!is.na(variant_info$CS95)) #the number of CS95 variants
n = nrow(variant_info) - m #the number of non-CS95 variants
cistrome_meta[,5:6] = foreach(id = cistrome_meta$Factor, .combine = "rbind") %dopar%
{
  n11 = sum(!is.na(variant_cistrome$CS95) & (variant_cistrome[[id]] == 1))
  n1 = sum(variant_cistrome[[id]])
  n01 = n1 - n11
  n10 = m - n11
  n00 = n - m - n01
  
  data.frame(OR = exp(log(n11) + log(n00) - log(n10) - log(n01)),
             phyper = phyper(n11 - 1, m, n, n1, lower.tail = F))
}
cistrome_meta$phyper.adjust = p.adjust(cistrome_meta$phyper)
write.table(cistrome_meta, "cistrome_meta.txt", sep = '\t', 
            row.names = F, col.names = T, quote = F)

cistrome_meta_Epi = cistrome_meta[cistrome_meta$type != "TF",]
cistrome_meta_TF = cistrome_meta[cistrome_meta$type == "TF",]

cistrome_meta_Epi$Factor = factor(
  cistrome_meta_Epi$Factor, levels = rev(cistrome_meta_Epi$Factor))
cistrome_meta_TF$Factor = factor(
  cistrome_meta_TF$Factor, levels = rev(cistrome_meta_TF$Factor))

p = list()
p[[1]] = ggplot(cistrome_meta_Epi, 
                aes(x = -log10(phyper.adjust), y = Factor)) +
  geom_bar(stat = "identity", aes(color = log(OR), fill = log(OR))) +
  theme_classic() + theme(aspect.ratio = 0.5) +
  scale_color_gradient(low = "white", high = "red", limits = c(0,2)) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0,2)) +
  labs(x = TeX("$-log_{10}(P_{adj})$"), y = "", color = "log(OR)")

p[[2]] = ggplot(cistrome_meta_TF[1:10,], 
                aes(x = -log10(phyper.adjust), y = Factor)) +
  geom_bar(stat = "identity", aes(color = log(OR), fill = log(OR))) +
  theme_classic() + theme(aspect.ratio = 0.5) +
  scale_color_gradient(low = "white", high = "red", limits = c(0,2)) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0,2)) +
  labs(x = TeX("$-log_{10}(P_{adj})$"), y = "", color = "log(OR)")

p12 = ggarrange(plotlist = p[1:2],
                nrow = 1, ncol = 2,
                labels = c("A","B"),
                font.label = list(size = 40/3, color = "black", face = "bold", 
                                  family = NULL),
                common.legend = T, legend = "bottom",
                align = "hv")

#2. linear model
variant_info_in = foreach(i = 1:22, .combine = "rbind") %dopar%
{
  read.delim(gzfile(sprintf(
    "results/variant_info_in_chr%d.txt.gz", i)), 
    stringsAsFactors = F, check.names = F)
}
variant_info_in$trans_genes_TCGA[is.na(variant_info_in$trans_genes_TCGA)] = ""

coding_meta = read.delim("coding_meta.txt", header = F, row.names = 1)

mcmc_n = 10000
data = data.frame(
  beta2_mean = (variant_info_in$beta_sd^2 + variant_info_in$beta_mean^2 * 
                  mcmc_n / (mcmc_n-1)) * (mcmc_n-1) / mcmc_n,
  MAF = log(variant_info_in$Freq1 * (1 - variant_info_in$Freq1)),
  LD_score = log(variant_info_in$LD_score))

tmp = separate_longer_delim(variant_info_in[,c(1,19)], gene_based_anno, ',')
tmp$value = 1
tmp = spread(tmp, gene_based_anno, value, fill = 0)[,-2]
rownames(tmp) = tmp$alternate_ids
tmp = tmp[variant_info_in$alternate_ids,]
data = cbind(data, tmp[,-1])

for(j in names(variant_info_in)[21:27])
{
  data[[j]] = variant_info_in[[j]] != ""
}

tmp = separate_longer_delim(variant_info_in[,c(1,28)], features, ',')
tmp$value = 1
tmp = spread(tmp, features, value, fill = 0)[,-2]
rownames(tmp) = tmp$alternate_ids
tmp = tmp[variant_info_in$alternate_ids,]
data = cbind(data, tmp[,-1])
data[,4:78] = data[,4:78] * 1

fmla = as.formula(paste("I(log(beta2_mean)) ~ ", 
                        paste(names(data)[c(2:78)], collapse = '+')))
lm = lm(fmla, data)

lm_coefs = data.frame(summary(lm)$coefficients, check.names = F)
lm_coefs$P.adjust = p.adjust(lm_coefs$`Pr(>|t|)`)

lm_coefs = data.frame(summary(lm)$coefficients, check.names = F)
lm_coefs$Annotation = rownames(lm_coefs)
lm_coefs$Annotation[2:3] = c("log(f(1-f))","log(LD score)")
lm_coefs$`Annotation type` = ""
lm_coefs = lm_coefs[,c(5:6,1:4)]
lm_coefs$Annotation[4:14] = coding_meta[lm_coefs$Annotation[4:14], 1]
lm_coefs$`Annotation type`[4:14] = "Gene-based"
lm_coefs$Annotation[15:21] = c("Cis eQTL (TCGA)", "Trans eQTL (TCGA)",
                                "Cis eQTL (GTEx)", "Enhancer (Hi-C,RWPE1)",
                                "Enhancer (Hi-C,C42B)", "Enhancer (Hi-C,22Rv1)",
                                "Enhancer (H3K27ac HiChIP, LNCaP)")
lm_coefs$`Annotation type`[15:17] = "eQTL"
lm_coefs$`Annotation type`[18:21] = "Enhancer"
lm_coefs$`Annotation type`[c(22:31,33:43,52:77)] = "TF ChIP-seq"
lm_coefs$`Annotation type`[c(44:51)] = "HM ChIP-seq"
lm_coefs$P.adjust = p.adjust(lm_coefs$`Pr(>|t|)`)
lm_coefs = arrange(lm_coefs, `Pr(>|t|)`)
write.table(lm_coefs, "lm_coefs.txt", sep = '\t',
            row.names = T, col.names = T, quote = F)

sub = lm_coefs[lm_coefs$P.adjust < 0.1,][-1,]
l = 1
u = -log10(sub$P.adjust[1])
sub$Annotation = factor(sub$Annotation, levels = sub$Annotation)

p[[3]] = ggplot(sub, aes(x = Annotation, y = Estimate, 
                         color = -log10(P.adjust), fill = -log10(P.adjust))) +
  geom_bar(stat = "identity") + theme_classic() + 
  theme(aspect.ratio = 0.2,
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_color_gradient(low = "white", high = "blue", limits = c(l,u),
                       trans = log_trans(),
                       breaks = log_breaks()) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(l,u),
                      trans = log_trans(),
                      breaks = log_breaks()) +
  guides(fill = "none") +
  labs(x = "", y = "Effect size", color = TeX("$-log_{10}(P_{adj})$"))

ggsave("variant_enrichment.pdf",
       ggarrange(p12, p[[3]],
                 nrow = 2, ncol = 1,
                 labels = c("","C"),
                 heights = c(1,1.2),
                 font.label = list(size = 40/3, color = "black", face = "bold", 
                                   family = NULL)),
       device = "pdf", width = 8, height = 6, units = "in")