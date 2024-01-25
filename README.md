# Fine-mapping causal variants of prostate cancer using h2-D2

This repository contains data and scripts to analyze fine-mapping results of prostate cancer (PrCa). The GWAS summary data is from a large meta-analysis of European ancestry (http://practical.icr.ac.uk/blog/?page_id=8164).

## LD_block

`LD_block.R`: Partition SNPs into nearly independent LD blocks for reference 
panel of UK Biobank White British individuals.

`ldetect-data/EUR`: LD blocks partitioned by LDetect based on reference panel 
1000 Genome EUR population.

`reblock.txt`: Coordinates of LD blocks.

## results

`variant_info_chr*.txt.gz`: Full list of variants in prostate cancer fine-mapping analysis. For each variant, we report the information contained in meta-analysis summary data, h2-D2 fine-mapping summary statistics, as well as detailed annotations.

`variant_info_in_chr*.txt.gz`: List of tag variants in prostate cancer fine-mapping analysis. For each tag variant, the annotations are aggregated 
across all variants tag with it.

`block_info_chr*.txt`: For each block, we reported the number of variants with
CL > 0.9, the number of variants with CL > 0.5, the number of detected 95% CSs,
the number of variants in 95% CSs, the number of tag variants in 95% CSs, risk
status (i.e. contain at least one SNP with $P < 5\times 10^{-8}$), posterior mean
and standard error of local heritability (estimated by $\beta^{\top} R \beta$).

`variant_info_CS95.csv`: Table S5 of manuscript.

## variant_annotation_enrichment

`variant_annotation_enrichment.R`: Code used to conduct functional enrichment analysis for credible causal variants.

`variant_annotation_enrichment.pdf`: Figure 2 of manuscript.

`cistrome_meta.txt`: Table S6 of manuscript.

`lm_coefs.txt`: Table S7 of manuscript.

`coding_meta.txt`: Description of gene-based annotations.

## gene_pathway_enrichment

`gene_pathway_enrichment.R`: Code used to conduct pathway enrichment analysis for genes associated with credible causal variants.

`CS95_genes.txt`: List of 369 putative target genes.

`enrich-input1-[GO_BP][WikiPathways].tsv`: Enrichment of putative target genes within pathways from GO Biological Process / WikiPathways. Enrichment analyses were conducted using [GeneCodis 4](https://genecodis.genyo.es/).

`gene_pathway_enrichment.txt`: Enrichment results of putative target genes after merging results of two databases and removing reduntant pathways.

`gene_pathway_enrichment.pdf`: Figure 4 of manuscript.

# Citation

Li, X., Sham, P. C., & Zhang, Y. D. (2023). A Bayesian fine-mapping model using a continuous global-local shrinkage prior with applications in prostate cancer analysis. The American Journal of Human Genetics.
https://doi.org/10.1016/j.ajhg.2023.12.007

# Contact the Author

Xiang Li (freddyl@connect.hku.hk)
