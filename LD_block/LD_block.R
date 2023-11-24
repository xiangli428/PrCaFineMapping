library(bigsnpr)
library(pheatmap)
library(Matrix)
library(foreach)
library(doParallel)

setwd("~/Documents/GWAS/data/UKBB/imputed_genotype_1001_unique_info-0.8_geno-0.01_maf-0.01_hwe-1e-6")

# Define function
split_cost = function(ld)
{
  ld = as.matrix(ld)
  M = nrow(ld)
  cost = rep(0, M-1)
  cost[1] = sum(ld[1,2:M]^2)
  for(j in 2:(M-1))
  {
    cost[j] = cost[j-1] - sum(ld[j,1:j-1]^2) + sum(ld[j,(j+1):M]^2)
  }
  cost = cost / (c(1:(M-1)) * c((M-1):1))
  return(cost)
}

recursive_break = function(ld, start, end)
{
  bps = c()
  
  if(end - start > 100)
  {
    cost = split_cost(ld[start:end,start:end])
    bp = which.min(cost[51:length(cost)]) + 50
    if(cost[bp] < 1e-3 && (length(cost) - bp >= 50))
    {
      bps = c(bps, bp + start - 1)
      bps = c(bps, recursive_break(ld, start, start+bp-1))
      bps = c(bps, recursive_break(ld, start+bp, end))
    }
  }
  return(sort(bps))
}

# Partition
i = 22 # ID of chromosome
data_rds = snp_attach(paste("genotype_rds/chr",i,".rds", sep = ''))

# Split the LD block.
breakpoint = read.delim(
  paste("ldetect-data/EUR/fourier_ls-chr",i,".bed",
        sep = ''), stringsAsFactors = F, check.names = F)
breakpoint$size = breakpoint$stop - breakpoint$start
breakpoint$start = breakpoint$start + 1

variant$block = NA
breakpoint$num_variant = NA
for(j in 1:nrow(breakpoint))
{
  variant$block[variant$position >= breakpoint[j,2] &
                  variant$position <= breakpoint[j,3]] = j
  breakpoint$num_variant[j] = sum(variant$block == j, na.rm = T)
}

# Reblock.
variant$reblock = NA
nblocks = nrow(breakpoint)
start = 1
end = breakpoint$num_variant[1]
j = 2
k = 1

data_cor = snp_cor(data_rds$genotypes,
                   ind.col = start:end,
                   size = end - start,
                   thr_r2 = 2.5e-3,
                   fill.diag = F,
                   ncores = 80)

while(T)
{
  cost = split_cost(data_cor)
  if(length(cost) > 50)
  {
    bp = which.min(cost[51:length(cost)]) + 50
  }
  
  if(length(cost) > 50 && cost[bp] < 1e-3 && (length(cost) - bp >= 50))
    # accept this breakpoint
  {
    bps = c(0, recursive_break(data_cor[1:bp,1:bp],1,bp), bp)
    
    for(l in 1:(length(bps)-1))
    {
      variant$reblock[(start + bps[l]):(start + bps[l+1] - 1)] = k
      dir.create(paste("LD_size-1000_thr-0.0025/chr",i,"/",k, sep = ''),
                 recursive = T)
      writeMM(data_cor[(bps[l]+1):bps[l+1], (bps[l]+1):bps[l+1]],
              paste("LD_size-1000_thr-0.0025/chr",i,"/",k,"/LD_0.05.mtx", sep = ''))
      system(paste("gzip ","LD_size-1000_thr-0.0025/chr",i,"/",k,"/LD_0.05.mtx", sep = ''))
      write.table(
        variant[(start + bps[l]):(start + bps[l+1] - 1),],
        paste("LD_size-1000_thr-0.0025/chr",i,"/",k,"/variant_reblock.txt", sep = ''),
        sep = '\t', row.names = F, col.names = T, quote = F)
      write.table(
        variant[(start + bps[l]):(start + bps[l+1] - 1),2],
        paste("LD_size-1000_thr-0.0025/chr",i,"/",k,"/variant.txt", sep = ''),
        sep = '\t', row.names = F, col.names = F, quote = F)
      k = k + 1
    }
    
    start = start + bp
    data_cor = data_cor[(bp+1):nrow(data_cor),(bp+1):nrow(data_cor)]
  } else if (j <= nblocks) {
    end = end + breakpoint$num_variant[j]
    j = j + 1
    data_cor = snp_cor(data_rds$genotypes,
                       ind.col = start:end,
                       size = end - start,
                       thr_r2 = 2.5e-3,
                       fill.diag = F,
                       ncores = 80)
  } else {
    variant$reblock[start:end] = k
    dir.create(paste("LD_size-1000_thr-0.0025/chr",i,"/",k, sep = ''),
               recursive = T)
    writeMM(data_cor,
            paste("LD_size-1000_thr-0.0025/chr",i,"/",k,"/LD_0.05.mtx", sep = ''))
    system(paste("gzip ","LD_size-1000_thr-0.0025/chr",i,"/",k,"/LD_0.05.mtx", sep = ''))
    write.table(
      variant[start:end,],
      paste("LD_size-1000_thr-0.0025/chr",i,"/",k,"/variant_reblock.txt", sep = ''),
      sep = '\t', row.names = F, col.names = T, quote = F)
    write.table(
      variant[start:end,2],
      paste("LD_size-1000_thr-0.0025/chr",i,"/",k,"/variant.txt", sep = ''),
      sep = '\t', row.names = F, col.names = F, quote = F)
    break
  }
}

reblock = data.frame("ID" = 1:k,
                     "start" = 0,
                     "stop" = 0,
                     "size" = 0)
for(j in 1:k)
{
  reblock[j,2] = variant$position[min(which(variant$reblock == j))]
  reblock[j,3] = variant$position[max(which(variant$reblock == j))]
}
reblock[,4] = reblock[,3] - reblock[,2] + 1
reblock$num_variants = as.integer(table(variant$reblock))

write.table(reblock,
            sprintf("LD_size-1000_thr-0.0025/chr%d/reblock.txt", i),
            sep = '\t', row.names = F, col.names = T, quote = F)