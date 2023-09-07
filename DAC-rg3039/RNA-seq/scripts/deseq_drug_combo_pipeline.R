args <- commandArgs(trailingOnly = TRUE)

PDIR <- args[1]
GTF <- args[2]
quantDir <- args[3]
outDir <- args[4]

setwd(PDIR)

suppressMessages(suppressWarnings(library (GenomicFeatures)))
suppressMessages(suppressWarnings(library (tximport)))
suppressMessages(suppressWarnings(library (tidyverse)))
suppressMessages(suppressWarnings(library (DESeq2)))
suppressMessages(suppressWarnings(library (BiocParallel)))

register(MulticoreParam(18))


source('../../scripts/util.R')

write_Result <- function(res, name_it, col=FALSE, row=FALSE){
    write.table(res,name_it, sep="\t", quote=FALSE, col.names=col, row.names=row)
}


gtf = rtracklayer::import(GTF)

files <- list.files(path=quantDir, pattern="quant.sf",full.names = TRUE, recursive=T)
names(files) <- gsub(paste0(quantDir,"/(\\S+)/quant.sf"),"\\1",files)

message('quant files:')
for (f in files){message(f)}

txi <- tximport(files, type = "salmon", txOut=T)

# Define sample sheet
# meta 
conds  <- factor(c(
    rep('Combination',2),rep('Decitabine',2),rep('DMSO',2),rep('rg3039',2)
), levels = c('DMSO','Decitabine','rg3039','Combination'))

reps <- factor(c(
    rep(c('rep1','rep2'),4)
),c('rep1','rep2'))

colData <- data.frame(
    cond=conds,
    reps=reps,
    row.names=colnames(txi$abundance),
    stringsAsFactors=FALSE
)


dds0 <- DESeqDataSetFromTximport(txi, colData, ~cond)
dds0 <- estimateSizeFactors(dds0)
ncu <- counts(dds0, normalized=TRUE) 

dds_simple <- DESeqDataSetFromTximport(
    txi, colData, ~0+cond
)
dds_simple <- DESeq(dds_simple, parallel=TRUE)

message(resultsNames(dds_simple))

saveRDS(
    dds_simple,
    paste(outDir,'dds.rds',sep='/')
)

# contrast design: combination treatment vs dmso 
res_comb_vs_dmso  = results(dds_simple, contrast=list(c('condCombination'),c('condDMSO')),listValues=c(1,-1))
message('res_comb_vs_dmso')
message(res_comb_vs_dmso %>% summary)

# contrast design: decitabine treatment vs dmso 
res_decitabine_vs_dmso  = results(dds_simple, contrast=list(c('condDecitabine'),c('condDMSO')),listValues=c(1,-1))
message('res_decitabine_vs_dmso')
message(res_decitabine_vs_dmso %>% summary)

# contrast design: combination treatment vs decitabine
res_comb_vs_decitabine  = results(dds_simple, contrast=list(c('condCombination'),c('condDecitabine')),listValues=c(1,-1))
message('res_comb_vs_decitabine')
message(res_comb_vs_decitabine %>% summary)

# contrast design: rg3039 treatment vs dmso
res_rg3039_vs_dmso  = results(dds_simple, contrast=list(c('condrg3039'),c('condDMSO')),listValues=c(1,-1))
message('res_rg3039_vs_dmso')
message(res_rg3039_vs_dmso %>% summary)

# contrast design: combination treatment vs rg3039 
res_comb_vs_rg3039  = results(dds_simple, contrast=list(c('condCombination'),c('condrg3039')),listValues=c(1,-1))
message('res_comb_vs_rg3039')
message(res_comb_vs_rg3039 %>% summary)

# merge results
RES = list(
    'comb_vs_dmso'=res_comb_vs_dmso %>% data.frame, 
    'comb_vs_decitabine'=res_comb_vs_decitabine %>% data.frame,
    'comb_vs_rg3039'=res_comb_vs_rg3039 %>% data.frame,
    'decitabine_vs_dmso'=res_decitabine_vs_dmso %>% data.frame,
    'rg3039_vs_dmso'=res_rg3039_vs_dmso %>% data.frame
    # 'rep2_vs_rep1'=res_rep2_vs_rep1
)

result_table <- cbind(
    RES[['comb_vs_dmso']] %>% select('log2FoldChange','pvalue') %>% 
        rename(comb_vs_dmso_log2FC=log2FoldChange,comb_vs_dmso_pvalue=pvalue),
    RES[['comb_vs_decitabine']] %>% select('log2FoldChange','pvalue') %>% 
        rename(comb_vs_decitabine_log2FC=log2FoldChange,comb_vs_decitabine_pvalue=pvalue),
    RES[['comb_vs_rg3039']] %>% select('log2FoldChange','pvalue') %>% 
        rename(comb_vs_rg3039_log2FC=log2FoldChange,comb_vs_rg3039_pvalue=pvalue),
    RES[['decitabine_vs_dmso']] %>% select('log2FoldChange','pvalue') %>% 
        rename(decitabine_vs_dmso_log2FC=log2FoldChange,decitabine_vs_dmso_pvalue=pvalue),
    RES[['rg3039_vs_dmso']] %>% select('log2FoldChange','pvalue') %>% 
        rename(rg3039_vs_dmso_log2FC=log2FoldChange,rg3039_vs_dmso_pvalue=pvalue)
    # RES[[6]] %>% select('log2FoldChange','pvalue') %>% 
    #     rename(rep2_vs_rep1_log2FC=log2FoldChange,rep2_vs_rep1_pvalue=pvalue)
) %>% drop_na

result_table <- result_table %>% replace(is.na(result_table),0)

# save results
write_Result(
    ncu,paste(outDir,'deseq2_norm.txt',sep='/'),col=T,row=T
)
write_Result(
    result_table,paste(outDir,'result_table.txt',sep='/'),col=T,row=T
)