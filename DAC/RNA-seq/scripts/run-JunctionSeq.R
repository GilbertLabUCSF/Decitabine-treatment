args <- commandArgs(trailingOnly = TRUE)

PDIR <- args[1]
samplesheetPATH = args[2]
COND <- args[3]
refCOND <- args[4]
name <- args[5]

setwd(PDIR)
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(JunctionSeq)))
message ('R libraries loaded!')
suppressMessages(suppressWarnings(dir.create('iso')))
suppressMessages(suppressWarnings(dir.create(paste0("iso/",name))))
suppressMessages(suppressWarnings(dir.create(paste0("iso/",name,"/jscs"))))
suppressMessages(suppressWarnings(dir.create(paste0("iso/",name,"/plots/"))))

loadRData <- function(fileName){
    ## loads an RData file, and returns it 
    # https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
    load(fileName)
    get(ls()[ls() != "fileName"])
}

if (file.exists(paste0("iso/",name,"/jscs/jscs.RData"))){
    jscs <- loadRData(paste0("iso/",name,"/jscs/jscs.RData"))

    message ('`jscs.RData` exists!')
    
} else {
    ###################################################################################################
    message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    message ('Part 1. Read genome annotations, count data and samplesheet')
    message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    message (date())

    GTF = '~/genomes/hg38/gencode.v38/gencode.v38.annotation.gtf'
    # annoFiles ='~/genomes/hg38/gencode.v38/JunctionSeq.flat.gff.gz'
    annoFiles ='/data_gilbert/home/aarab/Projects/Decitabine-treatment/RNA-seq/iso/annoFiles/JunctionSeq.flat.gff.gz'

    gtf <- rtracklayer::import(GTF)

    gene2name <- gtf[gtf$type == "gene"] %>% data.frame %>% 
        column_to_rownames('gene_id') %>% dplyr::select('gene_name') %>% 
        rownames_to_column('geneID')

    txpt2name <- gtf[gtf$type == "transcript"] %>% data.frame %>% 
        column_to_rownames('transcript_id') %>% 
        dplyr::select(c('transcript_name','gene_name'))

    message ('annotations loaded!')

    samplesheet = read.table(
        samplesheetPATH,header=TRUE,stringsAsFactors=FALSE,sep=','
    )

    countFiles <- paste0(
        "iso/rawCts/",
        samplesheet$Sample,
        "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz"
    )
    files = countFiles %>% as.list

    message ("samples:")
    message (samplesheet[,1])
    message ("conditions:")
    message (samplesheet[,1])

    ###################################################################################################
    message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    message ('Part 2. Run JunctionSeq')
    message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    message (date())

    jscs <- runJunctionSeqAnalyses(    
        sample.files = files,
        sample.names = samplesheet[,1],
        condition=relevel(factor(samplesheet[,COND]),ref=refCOND),
        # use.novel.junctions = FALSE,
        flat.gff.file = annoFiles,
        gene.names = gene2name,
        nCores = 15,
        analysis.type = "junctionsAndExons"
    )
    message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    message ('Part 3. Save results')
    message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    message (date())

    writeCompleteResults(
        jscs,
        outfile.prefix=paste0("iso/",name,"/jscs/"),
        save.jscs = TRUE
    )
}

message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
message ('Last Part. Draw plots')
message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
message (date())
buildAllPlots(
    jscs=jscs,
    outfile.prefix = paste0("iso/",name,"/plots/"),
    use.plotting.device = "svg",
    FDR.threshold = 0.01
)

message (date())
