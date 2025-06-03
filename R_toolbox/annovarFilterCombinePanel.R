#version:v20230521
#""" remove intervar filtering """
#""" combine variants after panel filtering """
#""" remove print() info """
library(optparse)
library(data.table)
library(tidyverse)

#-----------------------------------------------------------------------------------------------------#
# optparse variable
option_list <- list(
  make_option(c("-i", "--input"),      action = "store", type = "character", default = NULL, help = "annovar input file name"),
  make_option(c("-o", "--output"),     action = "store", type = "character", default = NULL, help = "prefix for output file name"),
  make_option(c("--gnomad"),           action = "store", type = "logical",   default = FALSE,help = "Genome aggregation database filter, T or F"),
  make_option(c("--twbiobank"),        action = "store", type = "logical",   default = FALSE,help = "Taiwan biobank database filter, T or F"),
  make_option(c("--gnomadCutoff"),     action = "store", type = "double",    default = 0.05, help = "cut-off for Genome aggregation database"),
  make_option(c("--twbiobankCutoff"),  action = "store", type = "double",    default = 0.05, help = "cut-off for Taiwan biobank database"),
  make_option(c("--regionRefGene"),    action = "store", type = "character", default = NULL, help = "Func.refGene: exonic:splicing:ncRNA:UTR5:UTR3:intronic:upstream:downstream:intergenic"),
  make_option(c("--remvSynonymousSNV"),action = "store", type = "logical",   default = FALSE,help = "remove synonymous SNV, T or F"),
  make_option(c("--keepPASS"),         action = "store", type = "logical",   default = FALSE,help = "select PASS, T or F"),
  make_option(c("--panelDir"),         action = "store", type = "character", default = NULL, help = "path for gene panel"),
  make_option(c("--clinvarPLP"),       action = "store", type = "logical",   default = FALSE,help = "output of ClinVar P/LP, T or F"),
  make_option(c("--spliceai"),         action = "store", type = "logical",   default = FALSE,help = "filter spliceai, T or F"),
  make_option(c("--spliceaiCutoff"),   action = "store", type = "character", default = NULL, help = "cut-off for filter spliceai, high_precision,	recommended, high_recall"),
  make_option(c("--sampleID"),         action = "store", type = "character", default = NULL, help = "sampleID file")
)
opt <- parse_args(OptionParser(option_list = option_list))

# function
#count variant for each step
variantCounter <- function(vci, dt, filterQuiteria){
  vci[length(vci)+1] <- nrow(dt) #vci: variant Count Info
  names(vci)[length(vci)] <- filterQuiteria
  return(vci)
}

#-----------------------------------------------------------------------------------------------------#
### read file ###
# read table
dt <- fread(opt$input)
# read sample names
if(!(is.null(opt$sampleID))){
  sampleID <- readLines(opt$sampleID)
  colnames(dt)[grep("^Otherinfo13", colnames(dt)):(grep("^Otherinfo13", colnames(dt))+length(sampleID)-1)] <- sampleID
}
# count variants - orignial dt
variantCountInfo <- c()
variantCountInfo <- variantCounter(variantCountInfo, dt, "total")

#-----------------------------------------------------------------------------------------------------#
### ClinVar P/LP ###
#dt -> ClinVar P/LP
if(opt$clinvarPLP){
  dt.clinvarPLP <- dt[grepl("Pathogenic|Likely_pathogenic", CLNSIG)]
  write.table(dt.clinvarPLP, paste0(opt$output, "_clinvar.PLP.txt"), sep = "\t", quote = F, row.names = F)
}

#-----------------------------------------------------------------------------------------------------#
### SpliceAI filtering ###
#dt -> spliceAI score -> gnomAD/twbiobank -> PASS -> panel
if(opt$spliceai){
  if(opt$spliceaiCutoff == "high_precision"){
    #dt <- dt[high_precision == "True"]
    dt <- dt[grepl("True|TRUE", high_precision)]
  }else if(opt$spliceaiCutoff == "recommended"){
    #dt <- dt[recommended == "True"]
    dt <- dt[grepl("True|TRUE", recommended)]
  }else if(opt$spliceaiCutoff == "high_recall"){
    #dt <- dt[high_recall == "True"]
    dt <- dt[grepl("True|TRUE", high_recall)]
  }
  variantCountInfo <- variantCounter(variantCountInfo, dt, "spliceai")
}

#-----------------------------------------------------------------------------------------------------#
### general variant filtering ###
#dt -> spliceAI score -> gnomAD/twbiobank -> PASS -> panel
# gnomAD
if(opt$gnomad){
  gnomadCutoff <- opt$gnomadCutoff
  dt <- dt[AF_gnomad30_genome == "." | as.numeric(AF_gnomad30_genome) < gnomadCutoff]
  variantCountInfo <- variantCounter(variantCountInfo, dt, "gnomAD")
}
# twbiobank
if(opt$twbiobank){
  twbiobankCutoff <- opt$twbiobankCutoff
  dt <- dt[`TaiwanBiobank-official_Illumina1000-AF` == "." | as.numeric(`TaiwanBiobank-official_Illumina1000-AF`) < twbiobankCutoff]
  variantCountInfo <- variantCounter(variantCountInfo, dt, "twbiobank")
}
# region
if(!(is.null(opt$regionRefGene))){
  regionRefGene <- opt$regionRefGene
  regionRefGene <- unlist(strsplit(regionRefGene, split = ":"))
  regiontype <- c("exonic", "splicing", "ncRNA", "UTR5", "UTR3", "intronic", "upstream", "downstream", "intergenic")
  regiontype.filter <- regiontype[!(regiontype %in% regionRefGene)]
  dt <- dt[!(grepl(paste(regiontype.filter, collapse="|"), Func.refGene))]
  variantCountInfo <- variantCounter(variantCountInfo, dt, paste0("region:", opt$regionRefGene))
}
# remove synonymous SNV
if(opt$remvSynonymousSNV){
  dt <- dt[ExonicFunc.refGene != "synonymous SNV"]
  variantCountInfo <- variantCounter(variantCountInfo, dt, paste0("remove synonymous SNV"))
}

# filter PASS
if(opt$keepPASS){
  dt <- dt[Otherinfo10 == "PASS"]
  variantCountInfo <- variantCounter(variantCountInfo, dt, "Otherinfo10:PASS")
}
# panel filter
if(!(is.null(opt$panelDir))){
  dt.geneset <- list()
  panel_files <- list.files(path = opt$panelDir)
  for (i in seq(panel_files)) {
    geneset <- fread(paste0(opt$panelDir, "/", panel_files[i]))
    refGene.list <- strsplit(dt$Gene.refGene, split = ";")
    refGene.list.seq <- rep(seq_along(refGene.list), sapply(refGene.list, length))
    dt.temp <- dt[unique(refGene.list.seq[unlist(refGene.list) %in% geneset$Gene])]
    dt.temp$Panel_filtering <- gsub(".txt", "", panel_files[i])
    dt.geneset[[i]] <- dt.temp
    variantCountInfo <- variantCounter(variantCountInfo, dt.geneset[[i]], gsub(".txt", "", panel_files[i]))
  }
  dt.geneset <- rbindlist(dt.geneset)
  cols_chosen <- colnames(dt.geneset)[colnames(dt.geneset) != "Panel_filtering"]
  dt.geneset <- dt.geneset[, list(Panel_filtering = paste(Panel_filtering, collapse = ";")), by = cols_chosen]
  variantCountInfo <- variantCounter(variantCountInfo, dt.geneset, "Combine all panel filtering")
  write.table(dt.geneset, paste0(opt$output, "_panel.filtering.txt"), sep = "\t", quote = F, row.names = F)
}


write.table(data.table(names(variantCountInfo), variantCountInfo), paste0(opt$output, "_variantCount.txt"),
            sep = "\t", quote = F, row.names = F)
