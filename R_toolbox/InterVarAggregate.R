#version:v20220809
library(optparse)
library(data.table)
library(tidyverse)
# optparse variable
#-----------------------------------------------------------------------------------------------------#
option_list <- list(
  make_option(c("-i", "--input"),      action = "store", type = "character", default = NULL, help = "intervar input file name"),
  make_option(c("-o", "--output"),     action = "store", type = "character", default = NULL, help = "output file name")
)
opt <- parse_args(OptionParser(option_list = option_list))


Convert_intervar_ACMG <- function(dt){
  col <- c("#Chr", "Start", "End", "Ref", "Alt", "Ref.Gene", "InterVar: InterVar and Evidence")
  dt <- dt[,..col]
  colnames(dt)[colnames(dt)=="#Chr"] <- "Chr"
  colnames(dt)[colnames(dt)=="Ref.Gene"] <- "Ref.Gene.InterVar"
  colnames(dt)[colnames(dt)=="InterVar: InterVar and Evidence"] <- "ACMG.rule.InterVar"
  
  ACMG <- dt$ACMG.rule.InterVar
  dt$ACMG.InterVar <- gsub("InterVar: ", "", sapply(strsplit(ACMG , split = " PVS1="), `[`, 1))
  
  ACMG <- sapply(strsplit(ACMG , split = " PVS1="), `[`, 2)
  ACMG <- gsub(" PS=\\["   , ", ", ACMG)
  ACMG <- gsub("\\] PM=\\[", ", ", ACMG)
  ACMG <- gsub("\\] PP=\\[", ", ", ACMG)
  ACMG <- gsub("\\] BA1="  , ", ", ACMG)
  ACMG <- gsub(" BS=\\["   , ", ", ACMG)
  ACMG <- gsub("\\] BP=\\[", ", ", ACMG)
  ACMG <- gsub("\\]"       , ""  , ACMG)
  ACMG <- strsplit(ACMG, split = ", ") %>% unlist()
  ACMG <-matrix(ACMG, nrow = length(dt$ACMG.rule.InterVar), ncol = 33, byrow = TRUE)
  colnames(ACMG) <-c("PVS1", paste0("PS",1:5), paste0("PM",1:7), paste0("PP",1:6), "BA1", paste0("BS",1:5), paste0("BP",1:8))
  dt$ACMG.rule.InterVar <- apply(ACMG, 1, function(x) paste(colnames(ACMG)[grepl("1", x)], collapse = ","))
  dt <- dt[, list(Ref.Gene.InterVar = paste(Ref.Gene.InterVar, collapse=";"),
                      ACMG.InterVar = paste(ACMG.InterVar,     collapse=";"),
                 ACMG.rule.InterVar = paste(ACMG.rule.InterVar,collapse=";")), by = list(Chr, Start, End, Ref, Alt)]
  return(dt)
}

ACMG.summary <- function(x) {
  x <- unique(x)
  if(length(x) == 2){
    if(all(x %in% c("Benign", "Likely benign"))){
      x <- "Benign/Likely benign"
    }else if(all(x %in% c("Pathogenic", "Likely pathogenic"))){
      x <- "Pathogenic/Likely pathogenic"
    }else{
      x <- "Conflicting interpretations of pathogenicity"
    }
  }else if(length(x) > 2){
    x <- "Conflicting interpretations of pathogenicity"
  }
  return(x)
}

ACMG.merge <- function(x) {
  x <- unique(x)
  x <- paste(x, collapse = ";")
  return(x)
}


# input data
dt <- fread(opt$input)
dt <- Convert_intervar_ACMG(dt)
# ACMG.intervar.summary
dt$ACMG.InterVar.summary <- sapply(strsplit(dt$ACMG.InterVar, split = ";"), ACMG.merge)
# add chromosome prefix
dt$Chr <- ifelse(grepl("^HLA", dt$Chr), dt$Chr, paste0("chr", dt$Chr))

write.table(dt, opt$output, sep = "\t", quote = F, row.names = F)
