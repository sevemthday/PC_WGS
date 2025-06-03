library(optparse)
library(data.table)
library(tidyverse)
# optparse variable
#-----------------------------------------------------------------------------------------------------#
option_list <- list(
  make_option(c("-i", "--input"),      action = "store", type = "character", default = NULL, help = "annovar input file name"),
  make_option(c("-o", "--output"),     action = "store", type = "character", default = NULL, help = "output file name")
)
opt <- parse_args(OptionParser(option_list = option_list))

# function
CLNSIGCONF.summary <- function(x){
  if(x[1] != "."){
  score <- as.integer(gsub(".*\\((.*)\\).*", "\\1", x))
  names(score) <- gsub("\\((.*)\\).*", "", x)
  
  score.type <- c()
  score.type[1] <- sum(score[names(score) %in% c("Benign", "Likely_benign")])
  score.type[2] <- sum(score[names(score) %in% c("Uncertain_significance")])
  score.type[3] <- sum(score[names(score) %in% c("Pathogenic", "Likely_pathogenic")])
  names(score.type) <- c("Benign/Likely_benign", "Uncertain_significance", "Pathogenic/Likely_pathogenic")
  
  score.type.major <- names(score.type)[score.type == max(score.type)]
  if(length(score.type.major) > 1){
    return(paste(score.type.major, collapse = ";"))
  }else{
    return(score.type.major)
  }
  }else{
    return("")
  }
}

# input data
dt <- fread(opt$input)
dt$ClinVar.CLNSIGCONF.major <- sapply(strsplit(dt$ClinVar_CLNSIGCONF, split = "\\|_"), CLNSIGCONF.summary)
# output data
write.table(dt, opt$output, sep = "\t", quote = F, row.names = F)
