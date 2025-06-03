#version:v20230802
library(optparse)
library(data.table)
library(tidyverse)
# optparse variable
#-----------------------------------------------------------------------------------------------------#
option_list <- list(
  make_option(c("-i", "--input"),      action = "store", type = "character", default = NULL, help = "annovar input file name"),
  make_option(c("-I", "--Input"),      action = "store", type = "character", default = NULL, help = "intervar input file name"),
  make_option(c("-o", "--output"),     action = "store", type = "character", default = NULL, help = "output file name"),
  make_option(c("--all"),              action = "store_true", default = FALSE, help = "keep all rows (variants)"),
  make_option(c("--all_x"),            action = "store_true", default = FALSE, help = "keep first data rows (variants)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# input data
dt1 <- fread(opt$input)
dt2 <- fread(opt$Input)
dt1[, UID := paste(Chr, Start, Ref, Alt, sep = "-")]
dt2[, UID := paste(Chr, Start, Ref, Alt, sep = "-")]
dt2[,c("Chr", "Start", "End", "Ref", "Alt") := NULL]

if(opt$all){
  dt <- merge(dt1, dt2, by = "UID", all = TRUE)
}else if(opt$all_x){
  dt <- merge(dt1, dt2, by = "UID", all.x = TRUE)
}else{
  dt <- merge(dt1, dt2, by = "UID")
}

write.table(dt, opt$output, sep = "\t", quote = F, row.names = F)