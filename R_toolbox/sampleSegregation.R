#version: 20220814
library(optparse)
library(data.table)
library(tidyverse)
# optparse variable
#-----------------------------------------------------------------------------------------------------#
option_list <- list(
  make_option(c("-i", "--input"),      action = "store", type = "character", default = NULL, help = "annovar input file name"),
  make_option(c("-o", "--output"),     action = "store", type = "character", default = NULL, help = "prefix for output file name"),
  make_option(c("--sampleID"),         action = "store", type = "character", default = NULL, help = "sampleID file")
)
opt <- parse_args(OptionParser(option_list = option_list))

# input data
dt <- fread(opt$input)
sampleID <- readLines(opt$sampleID)

# change sample name column
colnames(dt)[grep("^Otherinfo13", colnames(dt)):(grep("^Otherinfo13", colnames(dt))+length(sampleID)-1)] <- sampleID

# sample segregation
if(length(sampleID)==1){
  dt$Segregation <- sampleID
}else{
  dt_temp <- dt[, ..sampleID] %>% as.matrix()
  matrix_boolean_temp <- sapply(1:nrow(dt_temp), function(x) grepl(pattern = "1", sapply(strsplit(dt_temp[x,], split=":"), `[`, 1)))
  row.names(matrix_boolean_temp) <- sampleID
  dt$Segregation <- apply(matrix_boolean_temp, 2, function(x) paste(names(which(x)), collapse="," ))
}

write.table(dt, opt$output, sep = "\t", quote = F, row.names = F)