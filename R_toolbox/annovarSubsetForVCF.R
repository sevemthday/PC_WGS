#version:v20230815
library(optparse)
library(data.table)

# optparse variable
#-----------------------------------------------------------------------------------------------------#
option_list <- list(
  make_option(c("-I", "--input"),      action = "store", type = "character", default = NULL, help = "annovar input file name"),
  make_option(c("-O", "--output"),     action = "store", type = "character", default = NULL, help = "output file name")
)
opt <- parse_args(OptionParser(option_list = option_list))

# input data
dt <- fread(opt$input)
start = which(colnames(dt) == "Otherinfo4")
end = which(colnames(dt) == "Segregation") -1
dt <- dt[, c(start:end), with = FALSE]

write.table(dt, opt$output, sep = "\t", quote = F, row.names = F, col.names = F)