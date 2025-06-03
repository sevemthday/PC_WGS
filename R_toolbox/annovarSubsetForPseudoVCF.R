#version:v20241129
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
dt <- fread(opt$input, colClasses = "character")
start = which(colnames(dt) == "Otherinfo4")
end = which(colnames(dt) == "Otherinfo12")
dt <- dt[, c(start:end), with = FALSE]

colnames(dt) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

dt$SAMPLE <- ifelse(dt$FORMAT == "GT:AD:DP:GQ:PL", "0/1:100,100:200:99:100,0,150", paste0("0|1:100,100:200:99:0|1:" , dt$POS, "_", dt$REF, "_", dt$ALT, ":100,0,150:12345"))

write.table(dt, opt$output, sep = "\t", quote = F, row.names = F, col.names = T)