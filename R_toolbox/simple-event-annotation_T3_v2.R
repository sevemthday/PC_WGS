library(optparse)
# variable
option_list <- list(
  make_option(c("-i", "--input"),  action = "store", type = "character", default = NULL, help = "VCF input file name"),
  make_option(c("-o", "--output"), action = "store", type = "character", default = NULL, help = "output file name"),
  make_option(c("-g", "--genome"), action = "store", type = "character", default = FALSE,help = "genome reference builds. hg19 or hg38")
)
opt <- parse_args(OptionParser(option_list = option_list))



# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("StructuralVariantAnnotation")
#install.packages("stringr")
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)
library(data.table)
#' Simple SV type classifier
simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomosal
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
          "DUP")))))
}
# using the example in the GRIDSS /example directory
vcf <- readVcf(opt$input, opt$genome)
info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), 
                                    data.frame(row.names=c("SIMPLE_TYPE"),
                                               Number=c("1"),
                                               Type=c("String"),
                                               Description=c("Simple event type annotation based purely on breakend position and orientation.")),
                                    data.frame(row.names=c("SVLEN"),
                                               Number=c("1"),
                                               Type=c("Integer"),
                                               Description=c("SV_length"))), "DataFrame"))


gr <- breakpointRanges(vcf)
svtype <- simpleEventType(gr)
info(vcf[gr$sourceId])$SVTYPE <- svtype
info(vcf)$SVLEN <- NA_character_
info(vcf[gr$sourceId])$SVLEN <- gr$svLen
alt(vcf)[info(vcf)$SVTYPE != "BND"] <- paste0("<", info(vcf)$SVTYPE, ">")[info(vcf)$SVTYPE != "BND"]

#subset vcf
dt <- data.table(MATEID = row.names(vcf), CHROM = as.character(vcf@rowRanges@seqnames), POS = vcf@rowRanges@ranges@start)
dt$ID <- sapply(1:nrow(dt), function(x) substr(dt$MATEID[x], 1, nchar(dt$MATEID[x]) - 1))
dt[, POS.min := min(POS), by = ID]
vcf_subset <- vcf[(dt$POS == dt$POS.min) & (info(vcf)$SVTYPE != "BND")| (info(vcf)$SVTYPE == "BND")| (info(vcf)$SVTYPE == "CTX")] #BND  CTX  DEL  DUP  INS  INV
writeVcf(vcf_subset, opt$output)

