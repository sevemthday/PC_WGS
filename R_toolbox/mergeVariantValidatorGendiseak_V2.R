library(optparse)
library(data.table)

# optparse variable
#-----------------------------------------------------------------------------------------------------#
option_list <- list(
  make_option(c("--inputAnnovar"),          action = "store", type = "character", default = NULL, help = "annovar input file name"),
  make_option(c("--inputVariantValidator"), action = "store", type = "character", default = NULL, help = "annovar input file name"),
  make_option(c("--inputGendiseak"),        action = "store", type = "character", default = NULL, help = "annovar input file name"),
  make_option(c("-o", "--output"),          action = "store", type = "character", default = NULL, help = "output file name")
)
opt <- parse_args(OptionParser(option_list = option_list))

#Annovar
dt <- fread(opt$inputAnnovar)
dt[,UID2 := paste0(Otherinfo4, "-", Otherinfo5, "-", Otherinfo7, "-", Otherinfo8)]

#Gendiseak
dt.gendiseak <- fread(opt$inputGendiseak)
dt.gendiseak[,UID2 := paste0(`#AiViewer_chr`, "-", AiViewer_start, "-", AiViewer_ref , "-", AiViewer_alt)]
col_sel <- c("UID2", "ACMG_Classes", "ACMG_Rules")
dt.gendiseak <- dt.gendiseak[, ..col_sel]
dt <- merge(dt, dt.gendiseak, by = "UID2")

#VariantValidator
dt_vv <- fread(opt$inputVariantValidator, colClasses = "character")
#dt_vv <- dt_vv[Mane_Select == "True"]
#colnames(dt_vv) <- c("UID2", "Representative_transcript",	"Transcript_HGVS",	"P_HGVS_SLC",	"P_HGVS_TLC",	"Mane_Select",	"Gene_ID_vv",	"Symbol_vv", "Impact_vv")
dt <- merge(dt, dt_vv, by = "UID2", all.x = TRUE)

write.table(dt, opt$output, sep = "\t", quote = F, row.names = F)