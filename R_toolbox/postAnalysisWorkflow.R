library(data.table)
library(optparse)

# workflow functions script
function_path <- "/staging/biology/sevemthday77/software/Rscript_/functions/"
source(paste0(function_path, "parseGqDpAf_func.R"))
source(paste0(function_path, "parseRoleInCancer_func.R"))
source(paste0(function_path, "pvs1Annotation_func.R"))
source(paste0(function_path, "unipDomainAnnotation_func.R"))

option_list <- list(
  make_option(c("--input"),         action = "store", type = "character", default = NULL, help = "annovar input file name"),
  make_option(c("--output"),        action = "store", type = "character", default = NULL, help = "output file name"),
  make_option(c("--panel"),         action = "store", type = "character", default = NULL, help = "cancer gene panel with classification"),
  make_option(c("--maneInfo"),      action = "store", type = "character", default = NULL, help = "MANE select info"),
  make_option(c("--unipDomainBed"), action = "store", type = "character", default = NULL, help = "UniProt domain bed file")
)

opt <- parse_args(OptionParser(option_list = option_list))


# main function
runWorkflow <- function(input, output, RoleInCancer_file, mane_info_file, UniProt_Domain_Bed_file) {
  
  dt <- fread(input, colClasses = "character")
  
  # Step 1: parseGqDpAf
  dt <- parseGqDpAf(dt)
  
  # Step 2: parseRoleInCancer
  dt <- parseRoleInCancer(dt, RoleInCancer_file)
  
  # Step 3: pvs1Annotation
  dt <- pvs1Annotation(dt, mane_info_file)
  
  # Step 4: unipDomainAnnotation
  dt <- unipDomainAnnotation(dt, UniProt_Domain_Bed_file) 
  
  write.table(dt, output, quote = F, sep = "\t", row.names = F)
  
  message("post-analysis complete! Output written to ", output)
}


# run main function
runWorkflow(opt$input, opt$output, opt$panel, opt$maneInfo, opt$unipDomainBed)
