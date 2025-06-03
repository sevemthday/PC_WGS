library(data.table)
library(optparse)

# variable
option_list <- list(
  make_option(c("-I", "--input"),  action = "store", type = "character", default = NULL, help = "input file pattern"),
  make_option(c("-O", "--output"), action = "store", type = "character", default = NULL, help = "full path to output file"),
  make_option(c("--path"),         action = "store", type = "character", default = ".",  help = "full path to input dir"),
  make_option(c("--subdirectory"), action = "store", type = "logical",   default = FALSE,help = "list input files in subdirectory"),
  make_option(c("--addFileName"),  action = "store", type = "logical",   default = FALSE,help = "add file name to datatable"),
  make_option(c("--sep"),          action = "store", type = "character", default = NULL, help = "sep for file name strsplit split"),
  make_option(c("--dedup"),        action = "store_true",                default = FALSE,help = "remove duplicate rows"),
  make_option(c("--col_char"),     action = "store_true",                default = FALSE,help = "read data as character"),
  make_option(c("--sortColumn"),   action = "store", type = "character", default = NULL, help = "columns, which is to sort row. sep by colon; eg: column1:column2:column3")
)
opt <- parse_args(OptionParser(option_list = option_list))

# merge table
filelist <- list.files(path = opt$path, pattern = opt$input, full.names = T, recursive = opt$subdirectory)
dtlist <- list()
for (i in seq(filelist)) {
  
  if(opt$col_char){
    dt <- fread(filelist[i], colClasses = "character")
  }else{
    dt <- fread(filelist[i])
  }
  if(opt$addFileName){
    dt$Sample_name <- strsplit(filelist[i], split = opt$sep)[[1]][1]
  }
  dtlist[[i]] <- dt
  
}
dtlist <- rbindlist(dtlist)

if(opt$dedup) {
  dtlist <- unique(dtlist)
}

if(!is.null(opt$sortColumn)) {
  keycol <- unlist(strsplit(opt$sortColumn, split = ":"))
  dtlist <- setorderv(dtlist, keycol)
}

write.table(dtlist, file = opt$output, quote = F, sep = "\t", row.names = F)