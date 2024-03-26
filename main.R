
library(CL2)

dir_in <- Sys.getenv("INPUT_DIR")
dir_out <- Sys.getenv("OUTPUT_DIR")

files <- list.files(dir_in, pattern=".fcs")

### list channels to be used during analysis
cols <- c("CD45", "CD123", "CD19", "CD11c", "CD16",
          "CD56", "CD294", "CD14", "CD3", "CD20",
          "CD66b", "CD38", "HLA-DR", "CD45RA",
          "DNA1", "DNA2", "length")

### create subdirectory structure for output
create_dirs(dir_out)

### for now running sequentially
### later we can parallelize
lapply(files, analyze_cytof_file,
       dir_in=dir_in, dir_out=dir_out,
       defs=defs_major, defs_myel=defs_myel,
       defs_tcell=defs_tcell, cols=cols)




