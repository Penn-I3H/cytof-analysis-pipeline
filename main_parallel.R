
library(CL2)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

dir_in <- Sys.getenv("INPUT_DIR")
dir_out <- Sys.getenv("OUTPUT_DIR")

### list channels to be used during analysis
cols <- c("CD45", "CD123", "CD19", "CD11c", "CD16",
          "CD56", "CD294", "CD14", "CD3", "CD20",
          "CD66b", "CD38", "HLA-DR", "CD45RA",
          "DNA1", "DNA2", "Event_length")

### create subdirectory structure for output
create_dirs(dir_out)

file_to_process <- args[1]
print(file_to_process)
analyze_cytof_file(file_to_process, dir_in=dir_in, dir_out=dir_out, cols=cols)




