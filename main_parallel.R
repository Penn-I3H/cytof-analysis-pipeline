
library(CL2)
library(parallel)

dir_in <- Sys.getenv("INPUT_DIR")
dir_out <- Sys.getenv("OUTPUT_DIR")
max_cores <- Sys.getenv("MAX_CORES")

files <- list.files(dir_in, pattern=".fcs")

### list channels to be used during analysis
cols <- c("CD45", "CD123", "CD19", "CD11c", "CD16",
          "CD56", "CD294", "CD14", "CD3", "CD20",
          "CD66b", "CD38", "HLA-DR", "CD45RA",
          "DNA1", "DNA2", "Event_length")

### create subdirectory structure for output
create_dirs(dir_out)

### run analysis in parallel
nc <- detectCores()
nc <- min(nc, max_cores)
cl <- makeCluster(nc)

clusterExport(cl, c("dir_in", "dir_out", "cols"))
clusterEvalQ(cl, {
  library(CL2)
})

parLapply(cl=cl, files, analyze_cytof_file,
          dir_in=dir_in, dir_out=dir_out,
          cols=cols)

stopCluster(cl)



