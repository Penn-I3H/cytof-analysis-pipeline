
library(CL2)
library(parallel)

dir_in <- Sys.getenv("INPUT_DIR")
dir_out <- Sys.getenv("OUTPUT_DIR")
use_parallel <- Sys.getenv("USE_PARALLEL")

files <- list.files(dir_in, pattern=".fcs")

### list channels to be used during analysis
cols <- c("CD45", "CD123", "CD19", "CD11c", "CD16",
          "CD56", "CD294", "CD14", "CD3", "CD20",
          "CD66b", "CD38", "HLA-DR", "CD45RA",
          "DNA1", "DNA2", "Event_length")

### create subdirectory structure for output
create_dirs(dir_out)

if (use_parallel == 1) {
  ### run analysis in parallel
  nc <- detectCores()

  message(paste0("Detected ", nc, " cores."))

  log_file <- paste0(dir_out, "/log.txt")
  cl <- makeCluster(nc, outfile=log_file)

  clusterExport(cl, c("dir_in", "dir_out", "cols"))
  clusterEvalQ(cl, {
    library(CL2)
  })

  parLapply(cl=cl, files, analyze_cytof_file,
            dir_in=dir_in, dir_out=dir_out,
            cols=cols)

  stopCluster(cl)
} else {
  ### run serial analysis
  lapply(files, analyze_cytof_file,
         dir_in=dir_in, dir_out=dir_out,
         cols=cols)
}




