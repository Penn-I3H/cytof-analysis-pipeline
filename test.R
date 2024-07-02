dir_in <- Sys.getenv("INPUT_DIR")
dir_out <- Sys.getenv("OUTPUT_DIR")

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
print(dir_out)