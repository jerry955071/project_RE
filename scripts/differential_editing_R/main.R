# This script identifies differentially edited sites.
# 
# Dependencies:
#   Software:
#     - tabix
#   R-packages:
#     - doParallel
#     - foreach
#     - data.table
#     - magrittr
#     - rjson
#     - REDIT
# 
# Input:
#   params.json
#   
# Output:
#   pvalue.tsv
#   
# Workflow:
#   Perform statisitical test according to per run paramters

# Setup
# ==============================================================================
library(doSNOW)
library(magrittr)
library(rjson)
source("src/REDITs/REDIT_LLR.R")
# source("scripts/differential_editing/accessories/accessory.R")
# source("scripts/runREDITs.R")

# Get per run parameters
# params <- fromJSON(file = "configs/differential_editing/Egr.json")
args <- commandArgs(trailingOnly=TRUE)
params <- fromJSON(file = args[1])
message("Using parameters in config file: ", args[1])


# Setup run cluster
cl <- makeCluster(params$threads)
registerDoSNOW(cl)
message("Number of workers: ", getDoParWorkers())
# ==============================================================================


# Perform statistical tests at each target positions
# ==============================================================================
tmpdir1 <- params$tmpdir1
fileList <- list.files(tmpdir1)
dir.create(params$tmpdir2, showWarnings = FALSE, recursive = TRUE)
message("Number of tests to analyse: ", length(fileList))


# set progress bar object
pb <- txtProgressBar(max = length(fileList), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# foreach
system.time(
  foreach(i = 1:length(fileList),
          .packages = "magrittr",
          .options.snow = opts
    ) %dopar% {
    fname <- fileList[i]
    # read pre-processed basecount data
    df = read.csv(paste0(tmpdir1, fname), header = TRUE, row.names = 1)
    
    # get information from file name
    foo <- fname %>% sub(".csv$", "", .) %>% strsplit(split = "_") %>% unlist
    reg <- foo[1:c(length(foo)-2)] %>% paste(collapse = "_")
    pos <- foo[length(foo) - 1]
    subtype <- foo[length(foo)]
    non_edited <- substr(subtype, 1, 1)
    edited <- substr(subtype, 2, 2)
    
    # initiate output list
    output <- list(reg, pos, subtype)
    
    # do test for each pair
    for (pair in params$pairing) {
      foo2 <- strsplit(pair, split = ":")[[1]]
      group1 <- foo2[1]
      group2 <- foo2[2]
      # exception1: no RNA-seq record
      if (!sum(df$Tissue == group1) | !sum(df$Tissue == group2)) {
        # TODO: warning
        output %<>% append(NA)
      } else {
        # do test
        the_data <- rbind(
          df[df$Tissue == group1, c(edited, non_edited)],
          df[df$Tissue == group2, c(edited, non_edited)]
        ) %>% t %>% as.matrix
        groups <- c(
          rep(group1, sum(df$Tissue == group1)),
          rep(group2, sum(df$Tissue == group2))
        )
        output %<>% append(
          REDIT_LLR(data = the_data, groups = groups)$p.value %>% 
            sprintf("%.3f", .)
        )
      }
    }
    
    # write output
    write(
      paste(output, collapse = ","),
      file = paste0(params$tmpdir2, fname)
    )
    
    # update progress bar
    setTxtProgressBar(pb, i)
    return(1)
  }
)

close(pb)

message("Combining outputs in tmpdir: ", params$tmpdir2, " to ", params$outdir)
system(
  command = sprintf(
    "find %s | xargs cat > %s",
    params$tmpdir2,
    params$outdir
  )
)

message("Job done")

# Timing
# ==============================================================================
# Test1 using 36,984 files
# user  system elapsed
# 23.168   3.795  97.515
# >>> 0.0026 seconds/file

# Test2 using 66,092 files
# user  system elapsed 
# 48.667   7.223 190.846 
# >>> 0.0028875809477698 seconds/file


# # ABANDONED
# # Perform statistical tests at each target positions
# # ==============================================================================
# target_positions <- fread(
#   params$target_positions,
#   col.names = c("Region", "Position")
# )

# # partitioning the `target_position` data.table
# add_pos <- function(df) {
#   df$pos <- paste(df$Region, df$Position, sep = ":"); df
# }
# partitions <- split(
#   target_positions,
#   sort(1:nrow(target_positions) %% params$threads),
# ) %>% lapply(., add_pos) # -> list of data.table

# # import sample reditables as data.table
# reditables <- params$samples %>% 
#                 unlist %>% 
#                   lapply(., fread, col.names = params$redithead)

# # paste column `Region` and `Position` as new column `pos`
# reditables %<>% lapply(., add_pos)

# # naming
# names(reditables) <- basename(params$samples %>% unlist) 

# # partitioning
# # reditables$xylem %<>% lapply(
# #   .,
# #   function(df, p) {
# #     lapply(p, 
# #       function(p, df) {
# #         return(df[df$RP %in% p$RP])
# #       },
# #       df = df
# #     )
# #  },
# #  p = partitions
# # )


# # for each partition
# N <- 1
# foreach(N = 1:params$threads) %dopar% {
#   # position
#   pN <- partitions[[N]]
  
#   # create output matrix
#   output <- matrix(
#     nrow = nrow(pN),
#     ncol = length(params$pairing) + ncol(pN),
#     dimnames = list(NULL, c(colnames(pN), params$pairing))
#   ) %>% as.data.table
#   output[,1:3] <- pN
  
#   # for-for-loop
#   for (pos in pN$pos) {
#     sapply(reditables, function(x, pos){
#       return(x[x$pos == pos, ]),
#       pos = pos
#     }) %>% rbind
    
#     for (pair in params$pairing) {
#       # 
#       group_name1 <- strsplit(pair, ":")[[1]][1]
#       group_name2 <- strsplit(pair, ":")[[1]][2]
#       sample_names1 <- names(reditables[[group_name1]])
#       sample_names2 <- names(reditables[[group_name2]])
      
#       # REDIT inputs
#       the_data <- matrix(
#         nrow = 2,
#         ncol = length(sample_names1) + length(sample_names2)
#       )
#       groups <- c(
#         rep(group_name1, length(sample_names1)),
#         rep(group_name2, length(sample_names2))
#       )
#     }
     
      
#   }
#   return(output)
# }

# ==============================================================================