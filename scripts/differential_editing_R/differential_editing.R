"
This script identifies differentially edited sites.

Dependencies:
  Software:
    - tabix
  R-packages:
    - doParallel
    - foreach
    - data.table
    - magrittr
    - rjson
    - REDIT

Input:
  params.json
  
Output:
  pvalue.tsv
  
Workflow:
  Perform statisitical test according to per run paramters
"
# Setup
# ==============================================================================
library(doParallel)
library(data.table)
library(magrittr)
library(rjson)
source("src/REDITs/REDIT_LLR.R")
source("scripts/differential_editing/accessories/accessory.R")
# source("scripts/runREDITs.R")

# Get per run parameters
params <- fromJSON(file = "configs/R/differential_editing/Ptr.json")

# Setup run cluster
cl <- makeCluster(params$threads)
registerDoParallel(cl)
getDoParWorkers() # check parallel workers
# ==============================================================================


# Perform statistical tests at each target positions
# ==============================================================================



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