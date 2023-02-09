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

library(doParallel)
library(data.table)
library(magrittr)
library(rjson)
source("scripts/runREDITs.R")

# Get per run parameters
params <- fromJSON(file = "configs/R/differential_editing/Ptr.json")

# Perform statistical tests on each target position
target_positions <- fread(
  file = params$target_positions,
  col.names = c("Region", "Position")
)

foreach() %dopar% {}


