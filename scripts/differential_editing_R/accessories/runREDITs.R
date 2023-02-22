# Usage:
# Rscripts --vanilla threads edsites.txt sample1.txt.gz sample2.txt.gz ...

# load package
# library(doParallel)
# library(magrittr)


# cmd line args
# args = commandArgs(trailingOnly=TRUE)
# threads = args[2]
# pEdsites = args[3]
# pSamples = args[4:]

# for debugging
# threads = 60
# pEdsites = "shout/edsites/ptr.txt"
# pSamples = c(
#   "pyout/archive-20230104/ptr-SAMN04535869-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04535870-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04535936-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04535937-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04535946-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04535958-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04535959-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04535960-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04536034-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04536035-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04536036-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04536096-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04536100-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04536101-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04536102-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04536244-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04536245-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04536246-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04569083-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04569084-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04569085-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04569086-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04939121-sr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04939122-sr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04939123-sr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04939124-lr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04939125-lr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04939126-lr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04939127-pr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04939128-pr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04939129-pr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04939134-xr.txt.gz",
#   "pyout/archive-20230104/ptr-SAMN04939135-xr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp1-xr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp2-xr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp3-xr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp4-xr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp5-xr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp6-xr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp1-pr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp2-pr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp3-pr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp4-pr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp5-pr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp6-pr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp1-lr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp2-lr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp3-lr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp4-lr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp5-lr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp6-lr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp1-sr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp2-sr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp3-sr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp4-sr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp5-sr.txt.gz",
#   "pyout/archive-20230104/ptr-tsp6-sr.txt.gz"
# )
# 
# reditHeader <- c(
#   "Region",
#   "Position",
#   "Reference",
#   "Strand",
#   "Coverage",
#   "MeanQ",
#   "BaseCount",
#   "AllSubs",
#   "Frequency",
#   "gCoverage",
#   "gMeanQ",
#   "gBaseCount",
#   "gAllSubs",
#   "gFrequency"
# )

# # vars
# pREDIT_LLR <- "src/REDITs/REDIT_LLR.R"
# pREDIT_regression <- "src/REDITs/REDIT_regression.R"

# functions
source("src/REDITs/REDIT_LLR.R")
# source(pREDIT_regression)

# `tabix` caller
tabix <- function(fname, chr, start, end = NA) {
  end <- ifelse(is.na(end),  start, end)
  cmd <- sprintf("tabix %s %s:%s-%s", fname, chr, start, end)
  return(system(cmd, intern = TRUE))
}

# sample name --> group name
s2g <- function(sname) {
  list(
    "x" = "xylem",
    "p" = "phloem",
    "l" = "leaf",
    "s" = "shoot"
  ) %>%
    extract(substring(sname, nchar(sname) - 8, nchar(sname) - 8))
}

# get Data
getData <-
  function(reader, paths, colnames, chr, start, end = NA, para = FALSE) {
    # Read corresponding lines from `sample` using reader
    # output: data.frame
    if (para) {
      lines <- foreach(p=paths) %dopar% reader(p, chr, start, end)
    } else {
      lines <- lapply(paths, reader, chr, start, end)  
    }
    names(lines) <- paths
    lines <- unlist(lines)
    data <- read.table(
      text = lines,
      sep = "\t",
      header = FALSE,
      col.names = colnames,
      row.names = names(lines)
    )
    # add columns
    data["group"] <-  lapply(row.names(data), s2g) %>% unlist
    data[c("A", "C", "G", "T")] <-
      data$BaseCount %>%
      gregexpr("[0-9]+", .) %>%
      regmatches(data$BaseCount, .) %>%
      lapply(strtoi) %>%
      unlist %>%
      matrix(ncol = 4,
             byrow = TRUE,
             dimnames = list(NULL, c("A", "C", "G", "T")))
    return(data)
  }


# perform pairwise test using REDIT_LLR
pairwiseTest <- function(data, group1, group2) {
  ggdata <- rbind(data[which(data["group"] == group1),],
                  data[which(data["group"] == group2),])
  subtype <- unique(ggdata$AllSubs) %>%
    extract(., which(. != "-"))
  
  if (length(subtype) == 0) {
    return(list(p.value = NA))
  }
  else {
    type_group <- ggdata[["group"]]
    if (length(unique(type_group)) != 2) {
      return(list(p.value = NA))
    }
    else {
      type_data <- t(ggdata[, strsplit(subtype, "")[[1]]])
      colnames(type_data) <- type_group
      return(REDIT_LLR(type_data, type_group))
    }
    
  }
}
# 
# ########
# 
# # testing
# # configure doParalle workers
# cl <- makeCluster(threads)
# registerDoParallel(cl)
# getDoParWorkers()
# 
# edSites <- read.table(pEdsites)
# 
# allCombn <- list(
#   xp = c("xylem", "phloem"),
#   xl = c("xylem", "leaf"),
#   xs = c("xylem", "shoot"),
#   pl = c("phloem", "leaf"),
#   ps = c("phloem", "shoot"),
#   ls = c("leaf", "shoot")
# )
# 
# 
# # lapply + foreach(row)
# # for each generic position
# testn <- 345950
# pout <- foreach(nr = 1:testn,
#                 .packages = "magrittr",
#                 .combine = "rbind") %dopar% {
#                   # Output:
#                   #   A data.frame with dimension = (subtypes) * (comb + 1)
#                   #   example:
#                   #
#                   
#                   # read data from all files in `pSamples`
#                   data <- getData(tabix,
#                                   pSamples,
#                                   reditHeader,
#                                   edSites[nr, 1],
#                                   edSites[nr, 2]
#                                   )
#                   
#                   # exception 1: code 34
#                   # check nrow(data) > 1
#                   if (nrow(data) <= 1) {
#                     return(34)
#                   }
#                   
#                   # check total types of substitution at this gene position
#                   allSubs <- data["AllSubs"] %>%
#                     unique %>%
#                     lapply(strsplit, split = " ") %>%
#                     unlist %>%
#                     unique
#                   
#                   out <- data.frame(matrix(
#                     nrow = sum(allSubs != "-"),
#                     ncol = length(allCombn) + 1
#                   ))
#                   
#                   if (sum(allSubs != "-") > 1) {
#                     return(22)
#                   }
#                   
#                   pvals <- lapply(allCombn, function(x) {
#                     pairwiseTest(data, x[1], x[2])$p.value
#                   }) %>% unlist
#                   
#                 }
# 
# pout
# #
# 
# 
# #
# library(data.table)
# df <- as.data.table(pout)
# df <- cbind(edSites[1:testn,], df)
# # tmp <-
# #   df[df$xp < 0.05 & df$xl < 0.05 & df$xs < 0.05, 1:2] %>% na.exclude
# # 
# # for (r in 1:nrow(tmp)) {
# #   getData(tabix,
# #           pSamples,
# #           reditHeader,
# #           tmp[r, 1],
# #           tmp[r, 2])[, c(
# #             "Coverage",
# #             "AllSubs",
# #             "BaseCount",
# #             "Frequency",
# #             "gBaseCount",
# #             "gAllSubs",
# #             "gFrequency"
# #           )] %>%
# #     cbind(s2g(rownames(.)) %>% unlist) %>%
# #     View()
# #   readline(prompt = "Press [enter] to continue; [esc] to quit")
# # }
# # getData(tabix,
# #         pSamples,
# #         reditHeader,
# #         "Chr01",
# #         129111)[, c(
# #           "Coverage",
# #           "AllSubs",
# #           "BaseCount",
# #           "Frequency",
# #           "gBaseCount",
# #           "gAllSubs",
# #           "gFrequency"
# #         )] %>% View()
# # 
# # 
# # print(dim(df))
# # print(dim(edSites))
# # final <- cbind(edSites, df)
# write.table(
#   df,
#   file = "pval.tsv",
#   quote = FALSE,
#   sep = "\t",
#   row.names = FALSE,
#   na = "None"
# )
# 
# # dimnames(final_output)[[2]] <- names(allCombn)
# # final_output
# 
# 
# #    user  system elapsed
# #   3.316   0.402  16.349
# 
# 
# # # NOTE:
# # # 1. Testing `doParallel`
# # cl <- makeCluster(60)
# # registerDoParallel(cl)
# # length(samples)
# # # [1] 57
# 
# # tabix
# # # function(fname, chr, start, end = NA) {
# # #     end <- ifelse(is.na(end),  start, end)
# # #     cmd <- sprintf("tabix %s %s:%s-%s", fname, chr, start, end)
# # #     return(system(cmd, intern = TRUE))
# # # }
# # # <bytecode: 0x55753ab17860>
# 
# # system.time(
# #     {
# #         tmp <- lapply(samples, tabix, "Chr01", 13951)
# #     }
# # )[[3]]
# # # [1] 1.219
# # # There were 50 or more warnings (use warnings() to see the first 50)
# # system.time(
# #     {
# #         tmp <- foreach(s=samples) %dopar% {tabix(s, "Chr01", 13951)}
# #     }
# # )[[3]]
# # # [1] 0.118
# 
# # # COMMENT:
# # # The self-defined function `tabix` runs pretty slow (probabily due to IO latency). Using `doParallel` we boost the execution time to around 10 times faster using 60 process for 57 samples
# 
# 
# # # for-loop
# # system.time(
# #     {
# #         pvals <- list()
# #         for (nr in 1:200) {
# #             data <- getData(
# #                 tabix,
# #                 pSamples,
# #                 reditHeader,
# #                 edSites[nr, 1],
# #                 edSites[nr, 2]
# #             )
# #             out <- pairwiseTest(data, "xylem", "leaf")
# #             pval <- out %>% unlist(recursive = FALSE) %>% extract2("p.value")
# #             pvals %<>% append(pval)
# #         }
# #     }
# # )
# # # ^C
# # # Timing stopped at: 18.01 62.11 79.88
# 
# 
# # # foreach at tabix
# # system.time(
# #     {
# #         pvals <- foreach(nr=1:200,  .packages = "magrittr") %do% {
# #             data <- getDataParallel(
# #                 tabix,
# #                 pSamples,
# #                 reditHeader,
# #                 edSites[nr, 1],
# #                 edSites[nr, 2]
# #             )
# #             out <- pairwiseTest(data, "xylem", "leaf")
# #             pval <- out %>% unlist(recursive = FALSE) %>% extract2("p.value")
# #             pval
# #         }
# #     }
# # )
# # #    user  system elapsed
# # #  14.144   1.353  38.629
# 
# # # Un-used
# # getDataParallel <- function(reader, samples, colnames, chr, start, end = NA) {
# #     lines <- foreach(s=samples) %dopar% {reader(s, chr, start)}
# #     # lines <- lapply(samples, reader, chr, start)
# #     names(lines) <- samples
# #     data <- read.table(
# #         text = unlist(lines),
# #         sep = "\t",
# #         header = FALSE,
# #         col.names = colnames,
# #         row.names = names(unlist(lines))
# #     )
# #     # add columns
# #     data["group"] <-  lapply(row.names(data), s2g) %>% unlist
# #     data[c("A", "C", "G", "T")] <-
# #         data$BaseCount %>%
# #         gregexpr("[0-9]+", .) %>%
# #         regmatches(data$BaseCount, .) %>%
# #         lapply(strtoi) %>%
# #         unlist %>%
# #         matrix(
# #             ncol = 4,
# #             byrow = TRUE,
# #             dimnames = list(NULL, c("A", "C", "G", "T"))
# #         )
# #     return(data)
# # }
# 
# 
# # # Un-used2: Nested `foreach`
# # system.time(
# #     {
# #         pout <- foreach(
# #             nr=1:100,
# #             .packages = "magrittr",
# #             .combine = "rbind"
# #         ) %:% foreach(
# #             g=allCombn,
# #             .packages = "magrittr",
# #             .combine = "c"
# #         ) %dopar% {
# #             data <- getData(
# #                 tabix,
# #                 pSamples,
# #                 reditHeader,
# #                 edSites[nr, 1],
# #                 edSites[nr, 2]
# #             )
# #             if (length(unique(data["AllSubs"])) > 2) {
# #                 stop(nr)
# #                 NA
# #             } else {
# #                 pairwiseTest(data, g[1], g[2])$p.value
# #             }
# #             # tryCatch(
# #             #     { pairwiseTest(data, g[1], g[2])$p.value },
# #             #     error = function(e) {
# #             #         print(nr)
# #             #         print(g)
# #             #         return(data)
# #             #     }
# #             # )
# #         }
# #     }
# # )
# # #    user  system elapsed
# # #   3.953   0.382  38.558
