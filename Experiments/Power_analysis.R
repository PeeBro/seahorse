# ======================================= POWER ANALYSIS FOR OUTLIER REMOVAL ON SEAHORSE DATA
# Author: Matej Stevuliak

#library(data.table)
library(magrittr)   # Needed for pipe operator %>%
library(tools)
library(dplyr)
library(reshape2)
library(reshape)
library(ggplot2)
library(tidyr)      # Needed for tidy data functions like separate
library(tidyverse)
library(readxl)     # Read xlsx file
library(here)       # Used to set the path to the package dir. on the mashine
library(knitr)      # Used for tables
library(kableExtra)

# source the functions used in analysis
source('analysis_source_functions.R')


ocr.raw <- read.csv("OCR_Raw.csv")

# Choose samples used for analysis
# I take 4 random samples from our dataset you can change the number or add specific list
# of samples in filter function as c("","", ...)
subset.samples <- ocr.raw %>%
  filter(plate_id %in% sample(unique(ocr.raw$plate_id),4, replace = F)) %>%
  group_by(plate_id) %>%
  filter(sample_id == unique(sample_id)[1]) # only for intervals 1,2,4 for Int3 chenge to ==

######### IDENTIFY OUTLIERS optimization - run this first
idfy_sinleP_outlier_ <- function(DT, cut.point, x ) {
  # Identify point outliers based on MAD of the Interval
  # IN:
  # DT        = Seahorse data produced by read_xlsx_set() function
  # cut.point = Treashold for outier removal
  # x         = Varible used, one of: "OCR", "LOCR", "ECAR", "LECAR"
  # OUT: list containing data
  # dm_r      = data with boolean col. is.out.p

  DT <- DT %>% drop_na()
  dm <- DT
  dm$x <- dm[[x]] # Create column with regressed variable

  i <- 1
  keep <- T

  size <- nrow(dm)

  while (keep) {

    # fit simple categorical regression
    fit <- lm(x ~ -1 + sample_id + Interval + Well, data = dm)
    # add column with fitted
    dm$fitted <- fitted(fit)

    # get  of fitted in interval and squared orror of that
    dm <- dm %>%
      group_by(sample_id, Interval) %>%
      mutate(int_mean = mean(fitted, na.rm = T),
             residual = (x-int_mean),
             sq_err   = (x-int_mean)^2) %>%
      ungroup()

    # SINGLE point removal procedure
    dm <- dm %>%
      group_by(sample_id, Interval) %>%
      mutate(median_sqE = median(sq_err),
             mad_sqE    = mad(sq_err, na.rm=T),
             is.out.p   = median_sqE + cut.point * mad_sqE < sq_err) %>%
      ungroup()

    n_out_p <-  nrow(dm %>% filter(is.out.p == T))
    cat(i, " Point outliares: ", n_out_p, "--", n_out_p/size*100 ,  "% \n")

    dm <- dm %>% filter(is.out.p == F)
    if( n_out_p == 0) keep <-F # Stop when no more outliars found

    i <- i+1

  }

  # complete data
  dm_r <- DT %>%
    left_join(dm,
              by = c("Measurement", "Well", "Group", "Time", "plate_id", "Protocol", "Interval"),
              suffix  = c("",".y") ) %>%
    mutate(is.out.p = replace_na(is.out.p, T)) %>%
    select(-c(contains("y"), "median_sqE", "mad_sqE", "int_mean", "x", "fitted"))

  # Print summary
  cat("Tolat single point outliars: ", nrow(filter(dm_r, is.out.p == T))/size*100, "% \n" )

  return(dm_r)
}
#########

# here we sample different numbers of wells from before selected samples.(first loop)
# sampled wells are trimmed for outliers and sd of residuals is taken,
# 10 separate samplings are made for every Biological sample (different plates
# I take random wells in every interval
removed_all <- data.frame()
for (j in seq(4,14,1)) { # seq(2,8,2) for Int3// decrease max vector size if it throws error
  removed_samplings <- data.frame()
  for (i in seq(1,10)) {
    # Subsets of wells
    subset.wells <- subset.samples %>%
      group_by(sample_id, Interval) %>%
     # filter(Interval == "Int3") %>%
      filter(Well %in% sample(unique(Well),j, replace = F))

    removed.outiers <- idfy_sinleP_outlier_(subset.wells, cut.point = 6, x = "LOCR")
    # get sd of the residuals
    removed <- removed.outiers %>%
      filter(is.out.p == F) %>%
      group_by(sample_id, Interval) %>%
      summarise(sd_residuals = sd(residual), sampling = i, well_n = j)
    # adding samplings together
    removed_samplings <- bind_rows(removed, removed_samplings)

  }
  # adding samplings of different wells numbers together
  removed_all <- bind_rows(removed_all, removed_samplings)
}

# ploting the sd of residual from separate samplings
ggplot(removed_all)+
  ggtitle("Power Analysis of Outlier Removal method")+
  geom_point(aes(well_n,sd_residuals, col = sample_id))+
  facet_grid(Interval~.)+
  ylab("SD of residuals")+
  xlab("Number of wells")+
  geom_smooth(aes(well_n,sd_residuals),method = 'loess')+
  theme_bw()

# I think that the SD doesnt change much
