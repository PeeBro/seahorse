# ======================================= FUNCTIONS FOR SEAHORSE PIPELINE
# Author: Matej Stevuliak


idfy_sinleP_outlier <- function(DT, cut.point, x ) {
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
    fit <- lm(x ~ -1 + Interval + Well, data = dm)
    # add column with fitted
    dm$fitted <- fitted(fit)

    # get mean of fitted in interval and squared orror of that
    dm <- dm %>%
      group_by(sample_id, Interval) %>%
      mutate(int_mean = mean(fitted, na.rm = T),
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
    select(-c(contains("y"), "median_sqE", "mad_sqE", "sq_err", "int_mean", "x", "fitted"))

  # Print summary
  cat("Tolat single point outliars: ", nrow(filter(dm_r, is.out.p == T))/size*100, "% \n" )

  return(dm_r)
}
