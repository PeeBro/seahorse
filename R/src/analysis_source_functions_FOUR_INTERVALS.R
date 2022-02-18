# ======================================= FUNCTIONS USED IN ANALYSIS
# Author: Matej Stevuliak and Stine Østergaard

# ------------------- BACKGROUND CORRECTION

options(dplyr.summarise.inform = FALSE)

### Melt functions to make pretty data frames ###
melt_bioenergetics_sample <- function(d) {
  bio <- melt(d$bioenergetics)
  bio.se <- melt(d$standard.errors)
  colnames(bio.se)[3] <- "se"
  bio <- merge(bio,bio.se)
  colnames(bio)[3] <- "mean"
  return(bio)
}

melt_bioenergetics_replicate <- function(d) {
  bio <- melt(d$bioenergetics)
  bio.se <- melt(d$standard.errors)
  colnames(bio.se)[3] <- "se"
  bio <- merge(bio,bio.se)
  colnames(bio)[3] <- "mean"
  bio <- bio %>%  mutate(Date = substr(Sample, start = nchar(Sample)-15, stop = nchar(Sample)),
                                                 Group = substr(Sample, start = 1, stop = nchar(Sample)-19))
  return(bio)
}

melt_bioenergetics_well <- function(d) {
  bio <- melt(d$bioenergetics)
  bio.se <- melt(d$standard.errors)
  colnames(bio.se)[4] <- "se"
  bio <- merge(bio,bio.se)
  colnames(bio)[4] <- "mean"
  bio <- bio %>%  mutate(Date = substr(Sample, start = nchar(Sample)-15, stop = nchar(Sample)),
                         Group = substr(Sample, start = 1, stop = nchar(Sample)-19))
  return(bio)
}

### Plot functions to make plotting easier ###

plot1 <- function(well, replicate, sample){
  well %>%
    filter(variable != "Other") %>%
    ggplot(aes(x = "", y = mean))+
    geom_beeswarm(aes(color = Date, shape = Group, group = Group), cex = 3, dodge.width = 0.8)+
    geom_point(data = replicate, mapping = aes(x = "", y = mean, group = Group, fill = Date), 
               color = "black", shape = 21, size = 5, 
               position = position_dodge(width = 0.8))+
    geom_point(data = sample, mapping = aes(x = "", y = mean, group = Sample), 
               fill = "black", color = "red", shape = 21, size = 1.5,
               position = position_dodge(width = 0.8))+
    geom_errorbar(data = sample, aes(ymin=mean-se, ymax=mean+se, group = Sample), width=0.5, size = 0.8,
                  position=position_dodge(0.8))+
    scale_shape_manual(values= c(15,8, 17,3))+
    scale_color_brewer(palette = "Dark2")+
    scale_fill_brewer(palette = "Dark2")+
    facet_wrap(~variable, scales = "free")+
    xlab("Bio-Energetics")+
    theme_bw()
  
}

plot2 <- function(replicate, sample) {
  replicate %>%
    filter(variable != "Other") %>%
    ggplot(aes(x = variable, y = mean, group = Group))+
    geom_point(aes(fill = Group, shape = Date), size = 3, position = position_dodge(width = 0.6))+
    geom_point(data = sample, mapping = aes(x = variable, y = mean, group = Sample), color = "black", size = 1.5,
               position = position_dodge(width = 0.6))+
    geom_errorbar(data = sample, aes(ymin=mean-se, ymax=mean+se, group = Sample), width=0.5, size = 0.8,
                  position=position_dodge(0.6))+
    scale_shape_manual(values= c(21,22))+
    scale_color_brewer(palette = "Dark2")+
    scale_fill_brewer(palette = "Dark2")+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    xlab("Bio-Energetics")+
    theme_bw()
}

plot3 <- function(well) {
  well %>%
    filter(variable != "Other") %>%
    ggplot(aes(x = variable, y = mean))+
    geom_boxplot(aes(fill = Group, color = Group), alpha = 0.5, position = position_dodge(width = 0.8),)+
    geom_point(aes(color = Group), position = position_dodge(width = 0.8), size = 1)+
    scale_color_brewer(palette = "Dark2")+
    scale_fill_brewer(palette = "Dark2")+
    xlab("Bio-Energetics")+
    theme_bw()
}



correct_background_OCR <- function(d){
  # correct for outliars in background measurements and remove background
  # takes wells marked as Background and for each plate and time substracts median of those from OCR measurements
  # Each measurement has its own background (Mean_to_substr)

  background <- d %>%
    filter(Group == "Background") %>%
    group_by(Measurement) %>%
    mutate( Z_mod = abs((OCR - median(OCR))/mad(OCR)) )

  out        <- background %>% filter(Z_mod > 2.77)
  background <- background %>%
    filter(Z_mod < 2.77) %>%
    mutate(Mean_to_substr = mean(OCR))
  background <- background[!duplicated(background[,"Mean_to_substr"]),] # reduce the data before joining them

  d <- d %>%
    left_join(background, by = "Measurement", suffix  = c("",".y") ) %>%
    mutate(OCR = OCR - Mean_to_substr ) %>%
    select(-c(contains("y"),"Z_mod", "Mean_to_substr"))

  return(list(data = d, removed = out))
  rm(background,d)
}

correct_background_PER <- function(d){
# PER Normalizarion
  background <- d %>%
    filter(Group == "Background") %>%
    group_by(Measurement) %>%
    mutate( Z_mod = abs((PER - median(PER))/mad(PER)))


  out        <- background %>% filter(Z_mod > 3.26)
  background <- background %>%
    filter(Z_mod < 3.26) %>%
    mutate(Mean_to_substr = mean(PER)) # take mean to substract

  background <- background[!duplicated(background[,"Mean_to_substr"]),] # reduce the data before joining them

  d <- d %>%
    left_join(background, by = "Measurement", suffix  = c("",".y") ) %>%
    mutate(PER = PER - Mean_to_substr ) %>%
    select(-c(contains("y"),"Z_mod", "Mean_to_substr"))

  return(list(data = d, removed = out))
  rm(background,d)
}

get_intervals <- function(x, y){
  # Identifies Interval from measurement and accumulated sum of number of measurements in each interval. Used with mapply
  # INPUT x = Measurement, y = accumulated sum of number of measurements in each interval
  # OUTPUT = Interval
  interval = 1
  for(i in y){ # Iterate through all the intervals. i is now the value of the max measurement of the interval of each round of this for-loop.
    
    if(x <= i){
      return(sprintf("Int%1.0f", interval)) #If The measurement is smaller than the largest measurement number of the interval, assign that measurement to that interval
      
    }else{
      interval = interval + 1 # Otherwise go to next interval. (simplified explanation)
    }
  }
  return("Other")
}

normalize_based_on_cell_counts <- function(filename, d) {
  
  # Read cell counts (suppressMessages removes an irrelevant R message in the final output)
  cell_counts <- suppressMessages(read_xlsx(filename, 
            sheet = "Assay Configuration",
            range = "B84:N92",
            col_names = TRUE),  classes = "message")
  
  # Rename columns so the numbering is identical to that of the well names
  colnames(cell_counts)[2:10] <- c("01", "02", "03", "04", "05", "06", "07", "08", "09")
  
  # Convert the matrix into three columns: One with the well column number, one with the well letter and one with the number of cells in that well
  cell_counts <- cell_counts %>% pivot_longer(cols = !"...1", names_to = "well_col", values_to = "cells")
  colnames(cell_counts)[1] <- "well_row"
  
  # Make a new column with cells/1000 and make a new column with the well name (combine letter and number)
  cell_counts <- cell_counts %>% mutate(cells_1000 = cells/1000,
                                        Well = paste(well_row, well_col, sep = ""))
  
  # Remove unwanted columns
  cell_counts <- cell_counts[4:5]
  
  # Add the cell/1000 count values too the d table, and calculate the normalized OCR, ECAR and PER
  d <- merge(d, cell_counts)
  
  d <- d %>% mutate(OCR_norm = OCR/cells_1000,
                    PER_norm = PER/cells_1000,
                    ECAR_norm = ECAR/cells_1000) 
  
  # Remove old OCR, ECAR and PER columns, and rename the normalized columns to OCR, ECAR, PER
  d <- d %>% select(!c(OCR, ECAR, PER))
  colnames(d)[7:9] <- c("OCR", "PER", "ECAR")
      
  
  return(d)
}

report_removed_wells <- function(out) {
  out_wells_pretty <- out %>% group_by(Well, plate_id) %>% arrange(Measurement) %>% summarise(Measurement = toString(Measurement),
                                                                                                n_out = n())
  df <- data_frame(
    Plate = out_wells_pretty$plate_id,
    N_out = out_wells_pretty$n_out,
    Wells = out_wells_pretty$Well,
    Measurement = out_wells_pretty$Measurement)
  
  return(df)
}

# --------------------------------------------------------------------------------------CONSTRUCT DATA


read_xlsx_set <- function(path_, pattern_){
  # returns a combined dataframe from all .xlsx files in folder (path_).
  # can specify pattern to distinguish "group1" "group2"  treatment as: pattern_ = "*group1/group2.xlsx"
  
  # list of filenames
  files    <- list.files(path=path_, pattern=pattern_, full.names=TRUE, recursive=FALSE)
  merged_d <- data_frame()
  Hg_out   <- data_frame()
  PER_background_out <- data_frame()
  OCR_background_out  <- data_frame()
  Empty_out <- data_frame()
  No_cells_out <- data_frame()
  Bad_cells_out <- data_frame()
  
  # for every file create a data frame and merge into one
  for (x in files) {
    # --------- READ DATA --------- #
    # Read rates
    d <- read_xlsx(x,  sheet = "Rate")

    
    # Read assay configuration
    assay_name <- toString( read_xlsx(x,  
                                      sheet = "Assay Configuration", 
                                      range = "B3:B4",
                                      col_names = " ")[2,1])
    exper_date <- substr(toString(read_xlsx(x,  
                                            sheet = "Assay Configuration", 
                                            range = "B23:B24",
                                            col_names = " ")[2,1]),1,10)
    
    exper_time <- toString(read_xlsx(x,  
                                            sheet = "Assay Configuration", 
                                            range = "B24:B24",
                                            col_names = " "))
    

    # Check if PER column is empty. If yes stop the analysis with a message
    if (nrow(d) == sum(d$PER == 0)) {
      stop(sprintf("File with assay name '%s' has only zeroes in the PER column. Analysis stopped. Please check input files before trying again.",assay_name ))
      
    }
    
    
        
    #### Add plateID column
    d$plate_id <- rep(paste0(assay_name, " ", exper_date), times = nrow(d)) #? Why is this needed as repeats?
    
    
    # Read number of intervals + their length
    intervals <- read_xlsx(x,
                           sheet = "Assay Configuration", 
                           range = "C47:F47", col_names = c("1","2","3","4"))
    
    
    length_intervals <- as.integer(sub(".*: ", "", intervals))
    acc_length_intervals <- as.data.frame(cumsum(length_intervals))
    number_intervals <- length(length_intervals)
    
    message(sprintf("Name of assay: %s", assay_name), 
            sprintf("\nNumber of intervals : %1.0f",number_intervals),
            "\nLength of intervals : ", paste(length_intervals, collapse = ","))

    
    
    # Get measurement length (How many lines does each measurement have?. Used to load the data in "Read mmHg from raw data" )
    i = 1400
    while(TRUE){
      raw <- read_xlsx(x,  sheet = "Raw", range = sprintf("A1:A%1.0f",i))
      if(length(unique(raw$Measurement)) > 1) {
        measurement_length <- max(table(raw$Measurement))
        break
      } else {
        i = i + 500
      }
      
    }
    
    
    #Read mmHg from raw data
    raw <- read_xlsx(x,  sheet = "Raw", range = sprintf("A1:I%1.0f",measurement_length*length_intervals[1]+1)) # reading only columns 1-9 and first 3 measurements
    # find wells where mmHg is out of range 140 - 160
    out_Hg <- raw %>%
      group_by(Measurement) %>%
      filter(Tick == min(Tick)) %>% # Extracts data with the smallest Tick value in each measurement
      ungroup() %>% 
      group_by(Well) %>%
      summarise(average_mmHg = mean(`O2 (mmHg)`), # Calculates average Hg of the measurements (first tick of each) in interval 1 for each well
                out          = ifelse(average_mmHg > 160 | average_mmHg < 140, T, F )) %>%
      filter(out == T)
    
    # remove Unassigned wells
    d <- d %>%filter(Group != "Unassigned")
    
    # Check for missing PER values in background measurements and correct them
    background <- d %>% filter(Group == "Background")
    zero_PER <- nrow(background[background$PER == 0,]) > 0 # Looks in background table and check if any PER values are 0. If yes the value of zero_PER is "TRUE". If no zeroes are present, then zero_PER = "FALSE"
    zero_ECAR <- nrow(background[background$ECAR == 0,]) > 0 # Repeat for ECAR. If ECAR has no zeroes we can calculate PER from ECAR
    if(zero_PER == TRUE & zero_ECAR == FALSE) {
      conversion_factor <- d %>% filter(Group != "Background") %>% summarize(conversion = round(PER/ECAR, 4))
      conversion_factor <- unique(conversion_factor)
      if(nrow(conversion_factor) == 1) {
        conversion_factor <- as.double(conversion_factor)
        d <- d %>% mutate(PER = ifelse(PER == 0, ECAR*conversion_factor,PER))
        message(sprintf("Some or all PER background values are zero in file '%s'. Missing PER background values are calculated by: ECAR * %f (ratio of PER and ECAR)", assay_name,conversion_factor))
        
      } else {
        message(sprintf("Some or all PER background values are zero in file '%s'. Missing PER values could not be calculated, due to inconsistant ratios between ECAR and PER",assay_name))
      }
    }
    
    
    # Remove wells with OCR less than 10 in the first interval and which is NOT a background measurement
    
    out <- d %>% filter(Measurement <= acc_length_intervals[1,1]) %>% filter(!Group == "Background") %>% filter(OCR <= 10)
    less_than_10_OCR <- report_removed_wells(out) # Make table with removed wells
    d <- d %>% filter(!Well %in% out$Well) # Remove bad wells from the data
    Bad_cells_out <- rbind(Bad_cells_out, less_than_10_OCR) # Merge table with all other excel files seen so far
    
    
    
    #### remove Background from OCR  and PER ####
    is_norm_OCR  <- mean(filter(d, Group == "Background")$OCR) == 0
    is_norm_PER <- mean(filter(d, Group == "Background")$PER) == 0
    if (is_norm_OCR & is_norm_PER) {
      d <- d %>% filter(Group != "Background")
    } else {
      corr <- correct_background_PER(d)
      d    <- corr$data
      out  <- corr$removed
      # to report which measurements were removed
      removed_PER_background <- report_removed_wells(out)
      PER_background_out <- rbind(PER_background_out, removed_PER_background)
      #OCR
      corr <- correct_background_OCR(d)
      d    <- corr$data
      out  <- corr$removed
      # to report which measurements were removed      
      removed_OCR_background <- report_removed_wells(out)
      OCR_background_out <- rbind(OCR_background_out, removed_OCR_background)
      
      d <- d %>% filter(Group != "Background") %>% 
        drop_na()
    }
    
    #### Normalize based on cell counts (pr. 1000 cells) and remove wells where no cells are measured
    d <- normalize_based_on_cell_counts(x,d)

    out <- d %>% filter(cells_1000 == 0)
    no_cells_measured <- report_removed_wells(out)
    d <- d %>% filter(!cells_1000 == 0)
    No_cells_out<- rbind(No_cells_out, no_cells_measured)    
    
    
    # Remove wells where OCR or PER is 0
    out <- d %>% filter(OCR == 0 | PER == 0)
    zero_OCR_PER <- report_removed_wells(out)
    d <- d %>%filter(!(OCR == 0 | PER == 0))
    Empty_out<- rbind(Empty_out, zero_OCR_PER)
    
    # Filter out whole wells where average of first ticks from first three measurements
    # are exceeding interval 140 - 160 mmHg
    # Filtering has to be done after background norm. to avoid removing controll backgrounds
    d <- d %>% filter(! Well %in% out_Hg$Well)
    # make line for report matrix
    removed_Hg <- data_frame(
        Plate = unique(d$plate_id),
        N_out = nrow(out_Hg),
        Wells = paste0(out_Hg$Well, collapse = ", "))
    Hg_out <- rbind(Hg_out, removed_Hg)
    rm(raw, out_Hg)
    # ----------------
    
    ### Add Protocol column substring of Group column before "-" character
    d$Protocol <- substr(sub("-.*","", d$Group), start = 1, stop = 4)
    ### add project column
    d$Project  <- sub("#.*", "", sub(".*-","", d$Group))
    ### Add interval column
    d$Interval <- mapply(get_intervals, d$Measurement, y = acc_length_intervals)
    
    ### Add column with time on experiment scale
    d <- d %>% mutate(Time = ifelse(Interval == "Int1" | Interval == "Int2", Measurement,
                                    ifelse(Interval == "Int3" | Interval == "Int4", Measurement + 3, Measurement + 6)))
    
    # the entries with blanks should be filtered, Must be here !
    d <- d %>%filter(!is.na(Protocol))
    
    ### Add log OCR,  log PER
    d$LOCR  <- log(d$OCR)
    d$LPER <- log(d$PER)
    
    ### Add sample ID string Followed # sing in Group Column
    d <- d %>% mutate(sample_id = paste0(Project, " | ", exper_time))
    
    merged_d <- rbind(merged_d, d)
  }
  
  
  return(list(rates = merged_d, Hg_list = Hg_out, PER_background = PER_background_out, OCR_background = OCR_background_out, Zero_measurements = Empty_out, No_cells_measured = No_cells_out, Bad_Cells = Bad_cells_out))
}
# -------------------------------------------------------------- IDENTIFY SINGLE POINT OUTLIARS
# USED IN WORKING PIPELINE

idfy_sinleP_outlier <- function(DT, cut.point, x ) {
  # Identify point outliers based on MAD of the Interval
  # IN:
  # DT        = Seahorse data produced by read_xlsx_set() function
  # cut.point = Treashold for outier removal
  # x         = Varible used, one of: "OCR", "LOCR", "PER", "LPER"
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
    if(length(unique(dm$sample_id)) == 1){ # lm does not work if sample_id only has one level
      fit <- lm(x ~ -1 + Interval + Well, data = dm)
    } else {
      fit <- lm(x ~ -1 + sample_id + Interval + Well, data = dm)
    }
    
    # add column with fitted
    dm$fitted <- fitted(fit)

    # get  of fitted in interval and squared orror of that
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
    cat(i, " Point outliers: ", n_out_p, "--", n_out_p/size*100 ,  "% \n")

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
  cat("Total single point outliers: ", nrow(filter(dm_r, is.out.p == T))/size*100, "% \n" )

  return(dm_r)
}

# ---------------------------------------------------------------------- IDENTIFY OUTLIARS
# WELLS AND SINGLE POINT, COMBINED (Excluded from used Pipeline)
idfy_outlier <- function(DT, x, cut.well, cut.point ){


  dm <- DT
  dm$x <- dm[[x]] # Create column with regressed variable
  i <- 1
  keep <- T
  size <- nrow(dm)
  dm <- dm %>% drop_na()
  while (keep) {
    # --- Fitting linear regresion, separate wells in separate intervals (and sample)
    fit_all <- lm(x ~ -1 + sample_id + Interval + Well, data = dm) # can be one sample_ID
    dm$fitted <- fitted(fit_all)
    dm$residuals = residuals(fit_all)

    # get mean of fitted in interval
    dm <- dm %>%  group_by(sample_id, Interval) %>%
      mutate(int_mean = mean(fitted, na.rm = T),
             n= n(),
             sq_err   = (x-int_mean)^2) %>%
      ungroup()
    # get mean error in well
    dm <- dm %>%  group_by(sample_id, Well) %>%
      mutate(mean_err_sq = mean(sq_err, na.rm = T) ) %>%
      ungroup()
    # Done for Well
    dm <- dm %>% group_by(sample_id) %>%
      mutate(median_mean_sq_err = median(mean_err_sq, na.rm = T ),
             mad_mean_sqE       = mad(mean_err_sq, na.rm = T),   # of wells in particular interval (not all the wells)
             is.out.w           = median_mean_sq_err + cut.well * mad_mean_sqE < mean_err_sq) %>%
      ungroup()
    # Print
    n_out_w <-  nrow(dm %>% filter(is.out.w == TRUE))
    cat(i, " Well outliers: ", n_out_w, "--", n_out_w/size*100 ,  "% \n")

    dm <- dm %>% filter(is.out.w == F)
    if( n_out_w == 0) keep <-FALSE # Stop when no more outliars found
    i <- i+1
  }
  # complete data
  dm_r <- DT %>%  drop_na() %>%
    left_join(dm, by = c("Measurement", "Well", "Group", "Time", "plate_id", "Protocol", "Interval"), suffix  = c("",".y") ) %>%
    mutate(is.out.w = replace_na(is.out.w, T)) %>%
    select(-c(contains("y"), "mad_mean_sqE", "sq_err", "int_mean"))
  # ------------ Single Point Outliars

  i <- 1
  keep <- T
  size <- nrow(dm)
  dm <- dm %>% drop_na()
  while (keep) {

    fit_all <- lm(x ~ -1 + sample_id + Interval + Well, data = dm)
    dm$fitted <- fitted(fit_all)
    dm$residuals = residuals(fit_all)

    dm <- dm %>% group_by(sample_id, Interval) %>%
      mutate(int_mean = mean(fitted, na.rm = T),
             sq_err   = (x-int_mean)^2) %>%
      ungroup()

    # SINGLE
    dm <- dm %>%  group_by(sample_id, Interval) %>%
      mutate(median_sqE = median(sq_err),
             mad_sqE    = mad(sq_err, na.rm=T),
             is.out.p    = median_sqE + cut.point * mad_sqE < sq_err) %>%
      ungroup()

    n_out_p <-  nrow(dm %>% filter(is.out.p == T))
    cat(i, " Point outliers: ", n_out_p, "--", n_out_p/size*100 ,  "% \n")

    dm <- dm %>% filter(is.out.p == F)
    if( n_out_p == 0) keep <-F # Stop when no more outliars found

    i <- i+1

  }

  # complete data
  dm_r <- dm_r %>%
    left_join(dm, by = c("Measurement", "Well", "Group", "Time", "plate_id", "Protocol", "Interval"), suffix = c("",".y") ) %>%
    mutate(is.out.p = ifelse(is.out.w == T, NA, replace_na(is.out.p, T)),
           out      = ifelse(is.out.p == F & is.out.w == F, "NO",ifelse(is.out.w == T, "WELL", "SINGLE"))) %>%
    select(-c(contains("y"), "median_sqE", "mad_sqE", "sq_err", "int_mean"))

  # Print summary
  cat("Total well outliers: ", nrow(filter(dm_r, is.out.w == T))/nrow(dm_r)*100, "% \n" )
  cat("Total single point outliers: ", nrow(filter(dm_r, is.out.p == T))/nrow(dm_r)*100, "% \n" )

  return(dm_r)
}


# -------------------------------------------------------------------------- COMPUTE BIOENERGETICS FOR EACH REPLICATE
compute_bioenergetics_replicate <- function(dm_r, method) {
  dr <- dm_r %>%
    filter(is.out.p == FALSE)

  dr$x <- dr[[method]]
  # we are using median instead !!!
  estim_mean <- dr %>%
    group_by(sample_id, Interval) %>%
    summarise(mean = mean(x), SD = sd(x), SE = sd(x)/sqrt(n()), size = n(), CV = (SD/mean)*100 )

  # form datafames
  estimates  <- as.data.frame(cast(estim_mean, sample_id~Interval, value = "mean"))
  deviations <- as.data.frame(cast(estim_mean, sample_id~Interval, value = "SD"))
  SErrs      <- as.data.frame(cast(estim_mean, sample_id~Interval, value = "SE"))
  numbers    <- as.data.frame(cast(estim_mean, sample_id~Interval, value = "size"))
  CVs        <- as.data.frame(cast(estim_mean, sample_id~Interval, value = "CV"))

  # compute Bioenergetics according to the method used

  if (method == "OCR") {
    # difference based bioenergetics
    bio_e <- estimates %>%
      mutate(Sample           = sample_id,
             Basal.Resp       = Int1 - Int4,
             ATP.linked.Resp  = Int1 - Int2,
             Proton.Leak      = Int2 - Int4,
             Spare.Resp.Cpcty = Int3 - Int1,
             Maximal.Resp     = Int3 - Int4,
             Non.Mito.Resp    = Int4) %>%
      select(-c("Int1", "Int2", "Int3", "Int4", "sample_id"))
    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample           = sd.sample_id,
             Basal.Resp       = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int4^2)/n.Int4)),
             ATP.linked.Resp  = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int2^2)/n.Int2)),
             Proton.Leak      = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int4^2)/n.Int4)),
             Spare.Resp.Cpcty = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int1^2)/n.Int1)),
             Maximal.Resp     = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int4^2)/n.Int4)),
             Non.Mito.Resp    = sd.Int4/sqrt(n.Int4)) %>%
      select(c("Sample", "Basal.Resp", "ATP.linked.Resp", "Proton.Leak", "Spare.Resp.Cpcty", "Maximal.Resp", "Non.Mito.Resp"))


  } else if (method == "LOCR") {
    # Ratio based bioenergetics
    bio_e <- estimates %>%
      mutate(Sample               = sample_id,
             log.Basal.Resp       = Int1 - Int4,
             log.ATP.linked.Resp  = Int1 - Int2,
             log.Proton.Leak      = Int2 - Int4,
             log.Spare.Resp.Cpcty = Int3 - Int1,
             log.Maximal.Resp     = Int3 - Int4,
             log.Non.Mito.Resp    = Int4) %>%
      select(-c("Int1", "Int2", "Int3", "Int4", "sample_id"))
    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample               = sd.sample_id,
             log.Basal.Resp       = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int4^2)/n.Int4)),
             log.ATP.linked.Resp  = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int2^2)/n.Int2)),
             log.Proton.Leak      = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int4^2)/n.Int4)),
             log.Spare.Resp.Cpcty = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int1^2)/n.Int1)),
             log.Maximal.Resp     = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int4^2)/n.Int4)),
             log.Non.Mito.Resp    = sd.Int4/sqrt(n.Int4)) %>%
      select(c("Sample", "log.Basal.Resp", "log.ATP.linked.Resp", "log.Proton.Leak", "log.Spare.Resp.Cpcty",
               "log.Maximal.Resp", "log.Non.Mito.Resp"))

  } else if (method == "PER") {
    bio_e <- estimates %>%
      mutate(Sample            = sample_id,
             Basal.PER         = Int1,
             Max.PER           = Int2,
             Glyco.Rsrv.Cpcty  = Int2 - Int1) %>%
      select(c("Sample", "Basal.PER", "Max.PER", "Glyco.Rsrv.Cpcty"))

    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample            = sd.sample_id,
             Basal.PER         = sqrt(((sd.Int1^2)/n.Int1)),
             Max.PER           = sqrt(((sd.Int2^2)/n.Int2)),
             Glyco.Rsrv.Cpcty  = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int1^2)/n.Int1))) %>%
      select(c("Sample", "Basal.PER", "Max.PER", "Glyco.Rsrv.Cpcty"))
  } else if (method == "LPER") {
    bio_e <- estimates %>%
      mutate(Sample                = sample_id,
             log.Basal.PER         = Int1,
             log.Max.PER           = Int2,
             log.Glyco.Rsrv.Cpcty  = Int2 - Int1) %>%
      select(c("Sample", "log.Basal.PER", "log.Max.PER", "log.Glyco.Rsrv.Cpcty"))

    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample                = sd.sample_id,
             log.Basal.PER         = sqrt(((sd.Int1^2)/n.Int1)),
             log.Max.PER           = sqrt(((sd.Int2^2)/n.Int2)),
             log.Glyco.Rsrv.Cpcty  = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int1^2)/n.Int1))) %>%
      select(c("Sample", "log.Basal.PER", "log.Max.PER", "log.Glyco.Rsrv.Cpcty"))

  }
  return(list(bioenergetics = bio_e, standard.errors = st_errors, estimates = estim_mean ))
}

# -------------------------------------------------------------------------- COMPUTE BIOENERGETICS FOR EACH WELL
compute_bioenergetics_well <- function(dm_r, method) {
  dr <- dm_r %>%
    filter(is.out.p == FALSE)
  
  dr$x <- dr[[method]]
  # we are using median instead !!!
  estim_mean <- dr %>%
    group_by(sample_id, Interval, Well) %>%
    summarise(mean = mean(x), SD = sd(x), SE = sd(x)/sqrt(n()), size = n(), CV = (SD/mean)*100 )
  
  # form datafames
  estimates  <- as.data.frame(cast(estim_mean, sample_id + Well~Interval, value = "mean"))
  deviations <- as.data.frame(cast(estim_mean, sample_id + Well~Interval, value = "SD"))
  SErrs      <- as.data.frame(cast(estim_mean, sample_id + Well~Interval, value = "SE"))
  numbers    <- as.data.frame(cast(estim_mean, sample_id + Well~Interval, value = "size"))
  CVs        <- as.data.frame(cast(estim_mean, sample_id + Well~Interval, value = "CV"))
  
  # compute Bioenergetics according to the method used
  
  if (method == "OCR") {
    # difference based bioenergetics
    bio_e <- estimates %>%
      mutate(Sample           = sample_id,
             Basal.Resp       = Int1 - Int4,
             ATP.linked.Resp  = Int1 - Int2,
             Proton.Leak      = Int2 - Int4,
             Spare.Resp.Cpcty = Int3 - Int1,
             Maximal.Resp     = Int3 - Int4,
             Non.Mito.Resp    = Int4) %>%
      select(-c("Int1", "Int2", "Int3", "Int4", "sample_id"))
    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample           = sd.sample_id,
             Basal.Resp       = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int4^2)/n.Int4)),
             ATP.linked.Resp  = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int2^2)/n.Int2)),
             Proton.Leak      = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int4^2)/n.Int4)),
             Spare.Resp.Cpcty = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int1^2)/n.Int1)),
             Maximal.Resp     = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int4^2)/n.Int4)),
             Non.Mito.Resp    = sd.Int4/sqrt(n.Int4)) %>%
      select(c("Sample","n.Well", "Basal.Resp", "ATP.linked.Resp", "Proton.Leak", "Spare.Resp.Cpcty", "Maximal.Resp", "Non.Mito.Resp"))
    colnames(st_errors)[2] <- "Well"
    
  } else if (method == "LOCR") {
    # Ratio based bioenergetics
    bio_e <- estimates %>% na.omit() %>% 
      mutate(Sample               = sample_id,
             log.Basal.Resp       = Int1 - Int4,
             log.ATP.linked.Resp  = Int1 - Int2,
             log.Proton.Leak      = Int2 - Int4,
             log.Spare.Resp.Cpcty = Int3 - Int1,
             log.Maximal.Resp     = Int3 - Int4,
             log.Non.Mito.Resp    = Int4) %>%
      select(-c("Int1", "Int2", "Int3", "Int4", "sample_id"))
    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample               = sd.sample_id,
             log.Basal.Resp       = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int4^2)/n.Int4)),
             log.ATP.linked.Resp  = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int2^2)/n.Int2)),
             log.Proton.Leak      = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int4^2)/n.Int4)),
             log.Spare.Resp.Cpcty = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int1^2)/n.Int1)),
             log.Maximal.Resp     = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int4^2)/n.Int4)),
             log.Non.Mito.Resp    = sd.Int4/sqrt(n.Int4)) %>%
      select(c("Sample", "n.Well", "log.Basal.Resp", "log.ATP.linked.Resp", "log.Proton.Leak", "log.Spare.Resp.Cpcty",
               "log.Maximal.Resp", "log.Non.Mito.Resp"))
    colnames(st_errors)[2] <- "Well"
    
  } else if (method == "PER") {
    bio_e <- estimates %>% na.omit() %>% 
      mutate(Sample            = sample_id,
             Basal.PER         = Int1,
             Max.PER           = Int2,
             Glyco.Rsrv.Cpcty  = Int2 - Int1) %>%
      select(c("Sample", "Well", "Basal.PER", "Max.PER", "Glyco.Rsrv.Cpcty"))
    
    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample            = sd.sample_id,
             Basal.PER         = sqrt(((sd.Int1^2)/n.Int1)),
             Max.PER           = sqrt(((sd.Int2^2)/n.Int2)),
             Glyco.Rsrv.Cpcty  = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int1^2)/n.Int1))) %>%
      select(c("Sample","n.Well", "Basal.PER", "Max.PER", "Glyco.Rsrv.Cpcty"))
    colnames(st_errors)[2] <- "Well"
  } else if (method == "LPER") {
    bio_e <- estimates %>% na.omit() %>% 
      mutate(Sample                = sample_id,
             log.Basal.PER         = Int1,
             log.Max.PER           = Int2,
             log.Glyco.Rsrv.Cpcty  = Int2 - Int1) %>%
      select(c("Sample", "Well", "log.Basal.PER", "log.Max.PER", "log.Glyco.Rsrv.Cpcty"))
    
    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample                = sd.sample_id,
             log.Basal.PER         = sqrt(((sd.Int1^2)/n.Int1)),
             log.Max.PER           = sqrt(((sd.Int2^2)/n.Int2)),
             log.Glyco.Rsrv.Cpcty  = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int1^2)/n.Int1))) %>%
      select(c("Sample", "n.Well", "log.Basal.PER", "log.Max.PER", "log.Glyco.Rsrv.Cpcty"))
    colnames(st_errors)[2] <- "Well"
    
  }
  return(list(bioenergetics = bio_e, standard.errors = st_errors, estimates = estim_mean))
}
# -------------------------------------------------------------------------- COMPUTE BIOENERGETICS FOR EACH SAMPLE
compute_bioenergetics_sample <- function(dm_r, method) {
  dr <- dm_r %>%
    filter(is.out.p == FALSE)
  
  dr$x <- dr[[method]]
  # we are using median instead !!!
  estim_mean <- dr %>%
    group_by(Group, Interval) %>%
    summarise(mean = mean(x), SD = sd(x), SE = sd(x)/sqrt(n()), size = n(), CV = (SD/mean)*100 )
  
  # form datafames
  estimates  <- as.data.frame(cast(estim_mean, Group~Interval, value = "mean"))
  deviations <- as.data.frame(cast(estim_mean, Group~Interval, value = "SD"))
  SErrs      <- as.data.frame(cast(estim_mean, Group~Interval, value = "SE"))
  numbers    <- as.data.frame(cast(estim_mean, Group~Interval, value = "size"))
  CVs        <- as.data.frame(cast(estim_mean, Group~Interval, value = "CV"))
  
  # compute Bioenergetics according to the method used
  
  if (method == "OCR") {
    # difference based bioenergetics
    bio_e <- estimates %>%
      mutate(Sample           = Group,
             Basal.Resp       = Int1 - Int4,
             ATP.linked.Resp  = Int1 - Int2,
             Proton.Leak      = Int2 - Int4,
             Spare.Resp.Cpcty = Int3 - Int1,
             Maximal.Resp     = Int3 - Int4,
             Non.Mito.Resp    = Int4) %>%
      select(-c("Int1", "Int2", "Int3", "Int4", "Group"))
    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample           = sd.Group,
             Basal.Resp       = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int4^2)/n.Int4)),
             ATP.linked.Resp  = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int2^2)/n.Int2)),
             Proton.Leak      = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int4^2)/n.Int4)),
             Spare.Resp.Cpcty = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int1^2)/n.Int1)),
             Maximal.Resp     = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int4^2)/n.Int4)),
             Non.Mito.Resp    = sd.Int4/sqrt(n.Int4)) %>%
      select(c("Sample", "Basal.Resp", "ATP.linked.Resp", "Proton.Leak", "Spare.Resp.Cpcty", "Maximal.Resp", "Non.Mito.Resp"))
    
    
  } else if (method == "LOCR") {
    # Ratio based bioenergetics
    bio_e <- estimates %>%
      mutate(Sample               = Group,
             log.Basal.Resp       = Int1 - Int4,
             log.ATP.linked.Resp  = Int1 - Int2,
             log.Proton.Leak      = Int2 - Int4,
             log.Spare.Resp.Cpcty = Int3 - Int1,
             log.Maximal.Resp     = Int3 - Int4,
             log.Non.Mito.Resp    = Int4) %>%
      select(-c("Int1", "Int2", "Int3", "Int4", "Group"))
    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample               = sd.Group,
             log.Basal.Resp       = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int4^2)/n.Int4)),
             log.ATP.linked.Resp  = sqrt(((sd.Int1^2)/n.Int1)+((sd.Int2^2)/n.Int2)),
             log.Proton.Leak      = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int4^2)/n.Int4)),
             log.Spare.Resp.Cpcty = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int1^2)/n.Int1)),
             log.Maximal.Resp     = sqrt(((sd.Int3^2)/n.Int3)+((sd.Int4^2)/n.Int4)),
             log.Non.Mito.Resp    = sd.Int4/sqrt(n.Int4)) %>%
      select(c("Sample", "log.Basal.Resp", "log.ATP.linked.Resp", "log.Proton.Leak", "log.Spare.Resp.Cpcty",
               "log.Maximal.Resp", "log.Non.Mito.Resp"))
    
  } else if (method == "PER") {
    bio_e <- estimates %>%
      mutate(Sample            = Group,
             Basal.PER         = Int1,
             Max.PER           = Int2,
             Glyco.Rsrv.Cpcty  = Int2 - Int1) %>%
      select(c("Sample", "Basal.PER", "Max.PER", "Glyco.Rsrv.Cpcty"))
    
    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample            = sd.Group,
             Basal.PER         = sqrt(((sd.Int1^2)/n.Int1)),
             Max.PER           = sqrt(((sd.Int2^2)/n.Int2)),
             Glyco.Rsrv.Cpcty  = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int1^2)/n.Int1))) %>%
      select(c("Sample", "Basal.PER", "Max.PER", "Glyco.Rsrv.Cpcty"))
  } else if (method == "LPER") {
    bio_e <- estimates %>%
      mutate(Sample                = Group,
             log.Basal.PER         = Int1,
             log.Max.PER           = Int2,
             log.Glyco.Rsrv.Cpcty  = Int2 - Int1) %>%
      select(c("Sample", "log.Basal.PER", "log.Max.PER", "log.Glyco.Rsrv.Cpcty"))
    
    # standard errors of mean differences
    sd_n   <- cbind(sd = deviations, n = numbers)
    st_errors <- sd_n %>%
      mutate(Sample                = sd.Group,
             log.Basal.PER         = sqrt(((sd.Int1^2)/n.Int1)),
             log.Max.PER           = sqrt(((sd.Int2^2)/n.Int2)),
             log.Glyco.Rsrv.Cpcty  = sqrt(((sd.Int2^2)/n.Int2)+((sd.Int1^2)/n.Int1))) %>%
      select(c("Sample", "log.Basal.PER", "log.Max.PER", "log.Glyco.Rsrv.Cpcty"))
    
  }
  return(list(bioenergetics = bio_e, standard.errors = st_errors, estimates = estim_mean ))
}
