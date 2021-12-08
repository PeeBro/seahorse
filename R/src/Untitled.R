# ======================================= FUNCTIONS FOR SEAHORSE PIPELINE
# Author: Matej Stevuliak

dm <- dm$rates
sam <- dm %>%
  group_by(sample_id) %>%
  summarise(N = n(), Project = unique(Project), max(Measurement))

sam <- three_inj %>%
  group_by(sample_id) %>%
  summarise(N = n(), Project = unique(Project), max(Measurement))


samf <- d_OCR %>%
  group_by(sample_id) %>%
  summarise

be <- norm.be$bioenergetics

write_csv(sam, "final-sample_ids.csv")
write_csv(be, "../PostAnalysisSeahorse/BioeNormOCR.csv")


### READ The data
pd <- read.csv("PostAnalysisSeahorse/Age and gender_pipeline_R - Age-gender ALL.csv")
dm <- read.csv("PostAnalysisSeahorse/complete_dataset.csv")
be <- read.csv("PostAnalysisSeahorse/BioeNormOCR.csv")

### COMBINE INFORMATION FROM STUDIES

pd$ID <- paste0(pd$Code, "-", pd$Project, "-", pd$Blood.extraction.date)
be$Code <- sub("\\|.*", "", be$Sample)
be$Code <- gsub(" ", "", be$Code)
be$Project <- sub(".*-", "", be$Code)
be$ID   <- paste0(be$Code, "-", be$time)

pbmc$Code <- paste0(pbmc$Code, "-", pbmc$Project, "-", pbmc$Blood.extraction.date)
pbmc$ID <- pbmc$Code

overlap <- intersect(be$ID, pd$ID)
# glue data
p <- filter(pd, ID %in% overlap)
b <- filter(be, ID %in% overlap)

combined <- merge(p, b, by = "ID")

pNotOverlap <- filter(pd, ! ID  %in% overlap)
bNotOverlap <- filter(be, ! ID  %in% overlap)

write_csv(bNotOverlap, "FilesFromPipeline.csv")
write_csv(pNotOverlap, "FilesFromYourDescription.csv")
write_csv(combined, "SuccessfullyCombined.csv")

# glue data ECAR
be.ECAR <- log.be$bioenergetics
# finding overlap
overlap <- intersect(be.ECAR$Sample, combinedOCR_ECAR$Sample)

p <- filter(be.ECAR, Sample %in% overlap)
b <- filter(combinedOCR_ECAR, Sample %in% overlap)
combinedOCR_ECAR <- merge(p, b, by = "Sample")
write_csv(combinedOCR_ECAR, "PostAnalysisSeahorse/combinedOCR_ECAR_LOG.csv")

# glue data LOG

###### ESTIMATES

estm <- norm.be$estimates
estimates  <- cast(estm, sample_id~Interval, value = "mean")
write_csv(estimates, "estimates.csv")

#### 3 injections
dm <- dm %>% filter(Project != "TEST" & ! sample_id %in% list_ID  )

three <- dm %>% filter(Measurement > 9 )
list_ID   <- unique(three$sample_id)
three_inj <- dm %>% filter(sample_id %in% list_ID )
glyco <- three_inj %>% filter(Protocol == "Glyc" & Measurement > 3)
glyco <- glyco %>% mutate(Measurement = Measurement -3)

three_inj <- three_inj %>% filter(Measurement >3)

glyco$Interval <- mapply(get_intervals, glyco$Measurement, glyco$Protocol)

# fccp OLIGO

fccp_oligo <- three_inj %>% filter(Protocol != "Glyc" & Measurement < 9)
three_fixed <- rbind(fccp_oligo, glyco)

write_csv(three_fixed, "3_injections_fixed.csv")


dm <- rbind(three_fixed, dm)

write_csv(dm, "complete_data_seahorse.csv")

