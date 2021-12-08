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
write_csv(be, "PostAnalysisSeahorse/BioeNormOCR.csv")


### READ The data
pd <- read.csv("PostAnalysisSeahorse/All PBMC data Age_Gender_111119_2.csv")
be <- read.csv("PostAnalysisSeahorse/BioeNormOCR - BioeNormOCR.csv")

pd <- pd[,1:7]


### COMBINE INFORMATION FROM STUDIES
# put together log and normal from ecar ocar separately

beOCR <- merge(log.be$bioenergetics, norm.be$bioenergetics)

beECAR <-merge(log.be.ECAR$bioenergetics, norm.be.ECAR$bioenergetics)

be <- merge(beOCR, beECAR , by = "Sample", all = T)

pd$ID <- paste0(pd$Code, "-", pd$Project, "-", pd$Blood.extraction.date)

be$Code <- sub("\\|.*", "", be$Sample)
be$time <- sub(".*\\|", "", be$Sample)
be$Code <- gsub(" ", "", be$Code)
be$Project <- sub(".*-", "", be$Code)
be$ID   <- paste0(be$Code, "-", be$time)



overlap <- intersect(be$ID, pd$ID)
# glue data
p <- filter(pd, ID %in% overlap)
b <- filter(be, ID %in% overlap)

combined <- merge(p, b, by = "ID")

pNotOverlap <- filter(pd, ! ID  %in% overlap)
bNotOverlap <- filter(be, ! ID  %in% overlap)

write_csv(bNotOverlap, "FilesFromPipeline.csv")
write_csv(pNotOverlap, "FilesFromYourDescription.csv")


# glue data ECAR

write_csv(combined, "PostAnalysisSeahorse/combinedOCR_ECAR_LOG.csv")

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

write_csv(d, "complete_data_seahorse.csv")

d <- dm %>%  filter(! sample_id %in% unique(three_fixed$sample_id))

hf32 <- read_xlsx_set(here("Data/INPUT/Data for pipeline/a/"),".xlsx")
hf32 <- hf32$rates

dm<-filter(dm, sample_id != "32Ab- | 15-05-2018")
dm <- rbind(three_fixed, d)
dm <- dm %>% filter(Measurement < 10 )
write_csv(sumar, "filteredHgmm.csv")

###

two_b <- read_xlsx_set("Data/INPUT/Data for pipeline/temp/", ".xlsx")
two_b <- two_b$rates
dm <- dm %>% filter(! sample_id %in% c("3_1x-DB2 | 20-02-2018", "3_1y-DB2 | 20-02-2018"))
d <- rbind(dm, two_b)

dm<-d



