# ------------------------------INSTALL NECESSARY PACKAGES

packages <- c('dplyr', 'ggplot2', 'tidyr', 'tidyverse', 'readxl', 'here', 'reshape2', "reshape", 'magrittr', 'kableExtra', "shiny", "ggbeeswarm")
missing <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(missing)) install.packages(missing)