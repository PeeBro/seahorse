# ------------------------------INSTALL NECESSARY PACKAGES

packages <- c('dplyr', 'ggplot2', 'tidyr', 'tidyverse', 'readxl', 'here', 'reshape2', 'magrittr', 'kableExtra', "shiny")
missing <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(missing)) install.packages(missing)