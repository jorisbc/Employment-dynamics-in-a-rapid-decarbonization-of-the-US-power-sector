rm(list = ls())
library(readxl)
library(xtable)

tech.crosswalk <- as.data.frame( read_xlsx('data_other/match_tech_datasets.xlsx', sheet = 'table_match_latex') )

print( xtable(tech.crosswalk[order(tech.crosswalk$`Cambium technologies`),]), include.rownames = F)