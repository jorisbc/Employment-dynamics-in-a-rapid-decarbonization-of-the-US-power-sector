rm(list=ls())
library(reshape2)
library(readxl)

tech.crosswalk <- as.data.frame( read_xlsx('data_other/match_tech_datasets.xlsx') )

case.ref <- read.csv('D:/Science/data_other/energy/NREL Cambium 2021 scenarios/StdScen21_MidCase_annual_national/StdScen21_MidCase_annual_national.csv', skip = 4)
case.35 <- read.csv('D:/Science/data_other/energy/NREL Cambium 2021 scenarios/StdScen21_MidCase95by2035_annual_national/StdScen21_MidCase95by2035_annual_national.csv', skip = 4)
case.50 <- read.csv('D:/Science/data_other/energy/NREL Cambium 2021 scenarios/StdScen21_MidCase95by2050_annual_national/StdScen21_MidCase95by2050_annual_national.csv', skip = 4)

gen_name <- sort( colnames(case.ref)[77:94] )
# we don't include DAC [114] and Canada [95] imports
gen_name <- gen_name[-4]
# we exclude coal CCS (always zero)


# Generation data
genreflong <- cbind.data.frame( melt( case.ref[,c('t',gen_name)], id.vars = 't', value.name = 'generation', variable.name = 'technology' ), scenario = 'Reference')
gen35long <- cbind.data.frame( melt( case.35[,c('t',gen_name)], id.vars = 't', value.name = 'generation', variable.name = 'technology' ), scenario  = '95% by 2035')
gen50long <- cbind.data.frame( melt( case.50[,c('t',gen_name)], id.vars = 't', value.name = 'generation', variable.name = 'technology' ), scenario  = '95% by 2050')

genlong <- rbind(genreflong, gen35long, gen50long)
colnames(genlong)[1] <- 'years'
genlong$technology <- tech.crosswalk$scenarios[ match(genlong$technology, tech.crosswalk$gen_name) ]

# combining BECCS and Biomass
genlong <- aggregate( data.frame(generation = genlong$generation), by = list(years = genlong$years, technology = genlong$technology, scenario = genlong$scenario), FUN=sum )


# Capacity data
cap_name <- sort( colnames(case.ref)[96:113] )
cap_name <- cap_name[-4]
# we exclude coal CCS (always zero)

# CAPACITY data
capreflong <- cbind.data.frame( melt( case.ref[,c('t',cap_name)], id.vars = 't', value.name = 'capacity', variable.name = 'technology' ), scenario = 'Reference')
cap35long <- cbind.data.frame( melt( case.35[,c('t',cap_name)], id.vars = 't', value.name = 'capacity', variable.name = 'technology' ), scenario  = '95% by 2035')
cap50long <- cbind.data.frame( melt( case.50[,c('t',cap_name)], id.vars = 't', value.name = 'capacity', variable.name = 'technology' ), scenario  = '95% by 2050')

caplong <- rbind(capreflong, cap35long, cap50long)
colnames(caplong)[1] <- 'years'
caplong$technology <- tech.crosswalk$scenarios[ match(caplong$technology, tech.crosswalk$cap_name) ]

# combining BECCS and Biomass
caplong <- aggregate( data.frame(capacity = caplong$capacity), by = list(years = caplong$years, technology = caplong$technology, scenario = caplong$scenario), FUN=sum )


genlong$generation <- genlong$generation/10^3 # transform into GWh


write.csv(caplong, file = 'data_out/scenarios_cambium_capacity.csv')
write.csv(genlong, file = 'data_out/scenarios_cambium_generation.csv')