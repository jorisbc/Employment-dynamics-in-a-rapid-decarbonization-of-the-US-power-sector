rm(list=ls())
library(readxl)
library(ggplot2); theme_set(theme_bw())
library(reshape2)

### GETTING HISTORICAL DATA

## Installed capacity in MW
convcap.unmod <- as.data.frame(read_xlsx('../../data_other/energy/epa_2020/epa_04_08_a.xlsx', range = "A6:Q16", col_names = F ) )
convcap.raw <- convcap.unmod[ c(1, seq(2,17,2)) ]
convcap.cf <- convcap.unmod[ c(1, seq(3,17,2)) ]
colnames(convcap.raw) <- colnames(convcap.cf) <- c("years", "coal", "gas_cc", "gas_ct", "gas_st", "gas_ic", "oil_st", "oil_gas", "oil_ic")

## note:
# oil_gas and oil_ic have almost zero CF -> we don't include those
# gas_st + gas_ic + oil_st = 89GW just as indicated in Solar Futures Study for O-G-S (the three categories have very similar CF too)
# we don't have O-G-S price data, we treat them as GAS CT

## aggregate to basic tech categories
convcap <- convcap.raw[,c(1,2,3)]
convcap$gas_ct <- rowSums( convcap.raw[,c('gas_ct', 'gas_ic', 'gas_st', 'oil_st')] )  # we make aggregated GAS CT category

renewcap.unmod <- as.data.frame(read_xlsx('../../data_other/energy/epa_2020/epa_04_08_b.xlsx', range = "A6:S16", col_names = F ) )
renewcap.raw <- renewcap.unmod[ c(1, seq(2,19,2)) ]
renewcap.cf <- renewcap.unmod[ c(1, seq(3,19,2)) ]
colnames(renewcap.raw) <- colnames(renewcap.cf) <- c("years", "geo", "hydro", "nuclear", "bio_other", 'gas_other', "pv_util", "csp", "wind_on", "wood")

renewcap <- renewcap.raw[,c(1,2,3,4,7,8,9)]
renewcap$bio <- rowSums(renewcap.raw[,c('bio_other','gas_other')])  # bio is the sum of Other biomass and Other gas (to be consistent with the scenarios)

# small scale PV (we don't have time-averaged data here and no capacity factors)
renewcap$pv_small <- unlist(as.data.frame(read_xlsx('../../data_other/energy/epa_2020/epa_04_02_b.xlsx', range = "I7:I17", col_names = "pv_dist", col_types = 'numeric' ) ) )

# phs + battery
phsbatt.unmod <- as.data.frame(read_xlsx('../../data_other/energy/epa_2020/epa_04_08_c.xlsx', range = "B5:E12", col_names = F, col_types = 'numeric' ) )
renewcap$battery <- c(NA,NA,NA,phsbatt.unmod[,1])  # this battery data is fairly consistent with the july_generator2022.xlsx file (note the massive ramp up in 2021 and 2022)
renewcap$phs <- c(NA,NA,NA,phsbatt.unmod[,3])

# combine fossil and renewable data
dcap.w <- data.frame( years = convcap$years,
                    Coal = convcap$coal,
                    `Gas CC CCS` = 0,
                    `Gas CC` = convcap$gas_cc, 
                    `Gas CT` = convcap$gas_ct,
                    Nuclear = renewcap$nuclear,
                    Hydro = renewcap$hydro,
                    Biomass =renewcap$bio, 
                    Geothermal = renewcap$geo,
                    `Wind on` = renewcap$wind,
                    `Wind off` = 0,
                    CSP = renewcap$csp,
                    `PV utility` = renewcap$pv_util,
                    `PV dist` = renewcap$pv_small,
                    Batteries = renewcap$battery,
                    PHS = renewcap$phs, 
                    check.names = F # needed to make colnames correct
)
# dcap.w
# read.csv('data_out/historical_capacity_by_source.csv')[,-1] # check with old version



## historical CF data
ct.tech <- c('gas_ct', 'gas_ic', 'gas_st', 'oil_st')
gas_ct.CF <- rowSums( sweep(convcap.raw[,ct.tech], 1, rowSums(convcap.raw[,ct.tech]), '/') * convcap.cf[,ct.tech] ) # weighted average

dCF.w <- data.frame( years = convcap$years,
                      Coal = convcap.cf$coal,
                      `Gas CC CCS` = 0,
                      `Gas CC` = convcap.cf$gas_cc, 
                      `Gas CT` = gas_ct.CF,
                      Nuclear = renewcap.cf$nuclear,
                      Hydro = renewcap.cf$hydro,
                      Biomass =renewcap.cf$bio, 
                      Geothermal = renewcap.cf$geo,
                      `Wind on` = renewcap.cf$wind,
                      `Wind off` = 0,
                      CSP = renewcap.cf$csp,
                      `PV utility` = renewcap.cf$pv_util,
                      `PV dist` = NA,
                      Batteries =  c(NA,NA,NA,phsbatt.unmod[,2]),
                      PHS =  c(NA,NA,NA,phsbatt.unmod[,4]), 
                      check.names = F # needed to make colnames correct
)

# USE CF AND CAPACITY DATA TO MAKE NET GENERATION DATA
# round( dcap.w * dCF.w *8760/1000 )
# read.csv('data_out/historical_generation_by_source.csv')[,-1] # check with old version: very similar values (old version has errors in gas)

## GENERATION IN GWH
dgen.w <- dcap.w * dCF.w *8760/1000 
dgen.w$`PV dist` <- unlist(as.data.frame(read_xlsx('../../data_other/energy/epa_2020/epa_03_01_b.xlsx', range = "L6:L16", col_names = "pv_dist", col_types = 'numeric' ) ) )
dgen.w$years <- dcap.w$years

dcap <- melt(dcap.w, id.vars = 'years', variable.name = 'technology', value.name = 'capacity')
dgen <- melt(dgen.w, id.vars = 'years', variable.name = 'technology', value.name = 'generation')

write.csv(dgen.w, file = 'data_out/historical_generation_by_source_wide.csv')
write.csv(dgen, file = 'data_out/historical_generation_by_source_long.csv')

write.csv(dcap.w, file = 'data_out/historical_capacity_by_source_wide.csv')
write.csv(dcap, file = 'data_out/historical_capacity_by_source_long.csv')

# 
# 
# ### MERGE WITH NREL SCENARIOS TO VISUALIZE
# 
# ## CAPACITY
# 
# nrelcap <- read.csv('data_out/scenarios_cambium_capacity.csv')[,-1]
# capall <-  
#   rbind(
#     cbind.data.frame(dcap, scenario = '95% by 2035'),
#     cbind.data.frame(dcap, scenario = 'Reference'),
#     cbind.data.frame(dcap, scenario = '95% by 2050') ,
#     nrelcap[,c(1,2,4,3)] )
# 
# 
# ## GENERATION
# nrelgen <- read.csv('data_out/scenarios_cambium_generation.csv')[,-1]
# 
# genall <-  
#   rbind(
#     cbind.data.frame(dgen, scenario = '95% by 2035'),
#     cbind.data.frame(dgen, scenario = 'Reference'),
#     cbind.data.frame(dgen, scenario = '95% by 2050') ,
#     nrelgen[,c(1,2,4,3)] )
# 
# 
# 
# source('source_code/colorcodes.R')
# ## visualizations
# ggplot(capall, aes(x=years, color = technology, fill = technology, y=capacity/10^3)) +
#   theme_light() +
#   geom_area() +
#   scale_fill_manual(values = my_colors) +
#   scale_color_manual(values = my_colors) +
#   facet_wrap('scenario')  + 
#   #  scale_x_continuous(breaks = seq(2020,2050,10)) +
#   theme(
#     axis.title = element_text(size = 22),
#     axis.text = element_text(size = 22),
#     axis.ticks.length.x = unit(.3,'cm'),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.text = element_text(size = 22),
#     legend.title = element_blank(),
#     legend.background = element_blank(),
#     legend.key.size = unit(1,"cm"),
#     legend.position = 'bottom'
#   ) + xlab('') + ylab('Capacity GW')
# 
# 
# ggplot(genall, aes(x=years, color = technology, fill = technology, y=generation/10^3)) +
#   theme_light() +
#   geom_area() +
#   scale_fill_manual(values = my_colors) +
#   scale_color_manual(values = my_colors) +
#   facet_wrap('scenario')  + 
#   #  scale_x_continuous(breaks = seq(2020,2050,10)) +
#   theme(
#     axis.title = element_text(size = 22),
#     axis.text = element_text(size = 22),
#     axis.ticks.length.x = unit(.3,'cm'),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.text = element_text(size = 22),
#     legend.title = element_blank(),
#     legend.background = element_blank(),
#     legend.key.size = unit(1,"cm"),
#     legend.position = 'bottom'
#   ) + xlab('') + ylab('Generation TWh')