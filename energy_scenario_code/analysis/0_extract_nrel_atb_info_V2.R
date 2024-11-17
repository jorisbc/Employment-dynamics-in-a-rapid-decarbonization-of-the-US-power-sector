rm(list=ls())

library(readxl)
library(reshape2)
library(ggplot2) ; theme_set(theme_bw())

years <- 2019:2050


#### FUEL COSTS FOR GAS AND COAL ----------------------------------------------


# no fuel costs for gas and coal given in ATB 2021

# GAS FUEL given in $/MMBTU
## Heat rate gas [MMBtu/Mwh]: ATB 2020 
mmbtu.factor.Gas.CT <- 9.51
mmbtu.factor.Gas.CC	<- 6.40
mmbtu.factor.Gas.CC.CCS <-	7.53

# use 2020 data for fuels on GAS and COAL (not in 2021 version yet)
dd <- read.csv("D:/Science/data_other/energy/NREL_ATB/ATBe_2020.csv")

# coal fuel costs
d.coal.20 <- dd[dd$technology == 'Coal' & 
                  dd$techdetail == "newAvgCF" &                 # no CSS or IGCC assumption
                  dd$core_metric_variable %in% years &
                  dd$scenario == 'Moderate' &
                  dd$crpyears == 20 &                           # no difference across capex, opex and fuel
                  dd$core_metric_case == 'R&D',                 # cost changes only come from R&D (not from policy)
]

coal.fuel <- d.coal.20[ d.coal.20$core_metric_parameter == 'Fuel', 'value'] # COAL fuel correctly given $/MWh
names(coal.fuel) <- d.coal.20[ d.coal.20$core_metric_parameter == 'Fuel', 'core_metric_variable']


# gas fuel costs (we distinguish different gas types)
d.gas.20 <- dd[dd$technology == 'NaturalGas' & 
                 dd$core_metric_variable %in% years &
                 dd$scenario == 'Moderate' &
                 dd$crpyears == 20 &                           # no difference across capex, opex and fuel
                 dd$core_metric_case == 'R&D',                 # cost changes only come from R&D (not from policy)
]


# combined cycle
gasCC.fuel <- d.gas.20[d.gas.20$techdetail == "CCAvgCF" & d.gas.20$core_metric_parameter == 'Fuel', 'value'] * mmbtu.factor.Gas.CC
names(gasCC.fuel) <- d.gas.20[d.gas.20$techdetail == "CCAvgCF" & d.gas.20$core_metric_parameter == 'Fuel', 'core_metric_variable']

# combined cycle plus carbon capture
gasCCCCS.fuel <- d.gas.20[d.gas.20$techdetail == "CCCCSAvgCF" & d.gas.20$core_metric_parameter == 'Fuel', 'value'] * mmbtu.factor.Gas.CC.CCS
names(gasCCCCS.fuel) <- d.gas.20[d.gas.20$techdetail == "CCCCSAvgCF" & d.gas.20$core_metric_parameter == 'Fuel', 'core_metric_variable']

# combustion turbine
gasCT.fuel <- d.gas.20[d.gas.20$techdetail == "CTAvgCF" & d.gas.20$core_metric_parameter == 'Fuel', 'value'] * mmbtu.factor.Gas.CT
names(gasCT.fuel) <- d.gas.20[d.gas.20$techdetail == "CTAvgCF" & d.gas.20$core_metric_parameter == 'Fuel', 'core_metric_variable']


cbind.data.frame(coal.fuel, gasCC.fuel, gasCT.fuel, gasCCCCS.fuel)




#### COST DATA ATB 2021 ------------------------------------------------------


# use 2021 data for all other technologies
cost.factors <- c("CAPEX", "Fixed O&M", "Variable O&M", "Fuel")
d <- read.csv("D:/Science/data_other/energy/NREL_ATB/ATBe_2021.csv")
d <- d[,-c(1,2,3,14)]
d <- d[d$technology != "AEO",]
d <- d[d$scenario != "*",]
d <- d[d$techdetail != '*',]
d <- d[d$core_metric_case == 'R&D',]
d <- d[d$core_metric_parameter %in% cost.factors, ]
d <- d[d$core_metric_variable %in% years, ]

d <- d[,c('technology', 'techdetail', 'core_metric_variable', 'scenario', 'crpyears', 'core_metric_case', 'core_metric_parameter', 'units', "value")]
colnames(d)[3] <- 'years'
rownames(d) <- NULL

technologies_1 <- c( "Biopower", "CSP","Coal_FE",                      
                     "CommPV", "Geothermal", "Hydropower", "LandbasedWind",                
                     "NaturalGas_FE", "OffShoreWind", "Pumped Storage Hydropower",   
                     "ResPV", "Utility-Scale Battery Storage", "UtilityPV",                    
                     "Nuclear")


battery.l <- d[d$technology == 'Utility-Scale Battery Storage' & d$techdetail == "4Hr Battery Storage" & d$scenario == 'Moderate' & d$crpyears == 30,]
battery <- dcast(dat = battery.l, formula = years ~ core_metric_parameter, value.var = 'value')

bio.l <- d[d$technology == 'Biopower' & d$scenario == 'Moderate' & d$crpyears == 30,]
bio <- dcast(dat = bio.l, formula = years ~ core_metric_parameter, value.var = 'value')

coal.l <- d[d$technology == 'Coal_FE' & d$techdetail == "newAvgCF" & d$scenario == 'Moderate' & d$crpyears == 30,]
coal <- dcast(dat = coal.l, formula = years ~ core_metric_parameter, value.var = 'value')

csp.l <- d[d$technology == 'CSP' & d$techdetail == "Class3" & d$scenario == 'Moderate' & d$crpyears == 30,]
csp <- dcast(dat = csp.l, formula = years ~ core_metric_parameter, value.var = 'value')

pvcomm.l <- d[d$technology == 'CommPV' & d$techdetail == "Class5" & d$scenario == 'Moderate' & d$crpyears == 30,]
pvcomm <- dcast(dat = pvcomm.l, formula = years ~ core_metric_parameter, value.var = 'value')

pvres.l <- d[d$technology == 'ResPV' & d$techdetail == "Class5" & d$scenario == 'Moderate' & d$crpyears == 30,]
pvres <- dcast(dat = pvres.l, formula = years ~ core_metric_parameter, value.var = 'value')


# assume given empirical split (EIA EPA Table 6.1.B.):
# Residential PV: 21022.1 MW (2021)
# Total small scale (residential + commercial + industrial): 32972.3 MW (2021)
res.share <- 21022.1/32972.3
pvdist <- cbind.data.frame(years = pvcomm$years, pvcomm[,-1] * (1-res.share) + pvres[,-1] * res.share )

pvutil.l <- d[d$technology == 'UtilityPV' & d$techdetail == "Class5" & d$scenario == 'Moderate' & d$crpyears == 30,]
pvutil <- dcast(dat = pvutil.l, formula = years ~ core_metric_parameter, value.var = 'value')

gasct.l <- d[d$technology == 'NaturalGas_FE' &  d$techdetail == 'CTAvgCF' & d$scenario == 'Moderate' & d$crpyears == 30,]
gasct <- dcast(dat = gasct.l, formula = years ~ core_metric_parameter, value.var = 'value')

gascc.l <- d[d$technology == 'NaturalGas_FE' &  d$techdetail == 'CCAvgCF' & d$scenario == 'Moderate' & d$crpyears == 30,]
gascc <- dcast(dat = gascc.l, formula = years ~ core_metric_parameter, value.var = 'value')

gasccccs.l <- d[d$technology == 'NaturalGas_FE' &  d$techdetail == 'CCCCSAvgCF' & d$scenario == 'Moderate' & d$crpyears == 30,]
gasccccs <- dcast(dat = gasccccs.l, formula = years ~ core_metric_parameter, value.var = 'value')

hydro.l <- d[d$technology == 'Hydropower' &  d$techdetail == 'NPD1' & d$scenario == 'Moderate' & d$crpyears == 30,]
hydro <- dcast(dat = hydro.l, formula = years ~ core_metric_parameter, value.var = 'value')

phs.l <- d[d$technology == "Pumped Storage Hydropower" &  d$techdetail == 'Class 3' & d$scenario == 'Moderate' & d$crpyears == 30,]
phs <- dcast(dat = phs.l, formula = years ~ core_metric_parameter, value.var = 'value')

windon.l <- d[d$technology == "LandbasedWind" & d$techdetail == 'Class4' & d$scenario == 'Moderate' & d$crpyears == 30,]
windon <- dcast(dat = windon.l, formula = years ~ core_metric_parameter, value.var = 'value')

windoff.l <- d[d$technology == "OffShoreWind" & d$techdetail == 'Class3' & d$scenario == 'Moderate' & d$crpyears == 30,]
windoff <- dcast(dat = windoff.l, formula = years ~ core_metric_parameter, value.var = 'value')

geo.l <- d[d$technology == "Geothermal" & d$techdetail == 'HydroFlash' & d$scenario == 'Moderate' & d$crpyears == 30,]
geo <- dcast(dat = geo.l, formula = years ~ core_metric_parameter, value.var = 'value')

nuclear.l <- d[d$technology == "Nuclear" & d$techdetail == 'Nuclear' & d$scenario == 'Moderate' & d$crpyears == 30,]
nuclear <- dcast(dat = nuclear.l, formula = years ~ core_metric_parameter, value.var = 'value')

gascc$Fuel <- gasCC.fuel
gasct$Fuel <- gasCT.fuel
gasccccs$Fuel <- gasCCCCS.fuel
coal$Fuel <- coal.fuel


#### COMBINE DATA ------------------------------------------------------

techcost.l <- list(
  Batteries = battery,
  BECCS = bio,
  Biomass = bio,
  Coal = coal,
  CSP = csp,
  `PV dist` = pvdist,
  `Gas CC CCS` = gasccccs,
  `Gas CC` = gascc,
  `Gas CT` = gasct,
  Geothermal = geo,
  Hydro = hydro,
  Nuclear = nuclear,
  PHS = phs,
  `PV utility` = pvutil,
  `Wind off` = windoff,
  `Wind on` = windon )

battery$Fuel <- 0
bio <- bio[, c('years','CAPEX', 'Fixed O&M', 'Variable O&M','Fuel')]
csp$Fuel <- 0
pvdist$Fuel <- 0
geo$Fuel <- 0
hydro$Fuel <- 0
phs$Fuel <- 0
pvutil$Fuel <- 0
windoff$Fuel <- 0
windon$Fuel <- 0


techcost.d <- rbind.data.frame(
  cbind.data.frame(technology = 'Coal', coal),
  cbind.data.frame(technology = 'Gas CC', gascc),
  cbind.data.frame(technology = 'Gas CC CCS', gasccccs),
  cbind.data.frame(technology = 'Gas CT', gasct),
  cbind.data.frame(technology = 'Nuclear', nuclear),
  cbind.data.frame(technology = 'Hydro', hydro),
  cbind.data.frame(technology = 'PHS', hydro),   # we assume hydro for PHS since PHS data is incomplete
  cbind.data.frame(technology = 'Biomass', bio),
  cbind.data.frame(technology = 'Geothermal', geo),
  cbind.data.frame(technology = 'CSP', csp),
  cbind.data.frame(technology = 'PV utility', pvutil),
  cbind.data.frame(technology = 'PV dist', pvdist),
  cbind.data.frame(technology = 'Wind on', windon),
  cbind.data.frame(technology = 'Wind off', windoff),
  cbind.data.frame(technology = 'Batteries', battery)
)

colnames(techcost.d)[4:5] <- c('FOM','VOM')

saveRDS(techcost.d, file = 'data_out/techcost_df_curated.RDS' )
write.csv(techcost.d, file = 'data_out/techcost_df_curated.csv' )



#### VISUAL INSPECTION ------------------------------------------------------


## checking cost data disaggregate
restech <- unique(techcost.d$technology)
par(mfrow = c(2,2))
for (i in 1:length(restech)) {
  
  
  dplot <- techcost.d[techcost.d$technology == restech[i],]
  
  
  plot(dplot$years, dplot$CAPEX, ylim = range(dplot$CAPEX, na.rm = T), ylab = 'CAPEX', xlab = '', main= restech[i])
  plot(dplot$years, dplot$FOM, ylim = range(dplot$FOM, na.rm = T), ylab = 'FOM', xlab = '', main= restech[i])
  plot(dplot$years, dplot$VOM, ylim = range(dplot$VOM, na.rm = T), ylab = 'VOM', xlab = '', main= restech[i])
  plot(dplot$years, dplot$Fuel, ylim = range(dplot$Fuel, na.rm = T), ylab = 'Fuel', xlab = '', main= restech[i])

  
  
}
