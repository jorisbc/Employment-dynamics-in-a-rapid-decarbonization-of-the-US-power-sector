rm(list=ls())
library(reshape2)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(readxl)
library(zoo) # moving average function


#### BASIC DATA --------------------------------------------------------------


# the 3 scenarios considered
scenarios <- c('Reference', '95% by 2050', '95% by 2035')

# crosswalk data
tech.crosswalk <- as.data.frame( read_xlsx('data_other/match_tech_datasets.xlsx') )

# ATB cost data
techcost.d <- read.csv('data_out/techcost_df_curated.csv' )[,-1]

# cambium scenario data
capscenario <- read.csv('data_out/scenarios_cambium_capacity.csv')[,-1]   # MW
genscenario <- read.csv('data_out/scenarios_cambium_generation.csv')[,-1] # GWh
scenario.merge <- merge(capscenario,genscenario, all = T)

# historical data
caphist <- read.csv('data_out/historical_capacity_by_source_long.csv')[,-1]
genhist <- read.csv('data_out/historical_generation_by_source_long.csv')[,-1]

hist.merge <- merge(caphist, genhist, by = c('years','technology'), all = T)

# put together history and scenarios:
resfirst <- rbind.data.frame(
  cbind.data.frame(hist.merge, scenario = scenarios[1] )[,c('years', 'technology', 'scenario', 'capacity', 'generation')],
  cbind.data.frame(hist.merge, scenario = scenarios[2] )[,c('years', 'technology', 'scenario', 'capacity', 'generation')],
  cbind.data.frame(hist.merge, scenario = scenarios[3] )[,c('years', 'technology', 'scenario', 'capacity', 'generation')],
  scenario.merge
) # prelim results

## FIX CSP by hand (in references CSP stays constant at ~1GW instead of at 1.7GW) - same problem for generation
resfirst[resfirst$technology == 'CSP' & resfirst$years %in% 2014:2020,'capacity'] <- resfirst[resfirst$technology == 'CSP' & resfirst$years %in% 2022,'capacity']
resfirst[resfirst$technology == 'CSP' & resfirst$years %in% 2014:2020,'generation'] <- resfirst[resfirst$technology == 'CSP' & resfirst$years %in% 2022,'generation']

## FIX jump in hydro in scenarios (not so pronounced for generation)
resfirst[resfirst$technology == 'Hydro' & resfirst$years %in% 2022:2024,'capacity'] <- resfirst[resfirst$technology == 'Hydro' & resfirst$years %in% 2020,'capacity']

# retirement data
retire <- read.csv('data_out/retirement_data.csv')[,-1]

ressecond <- merge(resfirst, retire, by = c('years', 'technology', 'scenario'),  all = T) # prelim results v2
head(ressecond)


#### ANNUAL DATA INTERPOLATION IN SCENARIOS --------------------------------------------------------------

# interpolate cambium data (reported every second year)
# we linearly interpolate capacity and generation
# we distribute two-year retirements equally over both years
tech.seq <- as.character( unique( capscenario$technology) )
yearsfull <- 2010:2050

scenario.interpolated <- NULL
for (i in 1:length(tech.seq)) {
  for (j in 1:length(scenarios)) {

  subdat <- ressecond[ressecond$technology == tech.seq[i] & ressecond$scenario == scenarios[j],]
  capinter <- approx(x = subdat$years, y = subdat$capacity, xout=seq(2010,2050,by=1)) # linear interpolation as data is only given biannually
  geninter <- approx(x = subdat$years, y = subdat$generation, xout=seq(2010,2050,by=1)) # linear interpolation as data is only given biannually
  
  dattemp <- data.frame(years = capinter$x, technology = tech.seq[i], scenario = scenarios[j], capacity = capinter$y, generation = geninter$y, retirement = NA)
  
  # compute annual change of capacity
  dattemp$cap.change <- c(NA, diff( dattemp$capacity ) )
  
  # we split retirement evenly into single years
  dattemp$retirement <- subdat$retirement[ match(dattemp$years, subdat$years) ]
  
  for (t in 1:(length(yearsfull)-1)) {
    
    retire.temp <- dattemp[dattemp$years == yearsfull[t], 'retirement']
  
    if(is.na(retire.temp)) {
      dattemp[dattemp$years == yearsfull[t], 'retirement'] <- dattemp[dattemp$years == yearsfull[t+1], 'retirement']  # take following year of no value given
    } else {
      dattemp[dattemp$years == yearsfull[t], 'retirement'] <- dattemp[dattemp$years == yearsfull[t], 'retirement']
    }

  }

  dattemp$retirement <- ifelse(is.na(dattemp$retirement), 0, dattemp$retirement)
  
  # investment is the residual
  dattemp$investment <- pmax(dattemp$cap.change + dattemp$retirement, 0)
  
  # putting everything together
  scenario.interpolated <- rbind(scenario.interpolated, dattemp)

  }
  
}


#### THE MOST DISAGGREGATE DATA --------------------------------------------------------------

# merge the scenario data with the unit cost data
res <- merge(scenario.interpolated, techcost.d, by = c('years', 'technology'), all=T)
res$scenario <- factor( res$scenario, levels = c('Reference', '95% by 2050', '95% by 2035') )

# total cost estimates
res$annual.capex <- res$investment * res$CAPEX
res$annual.FOM <- res$capacity * res$FOM
res$annual.VOM <- res$generation * res$VOM
res$annual.Fuel <- res$generation * res$Fuel
res$annual.opex <- res$annual.FOM + res$annual.VOM + res$annual.Fuel

# save the data
write.csv(res, 'data_out/alldata_scenario_tech.csv')
writexl::write_xlsx(res, 'data_out/alldata_scenario_tech.xlsx')

head(res,20)



#### DATA AGGREGATION TO IO TECH CATEGORIES --------------------------------------------------------------


# aggregate technologies for IO
techagg <- tech.crosswalk$io[ match(res$technology, tech.crosswalk$scenarios) ]


# note that unit cost data cannot aggregated by summing (without weights)
resvis <- aggregate(res[,c('capacity', 'generation', 'retirement', 'cap.change', 'investment', 'annual.capex', 'annual.FOM', 'annual.VOM', 'annual.Fuel', 'annual.opex')], 
                    by = list(years = res$years, technology = techagg, scenario = res$scenario), function(x)sum(x))
write.csv(resvis,'data_out/alldata_io_tech_nosmoothing.csv' ) 
writexl::write_xlsx(resvis,'data_out/alldata_io_tech_nosmoothing.xlsx' ) 

head(resvis,20)



# #### SMOOTHING ANNUAL COST DATA --------------------------------------------------------------
# 
# io.tech <- unique(tech.crosswalk$io)
# relcols <- c('annual.capex','annual.opex', 'annual.FOM', 'annual.VOM', 'annual.Fuel')
# resio.raw <- resvis[resvis$years %in% 2019:2050,]
# resio <- NULL
# 
# for (i in 1:length(io.tech)) {
#   for (j in 1:length(scenarios)) {
# 
#     subdat <- resio.raw[resio.raw$technology == io.tech[i] & resio.raw$scenario == scenarios[j], ]
#     
#     restemp <- cbind.data.frame(subdat[,1:12], 
#                      rollapply(subdat[,relcols], width = 5, function(...) {mean(...)}, partial = TRUE))
# 
#     resio <- rbind.data.frame(resio, restemp)
#     
#   }
#   
#   
# }
# 
# write.csv(resio,'data_out/alldata_io_tech.csv' )
# writexl::write_xlsx(resio, 'data_out/alldata_io_tech.xlsx')
# 
# head(resio,20)


#### VISUAL INSPECTIONS ------------------------------------------------------


dexp.vis <- melt( aggregate(resvis[,relcols[-2]], list(years = resvis$years, scenario = resvis$scenario), function(x)sum(x)),
                    id.vars = c('years','scenario') )
dexp.res <- melt( aggregate(res[,relcols[-2]], list(years = res$years, scenario = res$scenario), function(x)sum(x)),
                  id.vars = c('years','scenario') )


ggplot(dexp.res, aes(x=years, color = variable, fill = variable, y=value/10^6)) +
  theme_light() +
  geom_area() +
  # scale_fill_manual(values = my_colors) +
  # scale_color_manual(values = my_colors) +
  facet_wrap('scenario')  +
  #  scale_x_continuous(breaks = seq(2020,2050,10)) +
  guides(color = guide_legend(reverse=TRUE), fill = guide_legend(reverse = T)) +
  theme(
    axis.title = element_text(size = 22),
    axis.text.x = element_text(size = 20, hjust = 0.6),
    axis.text.y = element_text(size = 22),
    axis.ticks.length.x = unit(.3,'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 22),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(1,"cm"),
    legend.position = 'bottom',
    panel.spacing = unit(6, "mm"),
    strip.background = element_rect(fill = 'white', color = 'grey'),
    strip.text = element_text(size = 22, color = 'black')
  ) + xlab('') + ylab('Expenses in bn USD')



## checking cost data disaggregate
restech <- unique(res$technology)
par(mfrow = c(2,2))
for (i in 1:length(restech)) {


  dplot <- res[res$technology == restech[i],]
  
  dplot.s1 <- dplot[dplot$scenario == scenarios[1],]
  dplot.s2 <- dplot[dplot$scenario == scenarios[2],]
  dplot.s3 <- dplot[dplot$scenario == scenarios[3],]
      
  plot(dplot.s1$years, dplot.s1$CAPEX, ylim = range(dplot$CAPEX, na.rm = T), ylab = 'CAPEX', xlab = '', main= restech[i])
  lines(dplot.s2$years, dplot.s2$CAPEX, col = 'blue' )
  lines(dplot.s2$years, dplot.s2$CAPEX, col = 'red' )
  abline(v=2017)
  abline(v=2019)
  
  plot(dplot.s1$years, dplot.s1$FOM, ylim = range(dplot$FOM, na.rm = T), ylab = 'FOM', xlab = '', main= restech[i])
  lines(dplot.s2$years, dplot.s2$FOM, col = 'blue' )
  lines(dplot.s2$years, dplot.s2$FOM, col = 'red' )
  abline(v=2017)
  abline(v=2019)
  
  
  plot(dplot.s1$years, dplot.s1$VOM, ylim = range(dplot$VOM, na.rm = T), ylab = 'VOM', xlab = '', main= restech[i])
  lines(dplot.s2$years, dplot.s2$VOM, col = 'blue' )
  lines(dplot.s2$years, dplot.s2$VOM, col = 'red' )
  abline(v=2017)
  abline(v=2019)
  
  
  plot(dplot.s1$years, dplot.s1$Fuel, ylim = range(dplot$Fuel, na.rm = T), ylab = 'Fuel', xlab = '', main= restech[i])
  lines(dplot.s2$years, dplot.s2$Fuel, col = 'blue' )
  lines(dplot.s2$years, dplot.s2$Fuel, col = 'red' )
  abline(v=2017)
  abline(v=2019)
  
    
}



## checking cap and gen data disaggregate
par(mfrow = c(2,2))
for (i in 1:length(restech)) {
  
  
  dplot <- res[res$technology == restech[i],]
  
  dplot.s1 <- dplot[dplot$scenario == scenarios[1],]
  dplot.s2 <- dplot[dplot$scenario == scenarios[2],]
  dplot.s3 <- dplot[dplot$scenario == scenarios[3],]
  
  plot(dplot.s1$years, dplot.s1$capacity, ylim = range(dplot$capacity, na.rm = T), ylab = 'capacity', xlab = '', main= restech[i])
  lines(dplot.s2$years, dplot.s2$capacity, col = 'blue' )
  lines(dplot.s3$years, dplot.s3$capacity, col = 'red' )
  abline(v=2017)
  abline(v=2019)
  
  plot(dplot.s1$years, dplot.s1$generation, ylim = range(dplot$generation, na.rm = T), ylab = 'generation', xlab = '', main= restech[i])
  lines(dplot.s2$years, dplot.s2$generation, col = 'blue' )
  lines(dplot.s3$years, dplot.s3$generation, col = 'red' )
  abline(v=2017)
  abline(v=2019)
  
  
  plot(dplot.s1$years, dplot.s1$retirement, ylim = range(dplot$retirement, na.rm = T), ylab = 'retirement', xlab = '', main= restech[i])
  lines(dplot.s2$years, dplot.s2$retirement, col = 'blue' )
  lines(dplot.s3$years, dplot.s3$retirement, col = 'red' )
  abline(v=2017)
  abline(v=2019)
  
  
  plot(dplot.s1$years, dplot.s1$cap.change, ylim = range(dplot$cap.change, na.rm = T), ylab = 'cap.change', xlab = '', main= restech[i])
  lines(dplot.s2$years, dplot.s2$cap.change, col = 'blue' )
  lines(dplot.s3$years, dplot.s3$cap.change, col = 'red' )
  abline(v=2017)
  abline(v=2019)
  
  
}


# resio$technology <- factor(resio$technology, levels = rev(c('Coal', 'Gas', 'Nuclear', 'Hydro', 'Geo', 'Bio', 'Solar', 'Wind', 'Batteries')) )
# my_colors <- rev(c('black', 'grey', 'purple', 'blue', 'brown', 'lightgreen', 'yellow', 'darkgreen', 'orange') )
# 
# 
# ggplot(resio, aes(x=years, color = technology, fill = technology, y=annual.opex/10^6)) +
#   theme_light() +
#   geom_area() +
#   scale_fill_manual(values = my_colors) +
#   scale_color_manual(values = my_colors) +
#   facet_wrap('scenario')  + 
#   #  scale_x_continuous(breaks = seq(2020,2050,10)) +
#   guides(color = guide_legend(reverse=TRUE), fill = guide_legend(reverse = T)) + 
#   theme(
#     axis.title = element_text(size = 22),
#     axis.text.x = element_text(size = 20, hjust = 0.6),
#     axis.text.y = element_text(size = 22),
#     axis.ticks.length.x = unit(.3,'cm'),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.text = element_text(size = 22),
#     legend.title = element_blank(),
#     legend.background = element_blank(),
#     legend.key.size = unit(1,"cm"),
#     legend.position = 'bottom',
#     panel.spacing = unit(6, "mm"),  
#     strip.background = element_rect(fill = 'white', color = 'grey'),
#     strip.text = element_text(size = 22, color = 'black')
#   ) + xlab('') + ylab('Opex bn USD')
# 
# ggplot(resio, aes(x=years, color = technology, fill = technology, y=annual.capex/10^6)) +
#   theme_light() +
#   geom_area() +
#   scale_fill_manual(values = my_colors) +
#   scale_color_manual(values = my_colors) +
#   facet_wrap('scenario')  + 
#   #  scale_x_continuous(breaks = seq(2020,2050,10)) +
#   guides(color = guide_legend(reverse=TRUE), fill = guide_legend(reverse = T)) + 
#   theme(
#     axis.title = element_text(size = 22),
#     axis.text.x = element_text(size = 20, hjust = 0.6),
#     axis.text.y = element_text(size = 22),
#     axis.ticks.length.x = unit(.3,'cm'),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.text = element_text(size = 22),
#     legend.title = element_blank(),
#     legend.background = element_blank(),
#     legend.key.size = unit(1,"cm"),
#     legend.position = 'bottom',
#     panel.spacing = unit(6, "mm"),  
#     strip.background = element_rect(fill = 'white', color = 'grey'),
#     strip.text = element_text(size = 22, color = 'black')
#   ) + xlab('') + ylab('Capex bn USD')
# 
# 
# techuni <- unique(res$technology)
# plotlist.gen <- plotlist.cap <- plotlist.retirement <- plotlist.invest <- list()
# for (i in 1:length(techuni)) {
# #  i <- 1
#   # dat <- resfirst[resfirst$technology == techuni[i],]
#   dat <- res[res$technology == techuni[i],]
#   
#   plotlist.cap[[i]] <- ggplot(dat, aes(x=years, y=capacity, color = scenario)) +
#     geom_line() +
#     theme(
#       legend.position = c(0,1),
#       legend.justification = c(0,1),
#       legend.background = element_blank()
#     ) + xlab("") + ylab(techuni[i])
#   
#   
#   plotlist.gen[[i]] <- ggplot(dat, aes(x=years, y=generation, color = scenario)) +
#     geom_line() +
#     theme(
#       legend.position = c(0,1),
#       legend.justification = c(0,1),
#       legend.background = element_blank()
#     ) + xlab("") + ylab(techuni[i])
#   
#   plotlist.retirement[[i]] <- ggplot(dat, aes(x=years, y=retirement, color = scenario)) +
#     geom_line() +
#     theme(
#       legend.position = c(0,1),
#       legend.justification = c(0,1),
#       legend.background = element_blank()
#     ) + xlab("") + ylab(techuni[i])
#   
#   # plotlist.invest[[i]] <- ggplot(dat, aes(x=years, y=invest, color = scenario)) +
#   #   geom_line() +
#   #   theme(
#   #     legend.position = c(0,1),
#   #     legend.justification = c(0,1),
#   #     legend.background = element_blank()
#   #   ) + xlab("") + ylab(techuni[i])
# }
# 
# do.call(grid.arrange, c(plotlist.cap, ncol=4) )
# do.call(grid.arrange, c(plotlist.gen, ncol=4) )
# do.call(grid.arrange, c(plotlist.retirement, ncol=4) )
# #do.call(grid.arrange, c(plotlist.invest, ncol=4) )
