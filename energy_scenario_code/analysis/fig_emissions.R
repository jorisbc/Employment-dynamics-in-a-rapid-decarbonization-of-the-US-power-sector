rm(list=ls())
library(reshape2)
library(readxl)
library(ggplot2); theme_set(theme_bw())

tech.crosswalk <- as.data.frame( read_xlsx('data_other/match_tech_datasets.xlsx') )

case.ref <- read.csv('D:/Science/data_other/energy/NREL Cambium 2021 scenarios/StdScen21_MidCase_annual_national/StdScen21_MidCase_annual_national.csv', skip = 4)
case.35 <- read.csv('D:/Science/data_other/energy/NREL Cambium 2021 scenarios/StdScen21_MidCase95by2035_annual_national/StdScen21_MidCase95by2035_annual_national.csv', skip = 4)
case.50 <- read.csv('D:/Science/data_other/energy/NREL Cambium 2021 scenarios/StdScen21_MidCase95by2050_annual_national/StdScen21_MidCase95by2050_annual_national.csv', skip = 4)

colnames(case.50$co)

gen_name <- sort( colnames(case.ref)[77:94] )

gen_name


# total_net_co2e # generation emissions - capturing
# total_gen_co2e # generation emissions


plot(case.35$total_net_co2e, ylim = c(0, max(case.35$total_net_co2)) )

netemit.w <- data.frame(
years = case.35$t,
"Reference" = case.ref$total_net_co2e,
"95% by 2035" = case.35$total_net_co2e,
"95% by 2050" = case.50$total_net_co2e,
check.names = FALSE
)

netemit <- melt(netemit.w, id.vars = 'years', variable.name = 'scenario', value.name = 'emission')
emit.hist <- unlist( read.csv("data_other/epa_emissions.csv", row.names = 1)['Total',] ) * 10^6 

# emit.hist.1 <- emit.hist.2 <- emit.hist.3 <- data.frame(years = 1990:2020, scenario = 'Reference', emission = emit.hist)
# emit.hist.2$scenario <- "95% by 2035"
# emit.hist.3$scenario <- "95% by 2050"

emit.hist.1 <- emit.hist.2 <- emit.hist.3 <- data.frame(years = 1990:2020, scenario = 'Empirical', emission = emit.hist)
emit.hist.1 <- rbind(emit.hist.1,  data.frame(years = 2022, scenario = 'Empirical', emission = netemit[1,3]))

emitfull <- rbind(
  emit.hist.1,
  netemit
)

emitfull$scenario <- factor(emitfull$scenario, levels = c('Empirical', "Reference", "95% by 2050", "95% by 2035"))





pdf('figs/emission_figs.pdf', width = 9)
ggplot(emitfull, aes(x=years, color = scenario, fill = scenario, y=emission/10^6)) +
  geom_line(size = 1.2) +
#  geom_hline(yintercept = emitfull[ emitfull$years == 2005, "emission"]/10^6 * 0.05) + 
  scale_color_manual(values = c('black', 'darkgreen', 'blue', "firebrick")) +
  scale_x_continuous(breaks = seq(1990,2050,10)) +
  scale_y_continuous(limits = c(0,max(emitfull$emission)/10^6)) +
  theme(
    axis.title = element_text(size = 22),
    axis.text.x = element_text(size = 20, hjust = 0.6),
    axis.text.y = element_text(size = 22),
    axis.ticks.length.x = unit(.3,'cm'),
    # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.text = element_text(size = 22),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(1,"cm"),
    legend.position = c(1,1),
    legend.justification = c(1,1),
    panel.spacing = unit(6, "mm"),  
    strip.background = element_rect(fill = 'white', color = 'grey'),
    strip.text = element_text(size = 22, color = 'black')
  ) + xlab('') + ylab(bquote(Power~sector~emissions~CO[2]~e~Mt) )
dev.off()
system('open "figs/emission_figs.pdf"')
