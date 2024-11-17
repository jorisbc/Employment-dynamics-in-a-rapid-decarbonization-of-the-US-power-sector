rm(list=ls())
library(reshape2)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)


res <- read.csv('data_other/alldata_io_tech.csv')[,-1]
res$technology <- factor(res$technology, levels = rev(c('Coal', 'Gas', 'Nuclear', 'Hydro', 'Geo', 'Bio', 'Solar', 'Wind', 'Batteries')) )
res$scenario <- factor( res$scenario, levels = c('Reference', '95% by 2050', '95% by 2035') )

my_colors <- rev(c('black', 'grey', 'purple', 'blue', 'brown', 'lightgreen', 'yellow', 'darkgreen', 'orange') )

ylim <- range( c(aggregate(res$annual.capex, by = list(res$years, res$scenario), sum)$x,
                 aggregate(res$annual.opex, by = list(res$years, res$scenario), sum)$x ) )/10^6

fig.opex <- ggplot(res, aes(x=years, color = technology, fill = technology, y=annual.opex/10^6)) +
  theme_light() +
  geom_area() +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
#  geom_vline(xintercept = 2020, linetype = 2) + 
  facet_wrap('scenario')  + 
  guides(color = guide_legend(reverse=TRUE), fill = guide_legend(reverse = T)) + 
# scale_y_continuous(limits = ylim ) +
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
    panel.spacing = unit(6, "mm"),  
    legend.position = 'bottom',
    strip.background = element_rect(fill = 'white', color = 'grey'),
    strip.text = element_text(size = 22, color = 'black')
  ) + xlab('') + ylab('Opex [bn USD]')

fig.capex <- ggplot(res, aes(x=years, color = technology, fill = technology, y=annual.capex/10^6)) +
  theme_light() +
  geom_area() +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
#  geom_vline(xintercept = 2020, linetype = 2) + 
  facet_wrap('scenario')  + 
#  scale_y_continuous(limits = ylim ) +
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
  ) + xlab('') + ylab('Capex [bn USD]')




fig1a <- fig.capex + theme(legend.position = 'none')
fig1b <- fig.opex + theme(strip.background = element_blank(),
                      strip.text = element_blank())


pdf('figs/cost_figs.pdf', width = 16, height = 11)
grid.arrange(fig1a, fig1b, ncol = 1, heights = c(1,1.15))
dev.off()
system('open "figs/cost_figs.pdf"')


resnew <- res[res$years %in% 2020:2050,]
aggregate(resnew$annual.capex, by = list(scenario = resnew$scenario), function(x) sum(x,na.rm=T)/10^6) # 1 to 1.9 trillion cumulatively (30-55 bn annually)
aggregate(resnew$annual.opex, by = list(scenario = resnew$scenario), function(x) sum(x,na.rm=T)/10^6) # 2.4 trillion cumulatively


# note much higher investments in Princeton net zero report -> but includes many more cateogries