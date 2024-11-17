rm(list=ls())
library(reshape2)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(readxl)

tech.crosswalk <- as.data.frame( read_xlsx('data_other/match_tech_datasets.xlsx') )

res <- read.csv('data_other/alldata_io_tech.csv')[,-1]
res$technology <- factor(res$technology, levels = rev(c('Coal', 'Gas', 'Nuclear', 'Hydro', 'Geo', 'Bio', 'Solar', 'Wind', 'Batteries')) )
res$scenario <- factor( res$scenario, levels = c('Reference', '95% by 2050', '95% by 2035') )

my_colors <- rev(c('black', 'grey', 'purple', 'blue', 'brown', 'lightgreen', 'yellow', 'darkgreen', 'orange') )

pcap <- ggplot(res, aes(x=years, color = technology, fill = technology, y=capacity/10^3)) +
  theme_light() +
  geom_area() +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  geom_vline(xintercept = 2020, linetype = 2) + 
  facet_wrap('scenario')  + 
  guides(color = guide_legend(reverse=TRUE), fill = guide_legend(reverse = T)) + 
  #  scale_x_continuous(breaks = seq(2020,2050,10)) +
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
  ) + xlab('') + ylab('Capacity GW')

pgen <- ggplot(res, aes(x=years, color = technology, fill = technology, y=generation/10^3)) +
  theme_light() +
  geom_area() +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  geom_vline(xintercept = 2020, linetype = 2) + 
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
  ) + xlab('') + ylab('Generation TWh')




fig1a <- pcap + theme(legend.position = 'none')
fig1b <- pgen + theme(strip.background = element_blank(),
             strip.text = element_blank())


pdf('figs/scenario_figs.pdf', width = 16, height = 11)
grid.arrange(fig1a, fig1b, ncol = 1, heights = c(1,1.15))
dev.off()
system('open "figs/scenario_figs.pdf"')
