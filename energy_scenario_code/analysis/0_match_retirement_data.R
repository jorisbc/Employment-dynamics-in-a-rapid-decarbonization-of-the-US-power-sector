rm(list=ls())

library(readxl)
library(ggplot2); theme_set(theme_bw())

# retirement data
fast<-read.csv('data_other/retirement_data/ret_ann_95by2035.csv')
slow<-read.csv('data_other/retirement_data/ret_ann_95by2050.csv')
no<-read.csv('data_other/retirement_data/ret_ann_midcase.csv')
# we ignore canada imports

# pvb is PV-Battery hybrids
set1 <- c("pvb1_2","pvb1_3","pvb1_4","pvb1_5", "pvb1_6", "pvb1_7")
plot( aggregate( no[no$Dim1 %in% set1,'Val'], by= list(year = no[no$Dim1 %in% set1,'Dim3']),sum)$x)

set2<- c("upv_2","upv_3","upv_4","upv_5", "upv_6", "upv_7")
plot( aggregate( no[no$Dim1 %in% set2,'Val'], by= list(year = no[no$Dim1 %in% set2,'Dim3']),sum)$x)

set3<- c("dupv_2","dupv_3","dupv_4","dupv_5", "dupv_6", "dupv_7")
plot( aggregate( no[no$Dim1 %in% set3,'Val'], by= list(year = no[no$Dim1 %in% set3,'Dim3']),sum)$x)

unique(fast$Dim1)
unique(slow$Dim1)
unique(no$Dim1)

sort( unique(fast$Dim2) ) # Dim2 is regional code (there should be no double counting, so simply aggregate)
table( data.frame(no$Dim1, substr(no$Dim2,1,1)) ) # only wind and csp are on the s-regions

no[no$Dim3 == 2030,]




par(mfrow=c(4,3))
aha <- aggregate(no$Val, by = list(no$Dim1, no$Dim3), sum)
aha.f <- aggregate(fast$Val, by = list(fast$Dim1, fast$Dim3), sum)
aha.s <- aggregate(slow$Val, by = list(slow$Dim1, slow$Dim3), sum)
ahatech <- unique(aha$Group.1)
for (i in 1:length(ahatech)) {
  #i <- 1
  tempdat <- aha[aha$Group.1 == ahatech[i],]
  tempdat.f <- aha.f[aha.f$Group.1 == ahatech[i],] # fast
  tempdat.s <- aha.s[aha.s$Group.1 == ahatech[i],] # slow
  
  x.range <- range(tempdat.f$Group.2,tempdat.s$Group.2, tempdat$Group.2)
  y.range <- range(tempdat.f$x,tempdat.s$x, tempdat$x)
  
  plot(tempdat$Group.2, tempdat$x, main = ahatech[i], xlab= "", ylab = 'retiremenet MW', type = 'l', xlim = x.range, ylim = y.range) # reference
  lines(tempdat.s$Group.2, tempdat.s$x, xlab= "", ylab = 'retiremenet MW', col = 'red', lty = 2) # slow
  points(tempdat.f$Group.2, tempdat.f$x, xlab= "", ylab = 'retiremenet MW', col = 'blue') # fast
  
}


# aggregating technologies
tech.crosswalk <- as.data.frame( read_xlsx('data_other/match_tech_datasets.xlsx', sheet = 'match_retirement') )


no$technology <- tech.crosswalk$scenarios[ match(no$Dim1, tech.crosswalk$Retirement) ]
slow$technology <- tech.crosswalk$scenarios[ match(slow$Dim1, tech.crosswalk$Retirement) ]
fast$technology <- tech.crosswalk$scenarios[ match(fast$Dim1, tech.crosswalk$Retirement) ]

retfast <- aggregate(data.frame(retirement=fast$Val), by = list(technology = fast$technology, years = fast$Dim3), sum)
retslow <- aggregate(data.frame(retirement=slow$Val), by = list(technology = slow$technology, years = slow$Dim3), sum)
retno <- aggregate(data.frame(retirement=no$Val), by = list(technology = no$technology, years = no$Dim3), sum)


retall <- rbind(
  cbind.data.frame(retfast, scenario = '95% by 2035'),
  cbind.data.frame(retslow, scenario = '95% by 2050'),
  cbind.data.frame(retno, scenario = 'Reference')
)


write.csv(retall, file = 'data_out/retirement_data.csv')


## aggregate further for figure
no$technology <- tech.crosswalk$io[ match(no$Dim1, tech.crosswalk$Retirement) ]
slow$technology <- tech.crosswalk$io[ match(slow$Dim1, tech.crosswalk$Retirement) ]
fast$technology <- tech.crosswalk$io[ match(fast$Dim1, tech.crosswalk$Retirement) ]

retfast <- aggregate(data.frame(retirement=fast$Val), by = list(technology = fast$technology, years = fast$Dim3), sum)
retslow <- aggregate(data.frame(retirement=slow$Val), by = list(technology = slow$technology, years = slow$Dim3), sum)
retno <- aggregate(data.frame(retirement=no$Val), by = list(technology = no$technology, years = no$Dim3), sum)


retall <- rbind(
  cbind.data.frame(retfast, scenario = '95% by 2035'),
  cbind.data.frame(retslow, scenario = '95% by 2050'),
  cbind.data.frame(retno, scenario = 'Reference')
)


ggplot(retall, aes(x=years, color = scenario, fill = scenario, y=retirement)) +
  geom_line(size = 1.2) +
  scale_x_continuous(breaks = seq(1990,2050,10)) +
#  scale_y_continuous(limits = c(0,max(retall$retirement)/10^6)) +
  facet_wrap('technology', scales = 'free_y') +
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
    legend.position = 'bottom',
#    legend.justification = c(1,1),
    panel.spacing = unit(6, "mm"),  
    strip.background = element_rect(fill = 'white', color = 'grey'),
    strip.text = element_text(size = 22, color = 'black')
  ) + xlab('') + ylab('MW')


# 
# 
# # ATB cost data
# techcost.d <- read.csv('data_out/techcost_df_curated.csv' )[,-1]
# 
# # cambium scenario data
# capscenario <- read.csv('data_out/scenarios_cambium_capacity.csv')[,-1]   # MW
# 
# 
# 
# unique(techcost.d$technology)
# unique(capscenario$technology)