rm(list=ls())

library(ggplot2); theme_set(theme_bw())

# take the 3 relevant scenarios
# data distinguishes intra- and inter-regional lines
allcsvs <- c("tran_mi_out_midcase.csv",  "tran_mi_out_95by2050.csv", "tran_mi_out_95by2035.csv",
             "tran_out_midcase.csv", "tran_out_midcase95by2050.csv", "tran_out_midcase95by2035.csv" )
scenarios <- c('Reference', '95% by 2050', '95% by 2035')

scenario.df <- data.frame(files =allcsvs, scenarios = rep(scenarios,2))


filelist <- list()
for (i in 1:length(allcsvs)) {
  
  filelist[[i]] <- read.csv(paste0('data_other/transmission data/', scenario.df[i,1]) )
  filelist[[i]]$scenario <- scenario.df[i,2]
  
}


df1 <- do.call(rbind, filelist[1:3]) # within region transmission lines
df2 <- do.call(rbind, filelist[4:6]) # across region transmission lines


# aggregate data
df1.agg <- aggregate( data.frame(capacity1=df1$Val), by = list(years=df1$Dim1, scenario = df1$scenario), sum ) # sum DC and AC,
df2.agg <- aggregate( data.frame(capacity2=df2$Val), by = list(years=df2$Dim3, scenario = df2$scenario), sum ) # sum DC and AC, also sum all p's

dfmerge <- merge(df1.agg, df2.agg, by = c('years', 'scenario'))

df <- cbind.data.frame(dfmerge[,c(1,2)], capacity = dfmerge$capacity1+dfmerge$capacity2)

write.csv(df, 'data_out/transmission_capacity.csv')


ggplot(df, aes(x=years, y=capacity/10^6, color= scenario)) +
  geom_line(size = 1.25) +
  guides(color = guide_legend(reverse=TRUE), fill = guide_legend(reverse = T)) +
  theme(
    axis.title = element_text(size = 22),
    axis.text.x = element_text(size = 20, hjust = 0.6),
    axis.text.y = element_text(size = 22),
    axis.ticks.length.x = unit(.3,'cm'),
    legend.text = element_text(size = 22),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(1,"cm"),
    legend.position = c(0,1),
    legend.justification = c(0,1),
    panel.spacing = unit(6, "mm"),
    strip.background = element_rect(fill = 'white', color = 'grey'),
    strip.text = element_text(size = 22, color = 'black') ) + 
  xlab('') + ylab('Transmission lines (TW-mi)') 
