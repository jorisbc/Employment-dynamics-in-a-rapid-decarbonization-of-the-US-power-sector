allcsvs <- list.files('data_other/transmission data/')[-1]

filelist <- list()
for (i in 1:length(allcsvs)) {
  
filelist[[i]] <- read.csv(paste0('data_other/transmission data/', allcsvs[i]) )

filelist[[i]]$file <- allcsvs[i]

}


k <- 1
head( filelist[[k]] )
allcsvs[k]
k <- k+1

sapply(filelist, nrow)


lapply(filelist[6:10], function(x) aggregate(x$Val, by = list(x$Dim3, x$Dim4), sum) )



library(ggplot2); theme_set(theme_bw())


df1 <- do.call(rbind, filelist[1:5])

ggplot(df1, aes(x=Dim1, y=Val/10^6, color= file)) +
  geom_line() +
  facet_wrap('Dim2', scales = 'free_y') +
  theme(
    legend.position = c(0,1),
    legend.justification = c(0,1),
    legend.background = element_blank()
  ) +
  ylab('TW-mi')

#sum DC and AC
df1.agg <- aggregate( data.frame(Val=df1$Val), by = list(Dim1=df1$Dim1, file = df1$file), sum )
df1.agg$Valnorm <- df1.agg$Val - 147383540 # subtracting 2020 value

ggplot(df1.agg, aes(x=Dim1, y=Valnorm/10^6, color= file)) +
  geom_line() +
  theme(
    legend.position = c(0,1),
    legend.justification = c(0,1),
    legend.background = element_blank()
  ) +
  ylab('TW-mi')


## the cases corresponding the the 3 fat dashed lines in Fig 8 of NREL Standard Scenarios Report 2021 seem to be 
## tran_mi_out_95by2035.csv
## tran_mi_out_95by2050.csv
## tran_mi_out_midcase.csv -> no new policy


## also look at the intra-regional data

df2 <- do.call(rbind, filelist[6:10])  # note these files also include p's for regions
df2.agg <- aggregate( data.frame(Val=df2$Val), by = list(Dim3=df2$Dim3, file = df2$file), sum ) # sum DC and AC, also sum all p's

df2.agg$Valnorm <- df2.agg$Val - 81722 # subtracting 2020 value

ggplot(df2.agg, aes(x=Dim3, y=Valnorm/10^6, color= file)) +
  geom_line() +
  theme(
    legend.position = c(0,1),
    legend.justification = c(0,1),
    legend.background = element_blank()
  ) +
  ylab('TW-mi')


## note that the intra-regional data is very small 0.5 TW-mi in 2020
## This data has to be summed up




