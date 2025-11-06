library(survSAKK)
library(survival)
library(openxlsx)
source('G:/湘雅/昼夜节律/ToD_20241030/sur_plot.R')
plot_dat <- read.xlsx('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/source_data_Figure S2.xlsx')
fit = survfit(Surv(DOR_months,DOR_status)~1,data = plot_dat)
cex_size = 1
re = summary(fit)
time_max = max(re$time,na.rm = T)
pdf('Figure S2.pdf',height =4.5,width = 5)
survSAKK::surv.plot(fit=fit,censoring.mark = T,
         xlab="Time since treatment start (months)",lwd = 2,
         ylab="Duration of Response\n(% of patients)",
         xlab.cex=cex_size,
         segment.cex = cex_size,
         censoring.cex = 0.8, 
         ylab.cex=cex_size,main = '',
         xlab.pos=5.5,
         y.unit = 'percent',conf.int = 0.95,
         col = c('#478bc0','#9e313c'),
         time.unit = "month",
         segment.quantile = 0.5,
         segment.annotation = "top",
         risktable.cex = cex_size,
         risktable.title.cex = cex_size,risktable.name.cex = cex_size,
         ylab.pos=3,legend.title.cex = 1,risktable.title ='No. at Risk' 
)
dev.off()
