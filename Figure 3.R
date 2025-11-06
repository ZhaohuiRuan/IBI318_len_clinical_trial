library(tidyverse)
library(Seurat)
library(openxlsx)
library(survSAKK)
library(survival)



#######-----------Figure 3B
plot_data <- read.xlsx('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/source_data_Figure 3B.xlsx')
fit = survfit(Surv(PFS_months,PFS_status)~1,data = plot_data)
pdf('Figure 3B.pdf',height =4.5,width = 5)
survSAKK::surv.plot(fit=fit,censoring.mark = T,
         xlab="Time since treatment start (months)",lwd = 2,
         ylab="Progression−free survival\n(% of patients )",
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


#######-----------Figure 3C
plot_data <- read.xlsx('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/source_data_Figure 3C.xlsx')
fit <- survfit(Surv(OS_months,OS_status)~1,data = plot_data)
cex_size <- 1
pdf('Figure 3C.pdf',height =4.5,width = 5)
survSAKK::surv.plot(fit=fit,censoring.mark = T,
                    xlab="Time since treatment start (months)",lwd = 2,
                    ylab="Overall survival\n(% of patients )",
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
