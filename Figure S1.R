library(tidyverse)
library(Seurat)
library(openxlsx)

dat <- read.xlsx("G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/source_data_Figure S1.xlsx")
dat$Best.tumor.response <- factor(dat$Best.tumor.response,levels = c('PR','SD','PD'))
pdf('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/Figure S1.pdf',height = 3,width = 5.5)
ggplot(data=dat,aes(x = patients,y = change))+
  geom_bar(stat = 'identity',aes(fill = Best.tumor.response, color = Best.tumor.response),
           color = 'black',size = 0.1)+
  theme_classic()+
  geom_text(aes(label = label))+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-30,20),color = 'grey',linetype = 'longdash')+
  theme(axis.text = element_text(color = 'black'),axis.text.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank())+
  labs(x = 'Patients',y = 'Change from baseline, %')+
  scale_fill_manual(values = c('#f4c56b','#a7ccce','#d297c2'))+
  scale_color_manual(values = c('#f4c56b','#a7ccce','#d297c2'))
dev.off()
