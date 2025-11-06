library(tidyverse)
library(Seurat)
library(openxlsx)

plot_dat <- read.xlsx('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/source_data_Figure 2A.xlsx')
plot_dat$ORR_12<- factor(plot_dat$ORR_12,levels = c('PR','SD','PD'))

pdf('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/Figure 2A.pdf',height = 3,width = 5.5)
ggplot(data=plot_dat[!(is.na(plot_dat$delta_12week)),],aes(x = factor(patient_order),y = delta_12week))+
  geom_bar(stat = 'identity',aes(fill = ORR_12, color = ORR_12),color = 'black',size = 0.1)+
  theme_classic()+
  geom_text(aes(label = label))+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-30,20),color = 'grey',linetype = 'longdash')+
  theme(axis.text = element_text(color = 'black'),axis.text.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank())+
  labs(x = 'Patients',y = 'Change from baseline, %')+
  scale_fill_manual(values = c('#f4c56b','#a7ccce','#d297c2'))+
  scale_color_manual(values = c('#f4c56b','#a7ccce','#d297c2'))
dev.off()

plot_dat<- read.xlsx('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/source_data_Figure 2B.xlsx')
pdf('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/Figure 2B.pdf',height = 3,width = 3)
ggplot(plot_dat, aes(x = '12-week ORR', y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Freq*100, 1), "%")),
            position = position_stack(vjust = 0.5), color = "black") +
  theme_classic()+
  theme(axis.text = element_text(color = 'black'),axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  labs(y = 'Percentage',x = '')+
  scale_fill_manual(values = c('#f4c56b','#a7ccce','#d297c2'),name = '12 week ORR')+
  scale_color_manual(values = c('#f4c56b','#a7ccce','#d297c2'),name = '12 week ORR')+scale_y_continuous(labels = scales::percent_format())
dev.off()
