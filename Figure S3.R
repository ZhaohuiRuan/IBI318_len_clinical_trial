library(tidyverse)
library(Seurat)
library(openxlsx)
library(forestploter)

get_forest_plot = function(relist){
  tm <- forest_theme(base_size = 9,
                     refline_lwd = 0.5,
                     refline_lty = "dashed",
                     refline_col = "grey20",
                     footnote_cex = 0.8,
                     footnote_fontface = "italic",
                     footnote_col = "grey30",
                     line_spacing = 0.5,
                     core = list(bg_params = list(fill = c("#FFFFFF","#f5f7f6"), col=NA))
  )
  relist$forest = '                                           '
  windowsFonts(myFont1 = windowsFont("Times New Roman"))
  p <- forest(relist[, c('subgroup_name',"subgroup", "No_of_patients",'forest', "HR_CI_label")],
              est = relist$HR %>% as.numeric(),
              lower = relist$HR_CI_lower %>% as.numeric(),
              upper = relist$HR_CI_upper %>% as.numeric(),
              ci_column =4,ref_line = 1,sizes = 0.8,
              xlim = c(0, 6),
              ticks_at = c(0, 1, 2,3,4,5,6),theme = tm)
  return(p)
}

uni_dat <- read.xlsx('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/source_data_Figure S3A.xlsx',sheet = 1)
uni_p <- get_forest_plot(uni_dat)
multi_dat <- read.xlsx('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/source_data_Figure S3A.xlsx',sheet = 2)
multi_p <- get_forest_plot(multi_dat)


pdf('Figure S3A.pdf',height = 15,width = 20)
cowplot::plot_grid(uni_p,multi_p) %>% print()
dev.off()


uni_dat <- read.xlsx('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/source_data_Figure S3B.xlsx',sheet = 1)
uni_p <- get_forest_plot(uni_dat)
multi_dat <- read.xlsx('G:/湘雅/IBI318临床文章/IBI318_小NC_R1_251015/Source data/source_data_Figure S3B.xlsx',sheet = 2)
multi_p <- get_forest_plot(multi_dat)


pdf('Figure S3B.pdf',height = 15,width = 20)
cowplot::plot_grid(uni_p,multi_p) %>% print()
dev.off()
