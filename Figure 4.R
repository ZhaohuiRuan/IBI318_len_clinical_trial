library(tidyverse)
library(Seurat)

#######------------Figure 4A
dat <- read_rds('dalei.rds')
color_used <- read_rds('color_used.rds')
dat$cluster_standard <- factor(dat$cluster_standard,levels = c("Cancer cells", "Club cells", 
                                                               "Alveolar epithelial cells", "Ciliated cells", "Endothelial cells", 
                                                               "Fibroblasts", "Mural cells", "B cells", "Plasma cells", "T cells", 
                                                               "Neutrophils", "Mast cells", "Macrophages", "Monocytes", "cDC1s", 
                                                               "cDC2s", "pDCs", 
                                                               "Osteoclasts")
)
png('main_celltype_dimplot.png',res = 300,height = 5,width = 7.3,units = 'in')
DimPlot(dat,group.by = 'cluster_standard',cols = color_used,raster = F)+NoAxes()+labs(title = '')
dev.off()

#######------------Figure 4B
plot_dat <- dat@meta.data[,c('cluster_standard','group')]
plot_dat$group <- plyr::mapvalues(plot_dat$group %>% as.character(),
                                  from = c("PD", "SD", "PR"),
                                  to = c("NR", "NR", "R"))
plot_dat$group <- factor(plot_dat$group %>% as.character(),levels = c('NR','R'))

summary_data <- plot_dat %>%
  group_by(group, cluster_standard) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

summary_data$cluster_standard <- factor(summary_data$cluster_standard %>% as.character(),
                                        levels = c("Cancer cells", "Club cells", 
                                                   "Alveolar epithelial cells", "Ciliated cells", "Endothelial cells", 
                                                   "Fibroblasts", "Mural cells", "B cells", "Plasma cells", "T cells", 
                                                   "Neutrophils", "Mast cells", "Macrophages", "Monocytes", "cDC1s", 
                                                   "cDC2s", "pDCs", 
                                                   "Osteoclasts"))
pdf('Figure 4B.pdf',height = 4,width = 3.5)
ggplot(summary_data, aes(x = group, y = proportion, fill = cluster_standard)) +
  geom_bar(stat = "identity") +
  theme_classic()+
  labs(x = NULL, y = "Proportion", fill = "Cell types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = 'black'))+
  theme(axis.text.y = element_text(color = 'black'))+
  scale_fill_manual(values = color_used)+
  scale_y_continuous(labels = scales::percent_format())
dev.off()

#######------------Figure 4C
Idents(dat) <- 'id'
plot_dat_cd274 <- DotPlot(dat,features = c('CD274'))$data
plot_dat_pdcd1 <- DotPlot(dat,features = c('PDCD1'))$data

plot_dat_cd274 <- plot_dat_cd274 %>% arrange(desc(pct.exp))
plot_dat_cd274$id <- factor(plot_dat_cd274$id %>% as.character(),levels = plot_dat_cd274$id %>% as.character() %>% rev)
plot_dat_pdcd1 <- plot_dat_pdcd1 %>% arrange(desc(pct.exp))
plot_dat_pdcd1$id <- factor(plot_dat_pdcd1$id %>% as.character(),levels = plot_dat_pdcd1$id %>% as.character() %>% rev)
plot_dat <- rbind(plot_dat_pdcd1,plot_dat_cd274)
plot_dat$id <- factor(plot_dat$id %>% as.character(),
                      levels = c("T cells",  "Myeloid cells", 
                                 "pDCs", 
                                 "Cancer cells", "Club cells","Alveolar epithelial cells", "Ciliated cells","Mural cells",
                                 "Fibroblasts", "Endothelial cells", 
                                 "B cells","Plasma cells" 
                      ) %>% rev)
library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(50) %>% rev()
pdf('Figure 4C.pdf',height = 3,width = 3.5)
ggplot(data = plot_dat,aes(x = features.plot,y = id))+
  geom_point(aes(size = pct.exp,color = avg.exp.scaled))+
  # scale_fill_manual()+
  theme_classic()+
  scale_color_gradientn(colors = cols)+
  theme(axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black',angle = 45,hjust = 1,vjust = 1)
  )+
  labs(x = '',y = '')
dev.off()

#######------------Figure 4D
tcell <- read_rds('TCells.rds')
color_used <- c("#2F70A1", "#3C83A9", "#4C96B0", "#65A6B4", "#82B5BB", "#A2C5C5", 
                "#AAD1D2", "#AADCE0",
                "#FFE6B7", "#FFD98D", "#FDCA6B", "#F9B45E", "#F4A053", "#F08E49", 
                "#EB784C", "#E76254",
                "#41045A", "#620366", "#840978", "#A62890", "#C14DA8", "#D171BB", 
                "#DC9BCD", "#E3C0DB"
)
png('Figure 4D.png',res = 300,height = 5,width =7,units = 'in')
DimPlot(my,group.by = 'cluster_standard',raster = F,cols = color_used)+labs(title = '')+NoAxes()
dev.off()

#######------------Figure 4E
my <- read_rds('MyeliodCells.rds')
color_used <- c("#D9636C", "#FF9898","#CFC734", 
                "#B4BF3A",  "#9EB438",  "#739F35", "#5E9432", "#4C8831", 
                "#3B7D31",
                "#EBCF2E", 
                "#7FCDBBFF", "#41B6C4FF",
                "#EDF8B1FF"
)
png('Figure 4E.png',res = 300,height = 5,width =7,units = 'in')
DimPlot(my,group.by = 'cluster_standard',raster = F,cols = color_used)+labs(title = '')+NoAxes()
dev.off()

#######------------Figure 4F
dat <- read.xlsx("source_data_Figure 4F.xlsx")
dat$sig <- ifelse((dat$p_val_adj<0.05) & (dat$avg_log2FC>0.2),'*','ns')
dat$cluster=as.character(dat$cluster)
dat$cluster = factor(dat$cluster %>% as.character(),
                     levels = c( "CancerCells","AlveolarEpi", "CiliatedCells", "ClubCells", "Fibroblasts", "MuralCells", "ECs", 
                                "BCells", "PlasmaCells", 
                                "CD4_NaiveT_SELL", "CD4_Tactive_CD69", "CD4_Tfh_CCL4", "CD4_Tfh_CXCL13","CD4_Treg_FOXP3", 
                                "CD4_Treg_TNFRSF4","CD4_IFI6", "CD4_ProliferatingT_MKI67", 
                                "CD8_Tactive_FOS", "CD8_Teff_GZMK", "CD8_Trm_ZNF683", "CD8_IFIT2",  "CD8_Tex_TNFRSF9","CD8_Tex_CTLA4",  "CD8_ProliferatingT_MKI67", 
                                "NK_KLRC1", "NKT_TRGC2","NK_CX3CR1","NK_FGFBP2", "ILC3_IL4I1", "NK_XCL1", "NK_TNF",
                                "Monocytes_FCN1",  "Macrophages_CXCL10", "Macrophages_CXCL1", 
                                "Macrophages_SELENOP", "Macrophages_SPP1","Macrophages_CHI3L1", "Macrophages_MKI67", 
                                "Macrophages_FABP4", "Neutrophils_VEGFA","Neutrophils_CXCR2",  "cDC1", "cDC2",  
                                "Osteoclasts") %>% rev)
pdf('Figure 4F.pdf',height = 4,width = 8)
ggplot(data = dat,aes(y = cluster,x = avg_log2FC))+
  geom_vline(xintercept = 0,linetype = 'dashed',color = 'black')+
  geom_segment(aes(xend = avg_log2FC,x = 0,yend = cluster,color = sig),size = 0.5)+
  geom_point(aes(color = sig))+
  theme_classic()+
  facet_wrap(gene~., strip.position = "top",nrow = 2)+
  theme(axis.text = element_text(color = 'black',size = 8),
        axis.text.x = element_text(color = 'black',size = 8,angle = 90,hjust = 1,vjust = 0.5),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        strip.background = element_blank())+
  scale_color_manual(values = c('#bebebe','#ad3d27') %>% rev)+
  coord_flip()+
  labs(y = '',x = 'avg_log2FC')
dev.off()
