#########-------------xgboost
library(xgboost)
library(tidyverse)
library(Seurat)
library(pROC)
library(openxlsx)
###--- 1. Data loading
load('data_xgb_model.RData')
dtrain <- xgb.DMatrix("IBI318_len_scRNA_seq_dtrain.buffer")
dtrain2 <- xgb.DMatrix("IBI318_len_scRNA_seq_dtrain_PD_L1.buffer")
###--- 2. XGBoost model based on single-cell RNA-seq data
seed_used = 1898
set.seed(seed_used)
cv <- xgb.cv(
  data = dtrain,
  objective = "binary:logistic",
  nrounds = 1000,
  eta = 0.05,
  max_depth = 3,
  subsample = 0.7,
  colsample_bytree = 0.8,
  nfold = 3,
  early_stopping_rounds = 20,
  verbose = 0
)
best_nrounds <- cv$best_iteration

set.seed(seed_used)
xgb_model <- xgboost(
  data = dtrain,
  objective = "binary:logistic",
  nrounds = best_nrounds,
  eta = 0.05,
  max_depth = 3,
  subsample = 0.7,
  colsample_bytree = 0.8,
  verbose = 0
)

###---3. XGBOOST based on PD-L1 TPS
set.seed(1898)
xgb_model_pdl1 <- xgboost(
  data = dtrain2,
  objective = "binary:logistic",
  nrounds = best_nrounds,
  eta = 0.05,
  max_depth = 4,
  subsample = 0.8,
  colsample_bytree = 0.8,
  verbose = 0
)

####-------Figure 5B
pred_prob <- predict(xgb_model, dtrain)
true_label <- getinfo(dtrain, "label")
roc_obj <- roc(response = true_label, predictor = pred_prob)

pred_prob_pdl1 <- predict(xgb_model_pdl1, dtrain2)
true_label_pdl1 <- getinfo(dtrain2, "label")
roc_obj_pdl1 <- roc(response = true_label_pdl1, predictor = pred_prob_pdl1)


pdf('Figure 5B.pdf',height =4,width = 4)
plot(roc_obj, print.auc = TRUE, col = "#1c61b6", main = "ROC Curve Comparison")
lines(roc_obj_pdl1, col = "#008600")
auc1 <- auc(roc_obj)
auc2 <- auc(roc_obj_pdl1)
legend("bottomright", 
       legend = c(paste0("XGBoost: AUC=", round(auc1, 3)),
                  paste0("PD-L1 model: AUC=", round(auc2, 3))),
       col = c("#1c61b6", "#008600"), lwd = 2)
dev.off()


####-------Figure 5C
load('data_for_xgboost_IBI318_len_scRNA_seq.RData')
shape_value <- SHAPforxgboost::shap.prep(xgb_model = xgb_model,X_train = x)
library(SHAPforxgboost)
p <- SHAPforxgboost::shap.plot.summary(shape_value)

shape_value_summary <- shape_value %>% select(variable,mean_value) %>% distinct() %>% as.data.frame()

library(viridis)
p2 <- ggplot(shape_value_summary,aes(y = factor(variable %>% as.character(),
                                                levels = levels(variable) %>% rev),
                                     x = mean_value,fill = mean_value))+
  theme_bw()+
  geom_bar(stat = 'identity')+
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        legend.position = "bottom", legend.title = element_text(size = 10), 
        legend.text = element_text(size = 8), axis.title.x = element_text(size = 10))+
  labs(x = 'Mean Shape value',y = '')+
  viridis::scale_fill_viridis()+
  scale_x_continuous(breaks = seq(min(shape_value_summary$mean_value),max(shape_value_summary$mean_value),length.out = 3) %>% round(2))

pdf('Figure 5C.pdf',height = 7,width = 7)
cowplot::plot_grid(p,p2,align = 'hv',rel_widths = c(1,0.8))
dev.off()

