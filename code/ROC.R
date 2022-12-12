
roc_pea <- roc(roc$group, roc$pearson_p,print.auc =TRUE,
               plot=TRUE,ci = TRUE,smooth = TRUE,
               print.thres="best")
auc_pea <- roc_pea$auc
AUC_CI_pea <- data.frame(ci(roc_pea, method="bootstrap")) 
ret <- c("specificity", "sensitivity", "npv", "ppv")
auc_pea <- paste0(round(auc_pea, 3), " (", round(AUC_CI_pea[1,], 3), ", ", round(AUC_CI_pea[3,], 3),")")

roc_spe <- roc(roc$group, roc$spearman,print.auc =TRUE,
               plot=TRUE,ci = TRUE,smooth = TRUE,
               print.thres="best")
auc_spe <- roc_spe$auc
AUC_CI_spe <- data.frame(ci(roc_spe, method="bootstrap")) 
ret <- c("specificity", "sensitivity", "npv", "ppv")
auc_spe <- paste0(round(auc_spe, 3), " (", round(AUC_CI_spe[1,], 3), ", ", round(AUC_CI_spe[3,], 3),")")

roc_zi <- roc(roc$group, roc$count,print.auc =TRUE,
              plot=TRUE,ci = TRUE,smooth = TRUE,
              print.thres="best")
auc_zi <- roc_zi$auc
AUC_CI_zi <- data.frame(ci(roc_zi, method="bootstrap")) 
ret <- c("specificity", "sensitivity", "npv", "ppv")
auc_zi <- paste0(round(auc_zi, 3), " (", round(AUC_CI_zi[1,], 3), ", ", round(AUC_CI_zi[3,], 3),")")

roc_mi <- roc(roc$group, roc$MI_v1,print.auc =TRUE,
              plot=TRUE,ci = TRUE,smooth = TRUE,
              print.thres="best")
auc_mi <- roc_mi$auc
AUC_CI_mi <- data.frame(ci(roc_mi, method="bootstrap")) 
ret <- c("specificity", "sensitivity", "npv", "ppv")
auc_mi <- paste0(round(auc_mi, 3), " (", round(AUC_CI_mi[1,], 3), ", ", round(AUC_CI_mi[3,], 3),")")

roc_mic <- roc(roc$group, roc$mic_pvalue1,print.auc =TRUE,
               plot=TRUE,ci = TRUE,smooth = TRUE,
               print.thres="best")
auc_mic <- roc_mic$auc
AUC_CI_mic <- data.frame(ci(roc_mic, method="bootstrap")) 
ret <- c("specificity", "sensitivity", "npv", "ppv")
auc_mic <- paste0(round(auc_mic, 3), " (", round(AUC_CI_mic[1,], 3), ", ", round(AUC_CI_mic[3,], 3),")")
ggroc(list( MI = roc_mi,MIC = roc_mic,Spearman = roc_spe,Pearson = roc_pea,ZI = roc_zi),  alpha = 1,aes=c("colour"), linetype = 1, size = 2) + 
  #ggtitle("ROC") +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), colour = "black",linetype=2, size = 1.5) +
  scale_color_jama() + 
  labs(colour = "Method")+
  theme_classic() + 
  # scale_shape_manual() +
  theme(panel.grid =element_blank()) +
  xlab("Specificity") + ylab("Sensitivity") +
  #scale_colour_brewer(palette="RdGy") +
  annotate("text", x=0.25, y=0.18, label=paste0("Spearman AUC: ", round(roc_spe$auc, 3), "(", round(roc_spe$ci[1], 3), ", ", round(roc_spe$ci[3], 3),")"),size=6, fontface="bold")+
  annotate("text", x=0.25, y=0.02, label=paste0("ZI AUC: ", round(roc_zi$auc, 3), "(", round(roc_zi$ci[1], 3), ", ", round(roc_zi$ci[3], 3),")"),size=6, fontface="bold")+
  annotate("text", x=0.25, y=0.1, label=paste0("Pearson AUC: ", round(roc_pea$auc, 3), "(", round(roc_pea$ci[1], 3), ", ", round(roc_pea$ci[3], 3),")"),size=6, fontface="bold")+
  annotate("text", x=0.25, y=0.26, label=paste0("MIC AUC: ", round(roc_mic$auc, 3), "(", round(roc_mic$ci[1], 3), ", ", round(roc_mic$ci[3], 3),")"),size=6, fontface="bold")+
  annotate("text", x=0.25, y=0.34, label=paste0("MI AUC: ", round(roc_mi$auc, 3), "(", round(roc_mi$ci[1], 3), ", ", round(roc_mi$ci[3], 3),")"),size=6, fontface="bold")+
  theme(axis.title.x =element_text(size=15,face = "bold"), 
        axis.title.y=element_text(size=15,face = "bold"),#坐标轴文字大小
        legend.text = element_text(face = "bold", size = 15),#图例文字
        legend.title = element_text("Method",face = "bold", size = 15),
        #axis.line = element_blank(),
        axis.line.y=element_line(linetype="solid",color="black",size=1.5),#坐标轴粗细和文字大小
        axis.line.x=element_line(linetype="solid",color="black",size=1.5),#坐标轴粗细和文字大小
        axis.text.y=element_text(size = 15,face = "bold",color="black"),#坐标轴粗细和文字大小
        axis.text.x=element_text(size = 15,face = "bold",color="black"),
        axis.ticks.x=element_line(color="black",size=1.5,lineend = "round"),
        axis.ticks.y=element_line(color="black",size=1.5,lineend = "round")
  )#坐标轴粗细和文字大小
ggsave("~/roc_linear_1005.png", height = 10, width = 12, dpi = 600)
