library("ggplot2")
library(gridExtra)
##############
###Figure 3###
##############

ggplotDataCombinedAdditiveNoise<-read.table("Data/ggplotData_AdditiveNoise.txt",header=T)
ggplotDataCombinedAdditiveNoise.Spearman<-ggplotDataCombinedAdditiveNoise[which(ggplotDataCombinedAdditiveNoise$expr.distance=="Spearman"),]
ggplotDataCombinedAdditiveNoise.Spearman$status<-factor(ggplotDataCombinedAdditiveNoise.Spearman$status,levels=c("0","10","20","30","40","50"))

fontsize=22
ggplotDataCombinedAdditiveNoise.Fig3<-ggplotDataCombinedAdditiveNoise.Spearman[which(ggplotDataCombinedAdditiveNoise.Spearman$ontology=="Cosine"),]
ggplotDataCombinedAdditiveNoise.Fig3<-ggplotDataCombinedAdditiveNoise.Fig3[which(ggplotDataCombinedAdditiveNoise.Fig3$similarity.measure=="Spearman"),]
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
fig3a<-ggplot2::ggplot(ggplotDataCombinedAdditiveNoise.Fig3,aes(x=status,y=correlation,fill=status,group=status))+
  geom_boxplot()+
  #geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  #stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  theme_bw(fontsize)+
  labs(fill=c("Noise[%]"))+
  xlab("Samples with gaussian noise [%]")+
  ylab("Spearman correlation")+
  theme(legend.position = "none")+
  ggtitle("(a)")+
  scale_fill_manual(values=cbPalette)
fig3a


ggplotDataMeanVariaion<-read.table("Data/ggplotData_Noise_Intensity.txt",header=T)
ggplotDataMeanVariaion.Spearman<-ggplotDataMeanVariaion[which(ggplotDataMeanVariaion$expr.distance=="Spearman"),]
ggplotDataMeanVariaion.Spearman$status<-factor(ggplotDataMeanVariaion.Spearman$status,levels=c("0","5","10","20","30"))

ggplotDataMeanVariaion.Spearman.Fig3<-ggplotDataMeanVariaion.Spearman[which(ggplotDataMeanVariaion.Spearman$ontology=="Cosine"),]
ggplotDataMeanVariaion.Spearman.Fig3<-ggplotDataMeanVariaion.Spearman.Fig3[which(ggplotDataMeanVariaion.Spearman.Fig3$similarity.measure=="Spearman"),]

cbPalette2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
fig3b<-ggplot2::ggplot(ggplotDataMeanVariaion.Spearman.Fig3,aes(x=status,y=correlation,fill=status))+
  geom_boxplot()+
  # geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  # stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  theme_bw(fontsize)+
  labs(fill=c("Noise[%]"))+
  xlab("Mean of gaussian noise")+
  ylab("Spearman correlation")+
  theme(legend.position = "none")+
  ggtitle("(b)")+
  scale_fill_manual(values=cbPalette2)
fig3b

pdf("Figure3.pdf",width=14,height=6)
grid.arrange(fig3a,fig3b,layout_matrix=cbind(1,2))
dev.off()

svg("Figure3.svg",width=14,height=6)
grid.arrange(fig3a,fig3b,layout_matrix=cbind(1,2))
dev.off()



########################################
###Corresponding supplementary figure###
########################################
fontsize=15
sup_fig3a<-ggplot2::ggplot(ggplotDataCombinedAdditiveNoise.Spearman,aes(x=status,y=correlation,fill=status))+
  # geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  # stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  geom_boxplot()+
  theme_bw(fontsize)+
  facet_grid(ontology ~similarity.measure)+
  labs(fill=c("Noise[%]"))+
  xlab("Percentage of added noise")+
  ylab("Correlation")+
  theme(legend.position="none")+
  ggtitle("(a)")+
  scale_fill_manual(values=cbPalette)
sup_fig3a

sup_fig3b<-ggplot2::ggplot(ggplotDataCombinedAdditiveNoise.Spearman,aes(x=tissues,y=correlation,fill=status))+
  # geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  # stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  geom_boxplot()+
  theme_bw(fontsize)+
  facet_grid(ontology ~similarity.measure)+
  labs(fill=c("Noise[%]"))+
  xlab("Tisues")+
  ylab("Correlation")+
  theme(legend.key.height = unit(1.5,"cm"))+
  ggtitle("(b)")+
  scale_fill_manual(values=cbPalette)
sup_fig3b

sup_fig3c<-ggplot2::ggplot(ggplotDataMeanVariaion.Spearman,aes(x=status,y=correlation,fill=status))+
  geom_boxplot()+
  # geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  # stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  theme_bw(fontsize)+
  facet_grid(ontology ~similarity.measure)+
  xlab("Mean of added gaussian noise")+
  ylab("Correlation")+
  theme(legend.position="none")+
  ggtitle("(c)")+
  scale_fill_manual(values=cbPalette2)
sup_fig3c

sup_fig3d<-ggplot2::ggplot(ggplotDataMeanVariaion.Spearman,aes(x=tissues,y=correlation,fill=status))+
  geom_boxplot()+
  # geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  # stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  theme_bw(fontsize)+
  facet_grid(ontology ~similarity.measure)+
  labs(fill=c("Mean of noise"))+
  xlab("Tisues")+
  ylab("Correlation")+
  theme(legend.key.height = unit(1.5,"cm"))+
  ggtitle("(d)")+
  scale_fill_manual(values=cbPalette2)
sup_fig3d

pdf("Supplementary_fig_for_main_fig_3.pdf",width=20,height=10)
grid.arrange(sup_fig3a,sup_fig3b,sup_fig3c,sup_fig3d,layout_matrix=cbind(c(1,1,1,3,3,3),c(2,2,2,4,4,4),c(2,2,2,4,4,4)))
dev.off()

svg("Supplementary_fig_for_main_fig_3.svg",width=20,height=10)
grid.arrange(sup_fig3a,sup_fig3b,sup_fig3c,sup_fig3d,layout_matrix=cbind(c(1,1,1,3,3,3),c(2,2,2,4,4,4),c(2,2,2,4,4,4)))
dev.off()
