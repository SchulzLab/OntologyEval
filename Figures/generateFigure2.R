library(ggplot2)
library(ggpubr)
library(gridExtra)

combinedData<-read.table("Data/ggplotData_Spearman_ontology_randomized.txt",header=T,stringsAsFactors = F)
combinedData$random.ontology[which(combinedData$random.ontology == "No")]<-"Original"
combinedData$random.ontology[which(combinedData$random.ontology == "Yes")]<-"Randomized"
combinedData$random.ontology[which(combinedData$random.ontology == "No")]<-"Fixed"
combinedData$random.ontology<-factor(combinedData$random.ontology,levels=c("Original","Randomized","Fixed"))
cbPalette <- c("#999999", "#E69F00", "#56B4E9")
##############
###Figure 2###
##############
fontsize=22
combinedData.Figure2<-combinedData[which(combinedData$ontology=="Cosine"),]
combinedData.Figure2<-combinedData.Figure2[which(combinedData.Figure2$similarity.measure=="Spearman"),]
comp<-list(c("Original","Randomized"),c("Original","Fixed"))
svg("Figure2.svg")
ggplot2::ggplot(combinedData.Figure2,aes(x=random.ontology,y=correlation,fill=random.ontology))+
  geom_boxplot(notch=T)+
  theme_bw(fontsize)+
  stat_compare_means(method="wilcox.test",label="p.signif",comparisons=comp,label.size=3)+
  xlab("")+
  ylab("Spearman correlation")+
  theme(legend.position = "none")+
  scale_fill_manual(values=cbPalette)
dev.off()

pdf("Figure2.pdf")
ggplot2::ggplot(combinedData.Figure2,aes(x=random.ontology,y=correlation,fill=random.ontology))+
  geom_boxplot(notch=T)+
  theme_bw(fontsize)+
  stat_compare_means(method="wilcox.test",label="p.signif",comparisons=comp,label.size=3)+
  xlab("")+
  ylab("Spearman correlation")+
  theme(legend.position = "none")+
  scale_fill_manual(values=cbPalette)
dev.off()

########################################
###Corresponding Supplementary Figure###
########################################
fontsize=20
supFiga<-ggplot2::ggplot(combinedData,aes(x=tissues,y=correlation,fill=random.ontology))+
  geom_boxplot()+theme_bw(fontsize)+
  facet_grid(ontology ~similarity.measure)+
  xlab("Tissue")+
  ylab("Correlation")+
  labs(fill="Ontology")+
  theme(legend.key.height = unit(2,"cm"))+
  ggtitle("(b)")+
  scale_fill_manual(values=cbPalette)
  
supFigb<-ggplot2::ggplot(combinedData,aes(x=random.ontology,y=correlation,fill=random.ontology))+
  facet_grid(ontology ~similarity.measure)+
  geom_boxplot(notch=T)+
  theme_bw(fontsize)+
  stat_compare_means(method="wilcox.test",label="p.signif",comparisons=comp,label.size=3)+
  xlab("")+
  ylab("Correlation")+
  theme(legend.position = "none")+
  ggtitle("(a)")+
  scale_fill_manual(values=cbPalette)

svg("Supplementary_fig_for_main_fig_2.svg",width=24,height=10)
grid.arrange(supFiga,supFigb,layout_matrix=cbind(2,2,2,1,1,1,1,1))
dev.off()

pdf("Supplementary_fig_for_main_fig_2.pdf",width=24,height=10)
grid.arrange(supFiga,supFigb,layout_matrix=cbind(2,2,2,1,1,1,1,1))
dev.off()
