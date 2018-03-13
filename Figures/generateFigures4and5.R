library(ggplot2)
library(gridExtra)

ggplotDataCombined<-read.table("Data/ggplotData_Noise_Experiment.txt",header=T)
ggplotDataCombined.Spearman<-ggplotDataCombined[which(ggplotDataCombined$expr.distance=="Spearman"),]
ggplotDataCombined.Spearman$status<-factor(ggplotDataCombined.Spearman$status,levels=c("Original","Noise","Combat"))
ggplotDataCombinedFig3a<-ggplotDataCombined.Spearman[which(ggplotDataCombined.Spearman$similarity.measure=="Spearman"),]
ggplotDataCombinedFig3a<-ggplotDataCombinedFig3a[which(ggplotDataCombinedFig3a$ontology=="Cosine"),]
ggplotDataNoLog<-read.table("Data/ggplotDataNoLog.txt",header=T,sep="\t",stringsAsFactors = F)
ggplotDataNoLog<-ggplotDataNoLog[which(ggplotDataNoLog$expr.distance=="Spearman"),]
fontsize=25

##############
###Figure 4###
##############
pcaOriginal<-read.table("Data/PCA_Original.txt",header=T)
pcaNoise<-read.table("Data/PCA_Noise.txt",header=T)
pcaNoise$Batch<-factor(pcaNoise$Batch,levels=c("Original","Noise"))
pcaCombat<-read.table("Data/PCA_Combat.txt",header=T)
pcaCombat$Batch<-factor(pcaCombat$Batch,levels=c("Original","Noise"))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
fontsize=20

pcaPlot12.Original<-ggplot2::ggplot(data=pcaOriginal,aes(x=PC1,y=PC2,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle("(a) Original")+
  theme(legend.position="none")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(20,20))

pcaPlot12.Original

pcaPlot12.noise<-ggplot2::ggplot(data=pcaNoise,aes(x=PC1,y=PC2,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle("(b) Noise")+
  scale_colour_manual(values=cbPalette)+
  theme(legend.position="none")+
  scale_shape_manual(values=c(4,20))
pcaPlot12.noise

pcaPlot12.combat<-ggplot2::ggplot(data=pcaCombat,aes(x=PC1,y=PC2,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle("(c) Combat")+
  labs(shape="Data")+
  theme(legend.position="none")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(4,20))
pcaPlot12.combat

pcaPlot34.Original<-ggplot2::ggplot(data=pcaOriginal,aes(x=PC3,y=PC4,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle("(d) Original")+
  theme(legend.position="none")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(20,20))

pcaPlot34.Original

pcaPlot34.noise<-ggplot2::ggplot(data=pcaNoise,aes(x=PC3,y=PC4,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle("(e) Noise")+
  scale_colour_manual(values=cbPalette)+
  theme(legend.key.height = unit(0.9,"cm"))+
  theme(legend.position = "bottom")+
  scale_shape_manual(values=c(4,20))
pcaPlot34.noise

pcaPlot34.combat<-ggplot2::ggplot(data=pcaCombat,aes(x=PC3,y=PC4,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle("(f) Combat")+
  labs(shape="Data")+
  theme(legend.position="none")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(4,20))
pcaPlot34.combat

blank<-ggplot2::ggplot()+geom_blank()+theme_void()

pdf("Figure4.pdf",height=12,width=18)
grid.arrange(pcaPlot12.Original,pcaPlot12.noise,pcaPlot12.combat,
             pcaPlot34.Original,pcaPlot34.noise,pcaPlot34.combat,blank,blank,
             layout_matrix=cbind(c(1,1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,7),
                                 c(2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5),
                                 c(3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6,8)))
dev.off()

svg("Figure4.svg",height=12,width=18)
grid.arrange(pcaPlot12.Original,pcaPlot12.noise,pcaPlot12.combat,
             pcaPlot34.Original,pcaPlot34.noise,pcaPlot34.combat,blank,blank,
             layout_matrix=cbind(c(1,1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,7),
                                 c(2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5),
                                 c(3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6,8)))
dev.off()


##############
###Figure 5###
##############
fontsize=18
Fig5<-ggplot2::ggplot(ggplotDataCombinedFig3a,aes(x=tissues,y=correlation,fill=status))+
  geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge",aes(group=status))+
  theme_bw(fontsize)+
  xlab("")+
  ylab("Ontology score (Spearman)")+
  labs(fill=" ")+
  theme(legend.key.height = unit(1.5,"cm"))+
  scale_fill_manual(values=cbPalette)

pdf("Figure5.pdf",height=6,width=9)
Fig5
dev.off()

svg("Figure5.svg",height=6,width=9)
Fig5
dev.off()




####################################
###Supplementary Figs for Figure 4###
####################################
SupFig4_Bar<-ggplot2::ggplot(ggplotDataCombined.Spearman,aes(x=tissues,y=correlation,fill=status))+
  geom_boxplot()+
  facet_grid(ontology ~similarity.measure)+
  theme_bw(fontsize)+
  xlab("Tissues")+
  ylab("Ontology score")+
  labs(fill="Data")+
  theme(legend.key.height = unit(1.5,"cm"))+
  scale_fill_manual(values=cbPalette)

pdf("Sup_Fig_4_a.pdf",width=16,height=10)
SupFig4_Bar
dev.off()

svg("Sup_Fig_4_a.svg",width=16,height=10)
SupFig4_Bar
dev.off()


pcaPlot12.original<-ggplot2::ggplot(data=pcaOriginal,aes(x=PC1,y=PC2,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle("(a) Original")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(20,20))
pcaPlot12.original

pcaPlot23.original<-ggplot2::ggplot(data=pcaOriginal,aes(x=PC2,y=PC3,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle(" ")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(20,20))
pcaPlot23.original

pcaPlot34.original<-ggplot2::ggplot(data=pcaOriginal,aes(x=PC3,y=PC4,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle(" ")+
  theme(legend.position = "none")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(20,20))
pcaPlot34.original

pcaPlot12.noise<-ggplot2::ggplot(data=pcaNoise,aes(x=PC1,y=PC2,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle("(b) Noise ")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(4,20))
pcaPlot12.noise

pcaPlot23.noise<-ggplot2::ggplot(data=pcaNoise,aes(x=PC2,y=PC3,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle(" ")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(4,20))
pcaPlot23.noise

pcaPlot34.noise<-ggplot2::ggplot(data=pcaNoise,aes(x=PC3,y=PC4,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle(" ")+
  theme(legend.position = "right")+
  theme(legend.key.height = unit(1.8,"cm"))+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(4,20))
pcaPlot34.noise

pcaPlot12.combat<-ggplot2::ggplot(data=pcaCombat,aes(x=PC1,y=PC2,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle("(c) Combat")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(4,20))
pcaPlot12.combat

pcaPlot23.combat<-ggplot2::ggplot(data=pcaCombat,aes(x=PC2,y=PC3,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle(" ")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(4,20))
pcaPlot23.combat

pcaPlot34.combat<-ggplot2::ggplot(data=pcaCombat,aes(x=PC3,y=PC4,color=tissueType_2,shape=Batch))+
  geom_point(size=3)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle(" ")+
  theme(legend.position = "none")+
  scale_colour_manual(values=cbPalette)+
  scale_shape_manual(values=c(4,20))
pcaPlot34.combat


blank<-ggplot2::ggplot()+geom_blank()+theme_void()

pdf("Sup_Fig_4_b.pdf",width=24,height=24)
grid.arrange(pcaPlot12.original,pcaPlot23.original,pcaPlot34.original,
             pcaPlot12.noise,pcaPlot23.noise,pcaPlot34.noise,
             pcaPlot12.combat,pcaPlot23.combat,pcaPlot34.combat,blank,blank,
             layout_matrix=cbind(c(1,4,7),
                                 c(1,4,7),
                                 c(1,4,7),
                                 c(1,4,7),
                                 c(2,5,8),
                                 c(2,5,8),
                                 c(2,5,8),
                                 c(2,5,8),
                                 c(3,6,9),
                                 c(3,6,9),
                                 c(3,6,9),
                                 c(3,6,9),c(10,6,11)))
dev.off()


svg("Sup_Fig_4_b.svg",width=24,height=24)
grid.arrange(pcaPlot12.original,pcaPlot23.original,pcaPlot34.original,
             pcaPlot12.noise,pcaPlot23.noise,pcaPlot34.noise,
             pcaPlot12.combat,pcaPlot23.combat,pcaPlot34.combat,blank,blank,
             layout_matrix=cbind(c(1,4,7),
                                 c(1,4,7),
                                 c(1,4,7),
                                 c(1,4,7),
                                 c(2,5,8),
                                 c(2,5,8),
                                 c(2,5,8),
                                 c(2,5,8),
                                 c(3,6,9),
                                 c(3,6,9),
                                 c(3,6,9),
                                 c(3,6,9),c(10,6,11)))
dev.off()


svg("Sup_Fig_GTEX_Logarithm.svg",width=16,height=9)
ggplot2::ggplot(ggplotDataNoLog,aes(x=tissues,y=correlation,fill=status))+
  geom_boxplot()+
  facet_grid(ontology ~similarity.measure)+
  theme_bw(fontsize)+
  xlab("Tissues")+
  ylab("Ontology score")+
  labs(fill="Data")+
  theme(legend.key.height = unit(1.5,"cm"))+
  scale_fill_manual(values=cbPalette)
dev.off()


pdf("Sup_Fig_GTEX_Logarithm.pdf",width=16,height=9)
ggplot2::ggplot(ggplotDataNoLog,aes(x=tissues,y=correlation,fill=status))+
  geom_boxplot()+
  facet_grid(ontology ~similarity.measure)+
  theme_bw(fontsize)+
  xlab("Tissues")+
  ylab("Ontology score")+
  labs(fill="Data")+
  theme(legend.key.height = unit(1.5,"cm"))+
  scale_fill_manual(values=cbPalette)
dev.off()
