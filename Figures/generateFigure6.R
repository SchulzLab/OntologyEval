ggplot2Package<-require("ggplot2")
gridPackage<-require("gridExtra")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
########################
###gtex TCGA use case###
########################
###############
###Load data###
###############
gtex_tcga.joinedMatrix.original<-read.table("Data/originalPCA.txt",header=T)
gtex_tcga.joinedMatrix.sva<-read.table("Data/svaPCA.txt",header=T)
gtex_tcga.joinedMatrix.combat<-read.table("Data/combatPCA.txt",header=T)
gtex_tcga.joinedMatrix.RUV<-read.table("Data/ruvPCA.txt",header=T)
gtex_tcga.ontologyScores<-read.table("Data/ggplotDataCombinedBatch.txt",header=T,row.names=1)

###################################
###Generate ontology score plots###
###################################
fontsize=22
ggplotDataCombined.Spearman<-gtex_tcga.ontologyScores[which(gtex_tcga.ontologyScores$expr.distance=="Spearman"),]
ggplotDataCombined.Spearman$status<-factor(ggplotDataCombined.Spearman$status,levels=c("Original","Combat","SVA","RUV"))
ggplotDataFig5.gtexTCGA<-ggplotDataCombined.Spearman[which(ggplotDataCombined.Spearman$similarity.measure=="Spearman"),]
ggplotDataFig5.gtexTCGA<-ggplotDataFig5.gtexTCGA[which(ggplotDataFig5.gtexTCGA$ontology=="Cosine"),]
 gtexTCGA.acrossTissues.Ontology<-ggplot2::ggplot(ggplotDataCombined.Spearman,aes(x=status,y=correlation,fill=status))+
  geom_boxplot()+
  theme_bw(fontsize)+
  facet_grid(ontology ~similarity.measure)+
  labs(fill="Data")+
  theme(legend.key.height = unit(1.5,"cm"))+
  xlab(" ")+ylab("Correlation")+
  ggtitle("(a)")+
  theme(legend.position="none")+
  scale_fill_manual(values=cbPalette)
gtexTCGA.acrossTissues.Ontology

gtexTCGA.tissue.specific.Ontology<-ggplot2::ggplot(ggplotDataCombined.Spearman,aes(x=tissues,y=correlation,fill=status))+
  geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  theme_bw(fontsize)+
  facet_grid(ontology ~similarity.measure)+
  labs(fill="Data")+
  theme(legend.key.width = unit(1.5,"cm"))+
  xlab(" ")+ylab("Correlation")+
  ggtitle("(b)")+
  theme(legend.position="bottom")+
  scale_fill_manual(values=cbPalette)+
  theme(legend.position="none")
gtexTCGA.tissue.specific.Ontology

ggplotDataCombined.Spearman.Spearman<-ggplotDataCombined.Spearman[which(ggplotDataCombined.Spearman$similarity.measure=="Spearman"),]
gtexTCGA.tissue.consortia.specific.Ontology<-ggplot2::ggplot(ggplotDataCombined.Spearman.Spearman,aes(x=tissues,y=correlation,fill=status))+
  geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  theme_bw(fontsize)+
  facet_grid(ontology ~ Batch)+
  labs(fill="Data")+
  theme(legend.key.width = unit(1.5,"cm"))+
  xlab(" ")+ylab("Correlation")+
  ggtitle("(c)")+
  theme(legend.position="bottom")+
  scale_fill_manual(values=cbPalette)
gtexTCGA.tissue.consortia.specific.Ontology

pdf("Ontology_gtex_tcga_Sup_Fig.pdf",width=12,height=21)
grid.arrange(gtexTCGA.acrossTissues.Ontology,gtexTCGA.tissue.specific.Ontology,gtexTCGA.tissue.consortia.specific.Ontology,
             layout_matrix=cbind(c(1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3)))
dev.off()

svg("Ontology_gtex_tcga_Sup_Fig.svg",width=12,height=21)
grid.arrange(gtexTCGA.acrossTissues.Ontology,gtexTCGA.tissue.specific.Ontology,gtexTCGA.tissue.consortia.specific.Ontology,
             layout_matrix=cbind(c(1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3)))
dev.off()

###################
###Generate PCAs###
###################

cbPalette2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
fontsize=22
pointsize=3
gtex_tcga.pcaPlot12.original<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.original,aes(x=PC1,y=PC2,color=tissueType_2,shape=studyName))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC1"))+ #,round((pca$sdev/sum(pca$sdev))[1]*100,2),"% variance explained"))+
  ylab(paste0("PC2"))+ #,round((pca$sdev/sum(pca$sdev))[2]*100,2),"% variance explained"))+
  ggtitle("(a) Original")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtex_tcga.pcaPlot12.original

gtex_tcga.pcaPlot23.original<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.original,aes(x=PC2,y=PC3,color=tissueType_2,shape=studyName))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC2"))+#,round((pca$sdev/sum(pca$sdev))[2]*100,2),"% variance explained"))+
  ylab(paste0("PC3"))+#,round((pca$sdev/sum(pca$sdev))[3]*100,2),"% variance explained"))+
  ggtitle(" ")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtex_tcga.pcaPlot23.original

gtex_tcga.pcaPlot34.original<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.original,aes(x=PC3,y=PC4,color=tissueType_2,shape=studyName))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  xlab(paste0("PC3 "))+#,round((pca$sdev/sum(pca$sdev))[3]*100,2),"% variance explained"))+
  ylab(paste0("PC4 "))+#,round((pca$sdev/sum(pca$sdev))[4]*100,2),"% variance explained"))+
  labs(color="Tissue")+
  ggtitle(" ")+
  theme(legend.position = "none")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))#+
#
gtex_tcga.pcaPlot34.original

gtex_tcga.pcaPlot12.combat<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.combat,aes(x=PC1,y=PC2,color=Tissue,shape=Batch))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle("(b) Combat")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtex_tcga.pcaPlot12.combat

gtex_tcga.pcaPlot23.combat<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.combat,aes(x=PC2,y=PC3,color=Tissue,shape=Batch))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle(" ")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtex_tcga.pcaPlot23.combat

gtex_tcga.pcaPlot34.combat<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.combat,aes(x=PC3,y=PC4,color=Tissue,shape=Batch))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle(" ")+
  theme(legend.position = "none")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtex_tcga.pcaPlot34.combat

gtex_tcga.pcaPlot12.sva<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.sva,aes(x=PC1,y=PC2,color=tissueType_2,shape=studyName))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle("(c) SVA")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtex_tcga.pcaPlot12.sva

gtex_tcga.pcaPlot23.sva<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.sva,aes(x=PC2,y=PC3,color=tissueType_2,shape=studyName))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle(" ")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtex_tcga.pcaPlot23.sva

gtex_tcga.pcaPlot34.sva<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.sva,aes(x=PC3,y=PC4,color=tissueType_2,shape=studyName))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle(" ")+
  theme(legend.position = "none")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtex_tcga.pcaPlot34.sva

gtex_tcga.pcaPlot12.RUV<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.RUV,aes(x=PC1,y=PC2,color=tissueType_2,shape=studyName))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle("(d) RUV")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtex_tcga.pcaPlot12.RUV

gtex_tcga.pcaPlot23.RUV<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.RUV,aes(x=PC2,y=PC3,color=tissueType_2,shape=studyName))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "bottom")+
  ggtitle(" ")+
  labs(shape="Data source")+
  labs(color="Tissue")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtex_tcga.pcaPlot23.RUV

gtex_tcga.pcaPlot34.RUV<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.RUV,aes(x=PC3,y=PC4,color=tissueType_2,shape=studyName))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  ggtitle(" ")+
  theme(legend.position = "none")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtex_tcga.pcaPlot34.RUV

blank<-ggplot2::ggplot()+geom_blank()+theme_void()

pdf("PCA_gtex_tcga_Sup_Fig.pdf",width=18,height=24)
grid.arrange(gtex_tcga.pcaPlot12.original,gtex_tcga.pcaPlot23.original,gtex_tcga.pcaPlot34.original,
             gtex_tcga.pcaPlot12.combat,gtex_tcga.pcaPlot23.combat,gtex_tcga.pcaPlot34.combat,
             gtex_tcga.pcaPlot12.sva,gtex_tcga.pcaPlot23.sva,gtex_tcga.pcaPlot34.sva,
             gtex_tcga.pcaPlot12.RUV,gtex_tcga.pcaPlot23.RUV,gtex_tcga.pcaPlot34.RUV,blank,blank,
             layout_matrix=cbind(c(1,1,1,1,1,1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,4,4,4,4,7,7,7,7,7,7,7,7,7,7,7,7,10,10,10,10,10,10,10,10,10,10,10,10,13),
                                 c(2,2,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5,5,5,8,8,8,8,8,8,8,8,8,8,8,8,11,11,11,11,11,11,11,11,11,11,11,11,11),
                                 c(3,3,3,3,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6,6,6,6,6,9,9,9,9,9,9,9,9,9,9,9,9,12,12,12,12,12,12,12,12,12,12,12,12,14)))
dev.off()

svg("PCA_gtex_tcga_Sup_Fig.svg",width=18,height=24)
grid.arrange(gtex_tcga.pcaPlot12.original,gtex_tcga.pcaPlot23.original,gtex_tcga.pcaPlot34.original,
             gtex_tcga.pcaPlot12.combat,gtex_tcga.pcaPlot23.combat,gtex_tcga.pcaPlot34.combat,
             gtex_tcga.pcaPlot12.sva,gtex_tcga.pcaPlot23.sva,gtex_tcga.pcaPlot34.sva,
             gtex_tcga.pcaPlot12.RUV,gtex_tcga.pcaPlot23.RUV,gtex_tcga.pcaPlot34.RUV,blank,blank,
             layout_matrix=cbind(c(1,1,1,1,1,1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,4,4,4,4,7,7,7,7,7,7,7,7,7,7,7,7,10,10,10,10,10,10,10,10,10,10,10,10,13),
                                 c(2,2,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5,5,5,8,8,8,8,8,8,8,8,8,8,8,8,11,11,11,11,11,11,11,11,11,11,11,11,11),
                                 c(3,3,3,3,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6,6,6,6,6,9,9,9,9,9,9,9,9,9,9,9,9,12,12,12,12,12,12,12,12,12,12,12,12,14)))
dev.off()

##############
###Figure 6###
##############
fontsize=18
gtexTCGA.Fig6a<-ggplot2::ggplot(ggplotDataFig5.gtexTCGA,aes(x=Batch,y=correlation,fill=status))+
  geom_boxplot()+
  theme_bw(fontsize)+
  labs(fill="Data")+
  theme(legend.key.height = unit(1.5,"cm"))+
  xlab("")+
  ylab("Spearman correlation")+
  theme(legend.position = "none")+
  scale_fill_manual(values=cbPalette)+
  ggtitle("(a)")
gtexTCGA.Fig6a

ggplotDataCombined.Spearman.Spearman.Cosine<-ggplotDataCombined.Spearman.Spearman[which(ggplotDataCombined.Spearman.Spearman$ontology=="Cosine"),]
gtexTCGA.Fig6b<-ggplot2::ggplot(ggplotDataCombined.Spearman.Spearman.Cosine,aes(x=tissues,y=correlation,fill=status))+
  geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  theme_bw(fontsize)+
  labs(fill=" ")+
  theme(legend.key.height = unit(0.75,"cm"))+
  xlab(" ")+ylab("Spearman correlation")+
  ggtitle("(b)")+
  theme(legend.position="right")+
  scale_fill_manual(values=cbPalette)+
  facet_wrap(~Batch)
gtexTCGA.Fig6b

gtexTCGA.Fig6c1<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.sva,aes(x=PC1,y=PC2,color=tissueType_2,shape=studyName))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  ggtitle("(c)")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))
gtexTCGA.Fig6c1


gtexTCGA.Fig6c2<-ggplot2::ggplot(data=gtex_tcga.joinedMatrix.sva,aes(x=PC3,y=PC4,color=tissueType_2,shape=studyName))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  labs(color="Tissue")+
  labs(shape="Consortium")+
  ggtitle(" ")+
  theme(legend.position = "right")+
  scale_colour_manual(values=cbPalette2)+
  scale_shape_manual(values=c(4,20))+
  theme(legend.key.height = unit(.75,"cm"))
gtexTCGA.Fig6c2


pdf("Figure6.pdf",width=17,height=10)
grid.arrange(gtexTCGA.Fig6a,gtexTCGA.Fig6b,gtexTCGA.Fig6c1,gtexTCGA.Fig6c2,layout_matrix=cbind(c(1,1,3,3,3),
                                                                                               c(1,1,3,3,3),
                                                                                               c(1,1,3,3,3),
                                                                                               c(1,1,3,3,3),
                                                                                               c(1,1,3,3,3),
                                                                                               c(2,2,3,3,3),
                                                                                               c(2,2,3,3,3),
                                                                                               c(2,2,4,4,4),
                                                                                               c(2,2,4,4,4),
                                                                                               c(2,2,4,4,4),
                                                                                               c(2,2,4,4,4),
                                                                                               c(2,2,4,4,4),
                                                                                               c(2,2,4,4,4),
                                                                                               c(2,2,4,4,4),
                                                                                               c(2,2,4,4,4)))
dev.off()
