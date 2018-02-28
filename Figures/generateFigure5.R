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

gtexTCGA.Fig5<-ggplot2::ggplot(ggplotDataFig5.gtexTCGA,aes(x=Batch,y=correlation,fill=status))+
  geom_boxplot()+
  theme_bw(fontsize)+
  labs(fill="Data")+
  theme(legend.key.height = unit(1.5,"cm"))+
  xlab("")+
  ylab("Spearman correlation")+
  theme(legend.position = "none")+
  scale_fill_manual(values=cbPalette)+
  ggtitle("(a) GTEx & TCGA")
gtexTCGA.Fig5

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

########################
###IHEC use case###
########################
###############
###Load data###
###############
ihec.joinedMatrix.original<-read.table("Data/IHEC_PCA_Data_Original_Outlier_Removed.txt",header=T,sep="\t")
ihec.joinedMatrix.sva<-read.table("Data/IHEC_PCA_SVA.txt",header=T,sep="\t")
ihec.joinedMatrix.combat<-read.table("Data/IHEC_PCA_Data_Combat.txt",header=T,sep="\t")
ihec.joinedMatrix.RUV<-read.table("Data/IHEC_PCA_RUV.txt",header=T,sep="\t")
ihec.ontologyScores<-read.table("Data/IHEC_RNA_seq_ontology.txt",header=T,sep="\t")

###################################
###Generate ontology score plots###
###################################
fontsize=22
ggplotDataCombined.Spearman<-ihec.ontologyScores
ihec.ontologyScores<-ihec.ontologyScores[c(which(ihec.ontologyScores$dataset=="Original"),which(ihec.ontologyScores$dataset=="Combat"),which(ihec.ontologyScores$dataset=="SVA"),which(ihec.ontologyScores$dataset=="RUV")),]
ihec.ontologyScores$dataset<-factor(ihec.ontologyScores$dataset,levels=c("Original","Combat","SVA","RUV"))
ihec.ontologyScores<-ihec.ontologyScores[c(which(ihec.ontologyScores$similarity=="Pearson"),which(ihec.ontologyScores$similarity=="Spearman")),]
ggplotDataFig5.ihec<-ihec.ontologyScores[which(ihec.ontologyScores$similarity=="Spearman"),]
ggplotDataFig5.ihec<-ggplotDataFig5.ihec[which(ggplotDataFig5.ihec$ontology=="Cosine"),]

ihec.Fig5<-ggplot2::ggplot(ggplotDataFig5.ihec,aes(x=Consortia,y=score,fill=dataset))+
  geom_boxplot()+
  theme_bw(fontsize)+
  labs(fill="Data")+
  theme(legend.key.height = unit(1.5,"cm"))+
  xlab("")+
  ylab("Spearman correlation")+
  theme(legend.position = "right")+
  scale_fill_manual(values=cbPalette)+
  ggtitle("(b) IHEC")
ihec.Fig5

pdf("Figure5.pdf",width=15,height=7)
grid.arrange(gtexTCGA.Fig5,ihec.Fig5,layout_matrix=cbind(1,1,1,1,1,1,2,2,2,2,2,2,2))
dev.off()

svg("Figure5.svg",width=15,height=7)
grid.arrange(gtexTCGA.Fig5,ihec.Fig5,layout_matrix=cbind(1,1,1,1,1,1,2,2,2,2,2,2,2))
dev.off()



ihec.acrossTissues.Ontology<-ggplot2::ggplot(ihec.ontologyScores,aes(x=dataset,y=score,fill=dataset))+
  geom_boxplot()+
  theme_bw(fontsize)+
  facet_grid(ontology ~similarity)+
  labs(fill="Data")+
  theme(legend.key.height = unit(1.5,"cm"))+
  xlab(" ")+ylab("Correlation")+
  theme(legend.position="none")+
  scale_fill_manual(values=cbPalette)
ihec.acrossTissues.Ontology

pdf("Ontology_IHEC_Sup_Fig_1.pdf",width=12,height=8)
ihec.acrossTissues.Ontology
dev.off()


svg("Ontology_IHEC_Sup_Fig_1.svg",width=12,height=8)
ihec.acrossTissues.Ontology
dev.off()

ihec.tissue.specific.Ontology<-ggplot2::ggplot(ggplotDataFig5.ihec,aes(x=dataset,y=score,fill=dataset))+
  geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  theme_bw(fontsize)+
  facet_wrap(~CellType..Condensed.,ncol=6)+
  labs(fill="Data")+
  theme(legend.key.width = unit(1.5,"cm"))+
  xlab(" ")+ylab("Correlation")+
  theme(legend.position="bottom")+
  scale_fill_manual(values=cbPalette)
ihec.tissue.specific.Ontology

pdf("Ontology_IHEC_Sup_Fig_2_Spearman_Cosine.pdf",width=28,height=40)
ihec.tissue.specific.Ontology
dev.off()

svg("Ontology_IHEC_Sup_Fig_2_Spearman_Cosine.svg",width=28,height=40)
ihec.tissue.specific.Ontology
dev.off()

###############
###PCA plots###
###############
ihec.pcaPlot12.original<-ggplot2::ggplot(data=ihec.joinedMatrix.original,aes(x=PC1,y=PC2,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC1"))+ 
  ylab(paste0("PC2"))+
  ggtitle("(a) Original")
ihec.pcaPlot12.original

ihec.pcaPlot23.original<-ggplot2::ggplot(data=ihec.joinedMatrix.original,aes(x=PC2,y=PC3,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC2"))+ 
  ylab(paste0("PC3"))+
  ggtitle(" ")
ihec.pcaPlot23.original

ihec.pcaPlot34.original<-ggplot2::ggplot(data=ihec.joinedMatrix.original,aes(x=PC3,y=PC4,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC3"))+ 
  ylab(paste0("PC4"))+
  ggtitle(" ")
ihec.pcaPlot34.original

ihec.pcaPlot12.combat<-ggplot2::ggplot(data=ihec.joinedMatrix.combat,aes(x=PC1,y=PC2,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC1"))+ 
  ylab(paste0("PC2"))+
  ggtitle("(b) Combat")
ihec.pcaPlot12.combat

ihec.pcaPlot23.combat<-ggplot2::ggplot(data=ihec.joinedMatrix.combat,aes(x=PC2,y=PC3,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC2"))+ 
  ylab(paste0("PC3"))+
  ggtitle(" ")
ihec.pcaPlot23.combat

ihec.pcaPlot34.combat<-ggplot2::ggplot(data=ihec.joinedMatrix.combat,aes(x=PC3,y=PC4,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC3"))+ 
  ylab(paste0("PC4"))+
  ggtitle(" ")
ihec.pcaPlot34.combat

ihec.pcaPlot12.sva<-ggplot2::ggplot(data=ihec.joinedMatrix.sva,aes(x=PC1,y=PC2,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC1"))+ 
  ylab(paste0("PC2"))+
  ggtitle("(c) SVA")
ihec.pcaPlot12.sva

ihec.pcaPlot23.sva<-ggplot2::ggplot(data=ihec.joinedMatrix.sva,aes(x=PC2,y=PC3,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC2"))+ 
  ylab(paste0("PC3"))+
  ggtitle(" ")
ihec.pcaPlot23.sva

ihec.pcaPlot34.sva<-ggplot2::ggplot(data=ihec.joinedMatrix.sva,aes(x=PC3,y=PC4,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC3"))+ 
  ylab(paste0("PC4"))+
  ggtitle(" ")
ihec.pcaPlot34.sva

ihec.pcaPlot12.RUV<-ggplot2::ggplot(data=ihec.joinedMatrix.RUV,aes(x=PC1,y=PC2,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC1"))+ 
  ylab(paste0("PC2"))+
  ggtitle("(d) RUV")
ihec.pcaPlot12.RUV

ihec.pcaPlot23.RUV<-ggplot2::ggplot(data=ihec.joinedMatrix.RUV,aes(x=PC2,y=PC3,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "bottom")+
  xlab(paste0("PC2"))+ 
  ylab(paste0("PC3"))+
  ggtitle(" ")+
  labs(shape="Consortia")+
  scale_colour_discrete(guide = FALSE)

ihec.pcaPlot34.RUV<-ggplot2::ggplot(data=ihec.joinedMatrix.RUV,aes(x=PC3,y=PC4,color=Celltype..Original.,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC3"))+ 
  ylab(paste0("PC4"))+
  ggtitle(" ")
ihec.pcaPlot34.RUV

pdf("PCA_iheca_Sup_Fig.pdf",width=20,height=24)
pca_ihec<-grid.arrange(ihec.pcaPlot12.original,ihec.pcaPlot23.original,ihec.pcaPlot34.original,
             ihec.pcaPlot12.combat,ihec.pcaPlot23.combat,ihec.pcaPlot34.combat,
             ihec.pcaPlot12.sva,ihec.pcaPlot23.sva,ihec.pcaPlot34.sva,
             ihec.pcaPlot12.RUV,ihec.pcaPlot23.RUV,ihec.pcaPlot34.RUV,blank,blank,
             layout_matrix=cbind(c(1,1,1,1,1,1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,4,4,4,4,7,7,7,7,7,7,7,7,7,7,7,7,10,10,10,10,10,10,10,10,10,10,10,10,13),
                                 c(2,2,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5,5,5,8,8,8,8,8,8,8,8,8,8,8,8,11,11,11,11,11,11,11,11,11,11,11,11,11),
                                 c(3,3,3,3,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6,6,6,6,6,9,9,9,9,9,9,9,9,9,9,9,9,12,12,12,12,12,12,12,12,12,12,12,12,14)))
dev.off()


svg("PCA_iheca_Sup_Fig.svg",width=20,height=24)
pca_ihec<-grid.arrange(ihec.pcaPlot12.original,ihec.pcaPlot23.original,ihec.pcaPlot34.original,
                       ihec.pcaPlot12.combat,ihec.pcaPlot23.combat,ihec.pcaPlot34.combat,
                       ihec.pcaPlot12.sva,ihec.pcaPlot23.sva,ihec.pcaPlot34.sva,
                       ihec.pcaPlot12.RUV,ihec.pcaPlot23.RUV,ihec.pcaPlot34.RUV,blank,blank,
                       layout_matrix=cbind(c(1,1,1,1,1,1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,4,4,4,4,7,7,7,7,7,7,7,7,7,7,7,7,10,10,10,10,10,10,10,10,10,10,10,10,13),
                                           c(2,2,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5,5,5,8,8,8,8,8,8,8,8,8,8,8,8,11,11,11,11,11,11,11,11,11,11,11,11,11),
                                           c(3,3,3,3,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6,6,6,6,6,9,9,9,9,9,9,9,9,9,9,9,9,12,12,12,12,12,12,12,12,12,12,12,12,14)))
dev.off()



###########################
###Correcting hepatocyte###
###########################

ihec.joinedMatrix.original.hep<-ihec.joinedMatrix.original
ihec.joinedMatrix.original.hep$CellType..Condensed<-as.character(ihec.joinedMatrix.original.hep$CellType..Condensed)
ihec.joinedMatrix.original.hep$CellType..Condensed[intersect(which(ihec.joinedMatrix.original$CellType..Condensed!="Kidney"),
                                                             which(ihec.joinedMatrix.original$CellType..Condensed!="Hepatocyte"))]<-"Other"
ihec.joinedMatrix.original.hep$CellType..Condensed<-factor(ihec.joinedMatrix.original.hep$CellType..Condensed)
ihec.pcaPlot12.Original.hep<-ggplot2::ggplot(data=ihec.joinedMatrix.original.hep,aes(x=PC1,y=PC2,color=CellType..Condensed,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC1"))+ 
  ylab(paste0("PC2"))+
  ggtitle("(a) Original")+
  labs(shape="Consortia")+
  scale_color_manual(values=c("green1","black","grey"))
ihec.pcaPlot12.Original.hep

ihec.pcaPlot23.Original.hep<-ggplot2::ggplot(data=ihec.joinedMatrix.original.hep,aes(x=PC2,y=PC3,color=CellType..Condensed,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC2"))+ 
  ylab(paste0("PC3"))+
  ggtitle(" ")+
  labs(shape="Consortia")+
  scale_color_manual(values=c("green1","black","grey"))
ihec.pcaPlot23.Original.hep

ihec.pcaPlot34.Original.hep<-ggplot2::ggplot(data=ihec.joinedMatrix.original.hep,aes(x=PC3,y=PC4,color=CellType..Condensed,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC3"))+ 
  ylab(paste0("PC4"))+
  ggtitle(" ")+
  labs(shape="Consortia")+
  scale_color_manual(values=c("green1","black","grey"))
ihec.pcaPlot34.Original.hep



ihec.ruv.perform.hep<-ihec.joinedMatrix.RUV
ihec.ruv.perform.hep$CellType..Condensed<-as.character(ihec.ruv.perform.hep$CellType..Condensed)
ihec.ruv.perform.hep$CellType..Condensed[intersect(which(ihec.joinedMatrix.RUV$CellType..Condensed!="Kidney"),
                                                   which(ihec.joinedMatrix.RUV$CellType..Condensed!="Hepatocyte"))]<-"Other"
ihec.ruv.perform.hep$CellType..Condensed<-factor(ihec.ruv.perform.hep$CellType..Condensed)
ihec.pcaPlot12.RUV.hep<-ggplot2::ggplot(data=ihec.ruv.perform.hep,aes(x=PC1,y=PC2,color=CellType..Condensed,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC1"))+ 
  ylab(paste0("PC2"))+
  ggtitle("(b) RUV ")+
  labs(shape="Consortia")+
  scale_color_manual(values=c("green1","black","grey"))
ihec.pcaPlot12.RUV.hep

ihec.pcaPlot23.RUV.hep<-ggplot2::ggplot(data=ihec.ruv.perform.hep,aes(x=PC2,y=PC3,color=CellType..Condensed,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "bottom")+
  xlab(paste0("PC2"))+ 
  ylab(paste0("PC3"))+
  ggtitle(" ")+
  labs(shape="Consortia")+
  labs(color="Celltype")+
  scale_color_manual(values=c("green1","black","grey"))
ihec.pcaPlot23.RUV.hep

ihec.pcaPlot34.RUV.hep<-ggplot2::ggplot(data=ihec.ruv.perform.hep,aes(x=PC3,y=PC4,color=CellType..Condensed,shape=Consortia))+
  geom_point(size=pointsize)+
  theme_bw(fontsize)+
  theme(legend.position = "none")+
  xlab(paste0("PC3"))+ 
  ylab(paste0("PC4"))+
  ggtitle(" ")+
  labs(shape="Consortia")+
  scale_color_manual(values=c("green1","black","grey"))
ihec.pcaPlot34.RUV.hep

svg("Sup_Fig_Hep_Kid.svg",width=18,height=12)
pca_ihec<-grid.arrange(ihec.pcaPlot12.Original.hep,ihec.pcaPlot23.Original.hep,ihec.pcaPlot34.Original.hep,
                       ihec.pcaPlot12.RUV.hep,ihec.pcaPlot23.RUV.hep,ihec.pcaPlot34.RUV.hep,blank,blank,
                       layout_matrix=cbind(c(1,1,1,1,1,1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,4,4,4,4,7),
                                           c(2,2,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5,5,5,5),
                                           c(3,3,3,3,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6,6,6,6,6,8)))
dev.off()


pdf("Sup_Fig_Hep_Kid.pdf",width=18,height=12)
pca_ihec<-grid.arrange(ihec.pcaPlot12.Original.hep,ihec.pcaPlot23.Original.hep,ihec.pcaPlot34.Original.hep,
                       ihec.pcaPlot12.RUV.hep,ihec.pcaPlot23.RUV.hep,ihec.pcaPlot34.RUV.hep,blank,blank,
                       layout_matrix=cbind(c(1,1,1,1,1,1,1,1,1,1,1,1,4,4,4,4,4,4,4,4,4,4,4,4,7),
                                           c(2,2,2,2,2,2,2,2,2,2,2,2,5,5,5,5,5,5,5,5,5,5,5,5,5),
                                           c(3,3,3,3,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6,6,6,6,6,8)))
dev.off()

