ggplot2Package<-require("ggplot2")
gridPackage<-require("gridExtra")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73")

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
ggplotDataFig7.ihec<-ihec.ontologyScores[which(ihec.ontologyScores$similarity=="Spearman"),]
ggplotDataFig7.ihec<-ggplotDataFig7.ihec[which(ggplotDataFig7.ihec$ontology=="Cosine"),]
ggplotDataFig7.ihec<-ggplotDataFig7.ihec[which(ggplotDataFig7.ihec$pca=="PC 1-4"),]

ihec.Fig7<-ggplot2::ggplot(ggplotDataFig7.ihec,aes(x=Consortia,y=score,fill=dataset))+
  geom_boxplot()+
  theme_bw(fontsize)+
  labs(fill=" ")+
  theme(legend.key.height = unit(1.5,"cm"))+
  xlab("")+
  ylab("Ontology score (Spearman)")+
  theme(legend.position = "none")+
  scale_fill_manual(values=cbPalette)+
  ggtitle("(a)")
ihec.Fig7

sampleTissuesSmallExample<-ggplotDataFig7.ihec[c(which(ggplotDataFig7.ihec$Celltype..Original.=="Hepatocyte"),which(ggplotDataFig7.ihec$Celltype..Original.=="Erythroblast")),]
ihec.tissue.specific.Ontology.Small<-ggplot2::ggplot(sampleTissuesSmallExample,aes(x=dataset,y=score,fill=dataset))+
  geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  theme_bw(fontsize)+
  facet_grid(CellType..Condensed.~.)+
  labs(fill=" ")+
  theme(legend.key.height = unit(1.5,"cm"))+
  xlab(" ")+ylab("Ontology score (Spearman)")+
  theme(legend.position="right")+
  scale_fill_manual(values=cbPalette)+
  ggtitle("(b)")
ihec.tissue.specific.Ontology.Small

pdf("Figure7.pdf",width=16,height=7)
grid.arrange(ihec.Fig7,ihec.tissue.specific.Ontology.Small,layout_matrix=cbind(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2))
dev.off()

svg("Figure7.svg",width=16,height=7)
grid.arrange(ihec.Fig7,ihec.tissue.specific.Ontology.Small,layout_matrix=cbind(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2))
dev.off()



ihec.acrossTissues.Ontology<-ggplot2::ggplot(ihec.ontologyScores,aes(x=dataset,y=score,fill=dataset))+
  geom_boxplot()+
  theme_bw(fontsize)+
  facet_grid(ontology ~similarity)+
  labs(fill="Data")+
  theme(legend.key.height = unit(1.5,"cm"))+
  xlab(" ")+ylab("Ontology score")+
  theme(legend.position="none")+
  scale_fill_manual(values=cbPalette)
ihec.acrossTissues.Ontology

pdf("Ontology_IHEC_Sup_Fig_1.pdf",width=12,height=8)
ihec.acrossTissues.Ontology
dev.off()


svg("Ontology_IHEC_Sup_Fig_1.svg",width=12,height=8)
ihec.acrossTissues.Ontology
dev.off()

ihec.tissue.specific.Ontology<-ggplot2::ggplot(ggplotDataFig7.ihec,aes(x=dataset,y=score,fill=dataset))+
  geom_bar(stat="summary",position="dodge",fun.y = "mean")+
  stat_summary(fun.data = mean_se, geom = "errorbar",position="dodge")+
  theme_bw(fontsize)+
  facet_wrap(~CellType..Condensed.,ncol=6)+
  labs(fill="Data")+
  theme(legend.key.width = unit(1.5,"cm"))+
  xlab(" ")+ylab("Ontology score")+
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

