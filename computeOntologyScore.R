library("tools")
ggplotAvailable<-require("ggplot2")

args<-commandArgs(TRUE)

###########################
###Function declarations###
###########################
############################
###Loading ontology table###
############################
load.Ontology<-function(filename=NULL){
  if (is.null(filename)){
    filename<-"Data/cosinesim.tsv"
  }
  ontology.scores<-read.table(filename)
  if (dim(ontology.scores)[2]!=3){
    stop("Ontology scores table must contain exactly three columns without a header:  <T1> tab <T2> tab <Cell Ontology Similarity>")
  }
  colnames(ontology.scores)<-c("T1","T2","Similarity")
  return(ontology.scores)
}

#############################
###Load ontology index map###
#############################
load.Idx2Term<-function(filename=NULL){
  if (is.null(filename)){
    filename<-"Data/idx2termid.tsv"
  }
  ontology.idx2termid<-read.table(filename)
  if (dim(ontology.idx2termid)[1]==0){
    stop("Index to Cell Ontology term mapping must provid at least one entry")
  }
  if (dim(ontology.idx2termid)[2]!=2){
    stop("Index to Cell Ontology Term mapping must contain only two columns without a header:  <Index> tab <Cell Ontology Term>")
  }
  colnames(ontology.idx2termid)<-c("Index","CLTermID")
  indices<-ontology.idx2termid$Index
  names(indices)<-ontology.idx2termid$CLTermID
  return(indices)
}

####################################
###Loading Sample -> Term Mapping###
####################################
load.Sample2Term<-function(filename=NULL){
  if (is.null(filename)){
    filename="Data/Example_Terms.txt"
  }
  ontology.sample2termid<-read.table(filename)
  if (dim(ontology.sample2termid)[1]==0){
    stop("Sample to Cell Ontology Term mapping must contain at least one entry.")
  }
  if (dim(ontology.sample2termid)[2]!=2){
    stop("Sample to Cell Ontology Term mapping must contain only two columns without a header:  <Sample ID> tab <Cell Ontology Term>")
  }
  colnames(ontology.sample2termid)<-c("SampleID","CLTermID")
  ontology.sample2termid$CLTermID<-as.character(ontology.sample2termid$CLTermID)
  return(ontology.sample2termid)
}

##########################
###Loading score matrix###
##########################
load.ObservedScoreMatrix<-function(filename=NULL){
  if (tools::file_ext(filename)=="rds"){
    ObservedScoreMatrix<-readRDS(filename)
  }else{
    ObservedScoreMatrix<-read.table(filename,header=T,row.names=1, stringsAsFactors = F)
  }
  if (dim(ObservedScoreMatrix)[2]==0){
    stop("No samples contained in data")
  }
  if (dim(ObservedScoreMatrix)[1]==0){
    stop("No genes contained in data")
  }
  return(ObservedScoreMatrix)
}

##########################################
###Computing distances on observed data###
##########################################
compute.Observed.Distance<-function(observed.data,Log2=T,Center=T,Scale=T,nPCA=4,method="spearman"){
  if (Log2){
    observed.data.pca<-prcomp(log2(observed.data+1),center=Center,scale. = Scale,retx= T)
  }else{
    observed.data.pca<-prcomp(observed.data,center=Center,scale. = Scale,retx= T)
  }
  #Compute correlation on PCA space
  coordinates<-cor(t(observed.data.pca$rotation[,1:nPCA]), use = "complete.obs", method = method)
  #Return a distance matrix (1.0-similarity)
  return(1.0-coordinates)
}

#########################################
###Generating ontology distance matrix###
#########################################
compute.Ontology.Distance<-function(similarities, idxMap, sampleMap){
  IndexMax<-length(unique(similarities$T1))-1
  ontology.matrix<-matrix(0,dim(sampleMap)[1],dim(sampleMap)[1])
  colnames(ontology.matrix)<-sampleMap$SampleID
  row.names(ontology.matrix)<-sampleMap$SampleID
  for (i in c(1:dim(sampleMap)[1])){
    #Retrieve CL term for sample i and map that to its index. Throws an error if the index can't be found in the term to index map
    firstIndex=tryCatch(idxMap[[sampleMap$CLTermID[i]]],error=function(cond){
                            print(paste0("Cell ontology Term",sampleMap$CLTermID[i]," could not be mapped. Original error:"));
                            print(cond);
                            return(NA)})
    if (is.na(firstIndex)){
      stop("Index could not be matched")
    }
    for (j in c(i:dim(sampleMap)[1])){
      #Retrieve CL term for sample j and map that to its index. Throws an error if the index can't be found in the term to index map
      secondIndex=tryCatch(idxMap[[sampleMap$CLTermID[j]]],error=function(cond){
        print(paste0("Cell ontology Term",sampleMap$CLTermID[j]," could not be mapped. Original error:"));
        print(cond);
        return(NA)})
      if (is.na(secondIndex)){
        stop("Index could not be matched")
      }
      #Due to the design of the ontology file, we need to make sure that firstIndex <= secondIndex
      if (firstIndex > secondIndex){
        temp=secondIndex
        secondIndex=firstIndex
        firstIndexT=temp
      }
      else{
        firstIndexT=firstIndex
      }
      #Retrieve the value from the similiarity matrix and convert to a distance matrix. The index computation is simply based on the gaussian sum formula
      value=1.0-similarities$Similarity[1+(firstIndexT*IndexMax)+secondIndex-((firstIndexT*(firstIndexT-1))/2)]
      ontology.matrix[i,j]=value
      ontology.matrix[j,i]=value
    }
  }
  return(ontology.matrix)
}  

###############################
###Generating ontology score###
###############################
compute.Ontology.Score<-function(expected.distances,observed.distances,sample2Term,method="spearman"){
  fixedNameOrdering<-colnames(expected.distances)
  reordered.observed.distances<-observed.distances[,match(fixedNameOrdering,colnames(observed.distances))]
  reordered.observed.distances<-reordered.observed.distances[match(fixedNameOrdering,row.names(reordered.observed.distances)),]
  resultVec<-c()
  #Pairwise score computation
  for (i in c(1:dim(reordered.observed.distances)[1])){
    resultVec<-rbind(resultVec,
                cbind(cor(reordered.observed.distances[,i], expected.distances[,i], use = "complete.obs", method = method),sample2Term[i,2]))
  }
  resultVec<-cbind(fixedNameOrdering,resultVec)
  colnames(resultVec)<-c("Sample","Score","TermID")
  resultVec<-data.frame(resultVec)
  resultVec[,2]<-as.numeric(as.character(resultVec[,2]))
  resultVec[,3]<-as.factor(resultVec[,3])
  return(resultVec)
}
 
################# 
###save output###
#################
save.ontology.score<-function(scoreMatrix,filename="Ontology_Score.txt"){
  if (tools::file_ext(filename)=="rds"){
    ObservedScoreMatrix<-saveRDS(scoreMatrix,file=filename)
  }else{
    ObservedScoreMatrix<-write.table(ontology.score,file=filename,sep="\t",quote=F,row.names=F,col.names = T)
  }
}

########################
###Generate a boxplot###
########################
generateBoxPlot<-function(ontology.score,fontsize=20,pngFile=NULL){
  library(ggplot2)
  p<-ggplot2::ggplot(ontology.score,aes(x=TermID,y=Score,fill=TermID))+
    geom_boxplot()+
    theme_bw(fontsize)+
    xlab("Cell Ontology Term")+
    ylab("Ontology Score")+
    theme(legend.key.height = unit(1.2,"cm"))+
    labs(fill=" ")
  if (!is.null(pngFile)){
  ggsave(pngFile,plot=p,width=3*length(unique(ontology.score$TermID)),height=7)
  }
  return(p)
}

##########
###main###
##########
###############################
###Parsing command arguments###
###############################
if(length(args) < 1) {
	cat("No arguments provided, running example. To show help, use the --help option")
}

if("--help" %in% args) {
          cat("
          OntologyEval computes an ontology score assessing sample similarity compared to an ontology based reference.\n
          Arguments:
          --Ontology: The ontology file to be used. Default is Data/cosinesim.tsv. \n
	--Idx2Term: The index to ontology term mapping file. Default is Data/idx2termid.tsv, \n
	--Sample2Term: The sampleID to ontology term mapping file. Default is the example file Data/Example_Terms.txt, \n
	--ObservedScoreMatrx: The matrix holding observed data, either in rds or txt format. Default is the example file Data/ExampleData.rds. \n
	--ObservedSimMethod: The method to assess similarity across the PCs on the observed data. Can be any of pearson, spearman, kendall. Default is spearman. \n
	--OntologySimMethod: The method to assess similarity across the distance vectors. Can be any of pearson, spearman, kendall. Default is spearman. \n
	--Output: Name of the output file holding the ontology scores. Default is Ontology_Score_Output.txt,\n
	--fontsize: fontsize to be used in a boxplot if generated. Default is 20. \n
	--Log2: TRUE (default) if observed data should be logarithmized, FALSE otherwise.\n
	--Center: TRUE (default) if observed data should be centered at 0, FALSE otherwise. \n
	--Scale: TRUE (default) if observed data should be scaled between 0 and 1, FALSE otherwise. \n
	--nPCA: Number of PC components to be considered to compute the distance on observed data. Default is 4.\n
	")
	stop("")
}

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsList <- as.list(as.character(argsDF$V2))
try(names(argsList) <- argsDF$V1)
cat("Parsing parameters")
if(is.null(argsList$Ontology)) {
  cat("No ontology file specified. Using default ontology file cosinesim.tsv\n")
  argsList$Ontology<-filename<-"Data/cosinesim.tsv"
}

if(is.null(argsList$Idx2Term)) {
  cat("No index to term mapping is specified. Using default mapping file idx2termid.tsv.\nThe file is supposed to be a 2 column file with row indices in the first, and Ontology IDs in the second column.\n")
  argsList$Idx2Term<-"Data/idx2termid.tsv"
}

if(is.null(argsList$Sample2Term)) {
  cat("No sample to term mapping is specified. Using example mapping file Example_Terms.txt.\nThe file is supposed to be a 2 column file with SampleIDs in the first, and Ontology IDs in the second column.\n")
  argsList$Sample2Term<-"Data/Example_Terms.txt"
}

if(is.null(argsList$OberservedScoreMatrix)) {
  cat("No sample to term mapping is specified. Using example mapping file ExampleData.rds.\n")
  argsList$OberservedScoreMatrix<-"Data/ExampleData.rds"
}

if(is.null(argsList$ObservedSimMethod)) {
  argsList$ObservedSimMethod<-"spearman"
}

if(is.null(argsList$OntologySimMethod)) {
  argsList$OntologySimMethod<-"spearman"
}

if(is.null(argsList$Output)) {
  cat("No output specified, using default Ontology_Score_Output.txt\n")
  argsList$Output<-"Ontology_Score_Output.txt"
}

if(is.null(argsList$fontsize)) {
  argsList$fontsize<-20
}

if(is.null(argsList$Log2)) {
  argsList$Log2<-TRUE
}

if(is.null(argsList$Center)) {
  argsList$Center<-TRUE
}

if(is.null(argsList$Scale)) {
  argsList$Scale<-TRUE
}

if(is.null(argsList$nPCA)) {
  argsList$nPCA<-as.numeric(4)
}
##########################
###Executing procedures###
##########################
cat("Loading ontology based similarities\n")
Ontology<-load.Ontology(filename=argsList$Ontology)
cat("Loading mapping files\n")
Idx2Term<-load.Idx2Term(filename=argsList$Idx2Term)
Sample2Term<-load.Sample2Term(filename=argsList$Sample2Term)
cat("Loading observed data matrix and converting it to a distance matrix\n")
OberservedScoreMatrix<-load.ObservedScoreMatrix(filename=argsList$OberservedScoreMatrix)
ontology.distance.matrix<-compute.Ontology.Distance(Ontology, Idx2Term, Sample2Term)
observed.distance.matrix<-compute.Observed.Distance(OberservedScoreMatrix,Log2=argsList$Log2,Center=argsList$Center,Scale=argsList$Scale,nPCA=argsList$nPCA,method=argsList$ObservedSimMethod)
cat("Computing ontology score\n")
ontology.score<-compute.Ontology.Score(ontology.distance.matrix, observed.distance.matrix, Sample2Term, method=argsList$OntologySimMethod)
save.ontology.score(ontology.score, filename = argsList$Output)
if (!is.null(argsList$pngFile)){
	if (ggplotAvailable){
		 p<-generateBoxPlot(ontology.score, fontsize= argsList$fontsize, pngFile = argsList$pngFile)
	}else{
	print("ggplot2 is not installed, boxplots are not generated")
	}
}


