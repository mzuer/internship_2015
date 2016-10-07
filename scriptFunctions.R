###############################################################
# Functions used in the script:                               #
# Script for analysis Mus/Danio/Drosophila/Xenopus/Caeno data #
# Microarray (Affymetrix) & RNA-seq data                      #
# Gene expression during embryogenesis                        #
###############################################################
# These functions run when called from "script.R"
# (cannot be run independtly)
#-----------------------------------------------------------------------------------------
# Draw expression data distribution
#-----------------------------------------------------------------------------------------
# avoid repetition : small function for the graphs of expression distribution
# takes the title of the graph in argument

drawExpressionDistribution<-function(y){
  my.col <- colorRampPalette(c("#FFFFFF", "black", "blue", "#FA8072","#00A2FF", "#00CC00", "#E0E0E0"))(7) 
  dev.new(noRStudioGD = TRUE) #avoid warning in Rstudio
  #dev.new(height=9, width=12)
  #to open a window with linux in Rstudio
  x11(height=9, width=12)
  par(cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
  palette(rev(rich.colors(nDevStages+2)))
  
  plot(density(orgRPKM[,devStageNames[1]],n=1000),     
       main ="Expression values among different developmental stages",xlab=y,col=(1), lwd=3)  
  for(i in c(2:length(devStageNames))){
    lines(density(orgRPKM[,devStageNames[i]],n = 1000), col=(i), lwd=3)
  }
  legend("topright",devStagePrintNames,col=(1:length(devStageNames)),lty="solid", lwd=3)
}

#-----------------------------------------------------------------------------------------
# Data normalization (Cytoscape part)
#-----------------------------------------------------------------------------------------
# function to normalize the dataset
# takes the dataset in argument

normalizeData <-function(x){
     minOmega <- min(geneData$Omega[geneData$Omega>0])
     x$Omega <- x$Omega + minOmega
     x$Omega <- log2(x$Omega)
     
     minLength <- min(x$Intron.Length[x$Intron.Length>0])
     x$Intron.Length <- x$Intron.Length + minLength
     x$Intron.Length <- log2(x$Intron.Length)
     
     if(organism != "Caeno"){
       minP1 <- min(x$P.1[x$P.1>0])
       x$P.1 <- x$P.1 + minP1
       x$P.1 <- sqrt(sqrt(x$P.1))
       
       x$LRT <- ifelse(x$LRT<0, 0, x$LRT)
       x$LRT <- sqrt(sqrt(x$LRT))
     }
     
     x$Intron.Number <- x$Intron.Number + 1
     x$Intron.Number <- log2(x$Intron.Number)
     
     x$Paralogs.Number <- x$Paralogs.Number + 1
     x$Paralogs.Number <- log2(x$Paralogs.Number)
     
     if(phylAge){
      x$Phyletic.Age <- x$Phyletic.Age + 1
     }
     
     minTau <- 1 - max(x$Tau)
     x$Tau <- x$Tau + minTau
     x$Tau <- log2(x$Tau)
     
     x$CDS.Length <- log2(x$CDS.Length)
     
     return(x)
}

#-----------------------------------------------------------------------------------------
# Data normalization (circlize part)
#-----------------------------------------------------------------------------------------
# function to normalize the dataset
# takes the dataset in argument

normalizeData2 <- function(x){
  minOmega <- min(x$Omega[x$Omega>0])
  x$Omega <- x$Omega + minOmega
  x$Omega <- log2(x$Omega)
  
  minLength <- min(x$Intron.Length[x$Intron.Length>0])
  x$Intron.Length <- x$Intron.Length + minLength
  x$Intron.Length <- log2(x$Intron.Length)
  
  if(organism != "Caeno"){
    x$LRT <- ifelse(x$LRT<0, 0, x$LRT)
    x$LRT <- sqrt(sqrt(x$LRT))
    
    minP1 <- min(x$P.1[x$P.1>0])
    x$P.1 <- x$P.1 + minP1
    x$P.1 <- sqrt(sqrt(x$P.1))
  }
  
  x$Intron.Number <- x$Intron.Number + 1
  x$Intron.Number <- log2(x$Intron.Number)
  
  x$Paralogs.Number <- x$Paralogs.Number + 1
  x$Paralogs.Number <- log2(x$Paralogs.Number)
  
  if(phylAge){
    x$Phyletic.Age <- x$Phyletic.Age + 1
  }
  
  x$CDS.Length <- log2(x$CDS.Length)
  
  return(x)
}

#-----------------------------------------------------------------------------------------
# Functions for Circlize (chord diagram)
#-----------------------------------------------------------------------------------------
# easier to draw chord diagramm with Circlize if we have a symmetric matrix
# here we build a symmetric matrix with the dataframe of the correlations
# between stages and the other parameters
# takes the dataframe in argument

dfToMatrix<-function(dataF){ #takes the data frame as argument
  ##### Various functions
  #the function findValue -> returns the corValue for given x1, x2 and dataframe
  findValue<-function(x1,x2,dataF){
    if(x1==x2)   #e.g. stage1 with stage1 [diagonale not used for the chord diagram]
      return (1)
    for(i in 1:nrow(dataF)){
      if((dataF[i,1]==x1 & dataF[i,2]==x2)|(dataF[i,1]==x2 & dataF[i,2]==x1))
        return(dataF[i,3])
    }
    return (0)  #e.g.  stage1 with stage2 (is not in the input data)
  }
  
  #construct the matrix
  col1Names<-as.vector(dataF[,1])
  col2Names<-as.vector(dataF[,2])
  allVar<-unique(append(col1Names,col2Names))
  mat<-matrix(nrow=length(allVar),ncol=length(allVar))
  rownames(mat)<-allVar
  colnames(mat)<-allVar
  
  #full the matrix with the values from correlation table
  #on the diagonal : 1   (because correlation with itself ==1 but later we need 0, could be done here)
  #full by rows  
  i<-1
  while(i<=length(allVar)){
    j<-i
    while(j<=length(allVar)){
      mat[i,j]<-findValue(colnames(mat)[i],colnames(mat)[j],dataF)
      mat[j,i]<-mat[i,j]
      j<-j+1
    }
    i<-i+1
  }
  return(mat)
}


#-----------------------------------------------------------------------------------------
# Plot correlations between % identity with orthologous genes
#-----------------------------------------------------------------------------------------
# function to plot the %identity of orthologous genes pairs  between mouse and another organism
# takes the name of the organism to compare with mouse (-> used for the title)

plotOrthologous <- function(ortholog){
  for(i in devStageNames){
    ct <- cor.test(expOrtho[,i],expOrtho$X..Identity.with.respect.to.query.gene, method="pearson")
    coeff <- ct$estimate
    p.val <- ct$p.value
    corrCoeff <- append(coeff, corrCoeff)
    p.values <- append(p.values, p.val)
  }
  
  orthoTable <- data.frame(stage.name=devStageNames, corrValue=corrCoeff, pValue=p.values)
  orthoTable$significant <- ifelse(orthoTable$pValue < 0.05, 1, 0)
  
  #plot orthoTable 
  orthoTable$significant <- as.factor(orthoTable$significant) #need as factor for the colour below
  
  #set parameters using "class" defined elsewhere in this script
  colLabel <- graphParamAll(organism)$colLabel
  orthoTable$stage.name <- factor(as.character(orthoTable$stage.name),levels=orthoTable$stage.name)
  speciesName <- graphParamAll(organism)$speciesName
  
  plot.title = paste(speciesName, "/", ortholog)
  plot.subtitle = "orthologs"
  
  x11()
  ggplot(orthoTable, aes(x=stage.name, y=corrValue,group=1,colour=significant)) +
    geom_point(size=4) + geom_line(colour="black")+
    scale_colour_manual(values=c("0"="chocolate","1"="green"),guide="none")+
    scale_y_continuous("Pearson's correlation coefficient")+
    scale_x_discrete(paste0("% of identity between orthologous pairs"))+
  ggtitle(bquote(atop(italic(.(plot.title)), atop(.(plot.subtitle), ""))))+ 
    theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
          axis.text.y = element_text(colour="black"),
          axis.title.x = element_text(face="bold", colour="#CD5B45", size=12),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
    theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
  
  ggsave(filename=paste0(folderAnalysis,"X.orthology.with",ortholog, ".png"))
}

#-----------------------------------------------------------------------------------------
# Plot correlations for genes that are broad/specific to developmmental stages
#-----------------------------------------------------------------------------------------
# plot correlation coefficient values across developmental stages for broad/specific genes only
# takes in argument: dataset, parameter, type (e.g. data, "Omega", "specific")

plotByStages <- function (data, x, type) {  
  #argument : the dataset ("data"), the parameter to draw("x"), type:"broad" or "specific"
  #e.g.: plotByStages(specificCorr, "Omega", "specific")
  
  data <- data[which(regexpr("stage", data$x1)>0 & regexpr(x, data$x2)>0),] 
  data$significant <- as.factor(data$significant) #need as factor for the colour below
  symbRNA <- graphParamAll(organism)$symbRNA
  speciesName <- graphParamAll(organism)$speciesName
  tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA))
  
  if(organism=="Mus"){
    colLabel <- c(rep("green",3), rep("orange",3), rep("blue", 1), rep("red",5), rep("black", 3))
  }else if(organism=="Danio"){
    colLabel <- c(rep("green",1), rep("orange", 2), rep("blue", 2), rep("red",4), rep("black", 5))
  }else if(organism=="Drosophila" & !FPKM){
    colLabel <- c(rep("green",1), rep("orange", 2), rep("blue",2), rep("red", 2),rep("black", 2))
  }else if(organism=="Caeno" & !FPKM){
    colLabel <- c(rep("orange", 1), rep("blue",1), rep("red", 1),rep("black", 2))
  }else if(organism=="Xenopus"){
    colLabel <- c(rep("green",5), rep("orange", 4), rep("blue",2), rep("red", 3))
  }else if(organism=="Drosophila" & FPKM){
    colLabel <- c(rep("green",1), rep("orange",1), rep("blue", 1), rep("red",9), rep("black", 15))
  } else if(organism=="Caeno" & FPKM){
    colLabel <- c(rep("green",2), rep("orange",8), rep("blue", 4), rep("red",10), rep("black", 11))
  }
  
  data$x1 <- factor(as.character(data$x1),levels=data$x1)
  x11()
  ggplot(data, aes(x=x1, y=corValue,group=1,colour=significant)) +
    geom_point(size=4) + geom_line(colour="black")+
    scale_colour_manual(values=c("0"="chocolate","1"="green"),guide="none")+
    scale_y_continuous("Spearmans' (partial) correlation coefficient")+
    scale_x_discrete(paste0(type))+
    ggtitle(tit)+
    theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
          axis.text.y = element_text(colour="black"),
          axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
    theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
  ggsave(paste0(folderAnalysis, organism, x,type,organism,".png"))
  #     dev.off()
}

#-----------------------------------------------------------------------------------------
# Partial correlations
#-----------------------------------------------------------------------------------------
# finally, not used
#       parCorrF<-function(varNames,varToUse,dataset){
#         for(i in varNames){
#           for(j in varToUse){      
#             par<-varToUse[varToUse!=j]      #paramaters that should appear on the right part of the formula
#             pred<-par[1]
#             par<-par[-1]
#             for(k in par){
#               pred<-paste(pred,k,sep="+")
#             }    
#             ct<-cor.test(resid(glm(paste(i,"~",pred,sep=""),data=dataset,na.action=na.exclude)),
#                       resid(glm(paste(j,"~",pred,sep=""),data=set,na.action=na.exclude)),method=corMethod)
#           }
#         }
#         return(ct)
#       }

#-----------------------------------------------------------------------------------------
# Class for graph parameters (1)
#-----------------------------------------------------------------------------------------
# "class" to set parameters for an organism, with stages taken until birth only
# takes the name of the organism in argument

graphParamModel <- function(organism){
  if(organism == "Mus"){  #scale : dpc
    me <- list(
    stageNum = c(0, 1, 2, 3, 4, 4.5, 7.5, 8.5, 9.5, 10.5, 12, 14, 16, 18),
    sc = "Days post conception",
    speciesName = "M. musculus"
    )
  } else if(organism == "Danio"){  #scale : hours of dvpmt at 28.5°C
    me = list(
    stageNum = c(0, 6, 8, 9, 10, 11.66, 16, 24, 30, 48),
    sc = "Hours of development",
    speciesName = "D. rerio"
    )
  } else if(organism == "Drosophila" & !FPKM){ #scale : minutes after fertilization
    me <- list(
    stageNum = c(45, 105, 150, 270, 290, 440, 1200),
    sc = "Minutes after fertilization",
    speciesName = "D. melanogaster"
    )
  } else if(organism == "Xenopus"){  #hours post-fertilization
    me <- list(
    stageNum = c(1.5, 2, 2.25, 2.75, 3, 5, 7, 8, 9, 17.5, 20.75, 32.5, 44.5, 66),
    sc = "Hours post-fertilization",
    speciesName = "X. tropicalis"
    )
  } else if(organism == "Drosophila" & FPKM){  #hours post-fertilization
    me <- list(
    stageNum = c(seq(1, 23,2))*60, 
    sc = "Minutes after fertilization",
    speciesName = "D. melanogaster"
    )
  } else if(organism == "Caeno" & FPKM){  #hours post-fertilization
    me <- list(
    stageNum = c(seq(0,240,30),seq(300, 720,30)),
    sc = "Minutes after fertilization",
    speciesName = "C. elegans"
    )
  }
  class(me) <- append(class(me), "graphParamModel")
  return(me)
}

#-----------------------------------------------------------------------------------------
# Class for graph parameters (2)
#-----------------------------------------------------------------------------------------
# "class" to set parameters for an organism, with stages also after birth 
# takes the name of the organism in argument

graphParamConti <- function(organism){  

    if(organism == "Mus"){  #scale : dpc
    me <- list(
      stageNum = c(0, 1, 2, 3, 4, 4.5, 7.5, 8.5, 9.5, 10.5, 12, 14, 16, 18, 42),
      sc = "Days post conception",
      speciesName = "M. musculus",
      colLabel=c(rep("green",3), rep("orange",3), rep("blue", 1), rep("red",7), rep("black", 1)),
      symbRNA = "M"
    )
  } else if(organism == "Danio"){  #scale : hours of dvpmt at 28.5°C
    me = list(
      stageNum = c(0, 6, 8, 9, 10, 11.66, 16, 24, 30, 48, 96, 120, 336, 720),
      sc = "Hours of development",
      speciesName = "D. rerio",
      colLabel = c(rep("green",1), rep("orange", 2), rep("blue", 2), rep("red",4), rep("black", 5)),
      symbRNA = "M"
    )
  } else if(organism == "Drosophila" & !FPKM){ #scale : minutes after fertilization
    me <- list(
      stageNum = c(45, 105, 150, 270, 290, 440, 1200, 2190, 3630),
      sc = "Minutes after fertilization",
      speciesName = "D. melanogaster",
      colLabel = c(rep("green",1), rep("orange", 2), rep("blue",2), rep("red", 2),rep("black", 2)),
      symbRNA = "M"
    )
  } else if(organism == "Xenopus"){  #hours post-fertilization
    me <- list(
      stageNum = c(1.5, 2, 2.25, 2.75, 3, 5, 7, 8, 9, 17.5, 20.75, 32.5, 44.5, 66),
      sc = "Hours post-fertilization",
      speciesName = "X. tropicalis",
      colLabel = c(rep("green",5), rep("orange", 4), rep("blue",2), rep("red", 3)),
      symbRNA = "R"
    )
  } else if(organism == "Drosophila" & FPKM){  #hours post-fertilization
    me <- list(
      stageNum = c(seq(1, 23,2),36.5, 60.5, 96)*60, 
      sc = "Minutes after fertilization",
      speciesName = "D. melanogaster",
      colLabel = c(rep("green",1), rep("orange",1), rep("blue", 1), rep("red",9), rep("black", 3)),
      symbRNA = "R"
    )
  } else if(organism == "Caeno" & FPKM){  #hours post-fertilization
    me <- list(
      stageNum = c(seq(0,240,30),seq(300, 720,30), 850),
      sc = "Minutes after fertilization",
      speciesName = "C. elegans",
      colLabel = c(rep("green",2), rep("orange",8), rep("blue", 4), rep("red",10), rep("black", 1)),
      symbRNA = "R"
    )
  }
  class(me) <- append(class(me), "graphParamConti")
  return(me)
}

#-----------------------------------------------------------------------------------------
# Class for graph parameters (3)
#-----------------------------------------------------------------------------------------
# takes the name of the organism in argument

graphParamAll <- function(organism){  
  
  if(organism == "Mus"){  #scale : dpc
    me <- list(
#       stageNum = c(0, 1, 2, 3, 4, 4.5, 7.5, 8.5, 9.5, 10.5, 12, 14, 16, 18, 42),
      sc = "Days post conception",
      speciesName = "M. musculus",
      colLabel=c(rep("green",3), rep("orange",3), rep("blue", 1), rep("red",7), rep("black", 1)),
      symbRNA = "M"
    )
  } else if(organism == "Danio"){  #scale : hours of dvpmt at 28.5°C
    me = list(
#       stageNum = c(0, 6, 8, 9, 10, 11.66, 16, 24, 30, 48, 96, 120, 336, 720),
      sc = "Hours of development",
      speciesName = "D. rerio",
      colLabel = c(rep("green",1), rep("orange", 2), rep("blue", 2), rep("red",4), rep("black", 5)),
      symbRNA = "M"
    )
  } else if(organism == "Drosophila" & !FPKM){ #scale : minutes after fertilization
    me <- list(
#       stageNum = c(45, 105, 150, 270, 290, 440, 1200, 2190, 3630),
      sc = "Minutes after fertilization",
      speciesName = "D. melanogaster",
      colLabel = c(rep("green",1), rep("orange", 2), rep("blue",2), rep("red", 2),rep("black", 2)),
      symbRNA = "M"
    )
  } else if(organism == "Xenopus"){  #hours post-fertilization
    me <- list(
#       stageNum = c(1.5, 2, 2.25, 2.75, 3, 5, 7, 8, 9, 17.5, 20.75, 32.5, 44.5, 66),
      sc = "Hours post-fertilization",
      speciesName = "X. tropicalis",
      colLabel = c(rep("green",5), rep("orange", 4), rep("blue",2), rep("red", 3)),
      symbRNA = "R"
    )
  } else if(organism == "Drosophila" & FPKM){  #hours post-fertilization
    me <- list(
#       stageNum = c(seq(1, 23,2),36.5, 60.5, 96)*60, 
      sc = "Minutes after fertilization",
      speciesName = "D. melanogaster",
      colLabel = c(rep("green",1), rep("orange",1), rep("blue", 1), rep("red",9), rep("black", 15)),
      symbRNA = "R"
    )
  } else if(organism == "Caeno" & FPKM){  #hours post-fertilization
    me <- list(
#       stageNum = c(seq(0,240,30),seq(300, 720,30), 850),
      sc = "Minutes after fertilization",
      speciesName = "C. elegans",
      colLabel = c(rep("green",2), rep("orange",8), rep("blue", 4), rep("red",10), rep("black", 11)),
      symbRNA = "R"
    )
  }
  class(me) <- append(class(me), "graphParamAll")
  return(me)
}

#-----------------------------------------------------------------------------------------
# Class for parameters (1)
#-----------------------------------------------------------------------------------------
# "class" to set parameters for an organism
# also contains the stage names for mouse and drosophila (used when checking normalization)
# takes the name of the organism in argument

setParamOrganism <- function(organism){  
  if(organism == "Mus"){  #scale : dpc
    me <- list(
      sc = "Days post conception",
      speciesName = "M. musculus",
      symbRNA = "M",
      expData = "Bgee",
      stagesName = c("zygote stage", "Theiler stage 02 (mouse)", "Theiler stage 03 (mouse)",
                     "Theiler stage 04 (mouse)", "Theiler stage 05 (mouse)", "Theiler stage 06 (mouse)",
                     "neurula stage", "Theiler stage 13 (mouse)", "Theiler stage 15 (mouse)", 
                     "Theiler stage 17 (mouse)", "Theiler stage 20 (mouse)", "Theiler stage 22 (mouse)",
                     "Theiler stage 24 (mouse)","Theiler stage 26 (mouse)","post-juvenile adult stage")
    )
  } else if(organism == "Danio"){  #scale : hours of dvpmt at 28.5°C
    me = list(
      sc = "Hours of development",
      speciesName = "D. rerio",
      symbRNA = "M",
      expData = "Bgee"
    )
  } else if(organism == "Drosophila" & !FPKM){ #scale : minutes after fertilization
    me <- list(
      #       stageNum = c(45, 105, 150, 270, 290, 440, 1200, 2190, 3630),
      sc = "Minutes after fertilization",
      speciesName = "D. melanogaster",
      symbRNA = "M",
      expData = "Bgee",
      stagesName = c("embryonic stage 2 (Drosophila)", "embryonic stage 4 (Drosophila)",
                     "embryonic stage 5 (Drosophila)", "early extended germ band stage (Drosophila)",
                     "embryonic stage 10 (Drosophila)", "late extended germ band stage (Drosophila)", 
                     "embryonic stage 17 (Drosophila)", "first instar larval stage (Drosophila)", 
                     "second instar larval stage (Drosophila)")
    )
  } else if(organism == "Xenopus"){  #hours post-fertilization
    me <- list(
      sc = "Hours post-fertilization",
      speciesName = "X. tropicalis",
      symbRNA = "R",
      expData = "BgeeRNA"
    )
  } else if(organism == "Drosophila" & FPKM){  #hours post-fertilization
    me <- list(
      sc = "Minutes after fertilization",
      speciesName = "D. melanogaster",
      symbRNA = "R",
      expData = "Bgee"
    )
  } else if(organism == "Caeno"){  #hours post-fertilization
    me <- list(
      sc = "Minutes after fertilization",
      speciesName = "C. elegans",
      symbRNA = "R",
      expData = "Bgee"
    )
  }
  class(me) <- append(class(me), "setParamOrganism")
  return(me)
}

#-----------------------------------------------------------------------------------------
# Import after ComBat normalization (if 2 datasets only: for mouse and drosophila)
#-----------------------------------------------------------------------------------------
# function called for "orgExpression" when 2 datasets are used
# function loads the 2 datasets, performs ComBat normalization and returns
# the dataset in the form compatible with the rest of the script
# takes the name of the organism (-> "Mus" or "Drosophila" only) in argument

importAfterComBat <- function(organism){
  if(organism=="Mus"){
    tab1 <- read.delim("musEMEXP51.tsv") 
    tab2 <- read.delim("musEMTAB368.tsv")
  }else if(organism=="Drosophila"){
    tab1 <- read.delim("drosoEMTAB379.tsv")
    tab2 <- read.delim("drosoGSE3955.tsv")
  }
  data <- rbind(tab1, tab2)
  data <- na.omit(data)
  temp0 <- data 
  stageNames1 <- levels(tab1$Stage.name)  
  stageNames2 <- levels(tab2$Stage.name)  
  stageNames <- levels(data$Stage.name)  
  ## Reshape the data
  data <- data[,c("Gene.ID", "Stage.name",  "Chip.ID", "Log.of.normalized.signal.intensity")]
  colnames(data) <- c("Ensembl.Gene.ID", "Stage.name",  "Chip.ID", "Log.of.normalized.signal.intensity")
  data$Stage.ID <- paste(data$Stage.name, data$Chip.ID, sep=".")
  data <- data[,c("Ensembl.Gene.ID", "Stage.ID", "Log.of.normalized.signal.intensity")]
  data <- aggregate(Log.of.normalized.signal.intensity~Ensembl.Gene.ID+Stage.ID, FUN="mean",data=data)
  data <- dcast(data, Ensembl.Gene.ID~Stage.ID, value.var="Log.of.normalized.signal.intensity")
  temp1 <- data
  data <- na.omit(data)
  genes <- data$Ensembl.Gene.ID
  data <- as.matrix(data[,-1])

  # Expression matrix
  expData <- data
  
  # Design matrix
  func <- function(x){
    devS <-sub("(\\..*)", "",x)
    ifelse(devS %in% stageNames1, x <- 1, x <- 2)
  }
  set <-unlist(lapply(colnames(expData), func))
  nStages <- c(1:ncol(expData))
  design <- data.frame(columns=nStages, condition=set)
  adjData <- as.matrix(design) 
  
  # ComBat normalization : 
  comBatObj <- ComBat(expData, adjData[,2])
  comBatObj <- as.data.frame(comBatObj)
  
  stages <- unique(sub("(\\..*)", "",colnames(comBatObj)))
  stagesT <- sub("(\\(.*)","", stages) #e.g. if stage1 and stage11
  comBatData <- sapply(stagesT, function(x) rowMeans(comBatObj [, grep(x, names(comBatObj))]))
  colnames(comBatData) <- stages
  comBatData <- as.data.frame(comBatData)
  comBatData$Ensembl.Gene.ID <- genes
  ord <-  c(colnames(comBatData)[length(colnames(comBatData))], colnames(comBatData)[1:(length(colnames(comBatData))-1)])
  comBatData <- comBatData[,c(ord)]
  
 return (comBatData)

}


#-----------------------------------------------------------------------------------------
# PCA, heatmap and boxplot analysis when 1 dataset
#-----------------------------------------------------------------------------------------
# function used for fly and mouse microarray datasets
# PCA plot and heatmaps for : without, ComBat, ComBat+QN, QN (on each set separately)+ComBat normalization
# takes the name of the organism in argument

plotHeatmapPCAoneSet <- function(organism){
#   reshape2OLD <- "http://cran.r-project.org/src/contrib/Archive/reshape2/reshape2_1.4.tar.gz"
#   install.packages(reshape2OLD, repos=NULL, type="source")
# -> problem can occur with the melt function: it produces a dataframe with X1, X2 or Var1, Var2  
# depending of the version of reshape2 package. If you get error like "X1" not found, ensure that
# you have the latest version installed and/or that the function is not masked by a function from
# another package
#normally defined in "script.R"  
#   folder <- getwd()
#folderAnalysis : can be changed according to your folder's path 
#                 but already defined if "script.R" run from the beginning
  if(organism=="Danio"){
    data <- read.delim("danioETABM33.tsv")
    expData <- "Bgee"
    FPKM <-FALSE
#     folderAnalysis <- paste0(folder, "/zebrafish/") #can be changed according to your folder's path
  }else if(organism=="Xenopus"){
    data <- read.delim("xenoGSE37452.tsv")
    expData <- "BgeeRNA"
    FPKM <-FALSE
#     folderAnalysis <- paste0(folder, "/xenopus/") #can be changed according to your folder's path
  }else if(organism=="Drosophila" & FPKM){ #FPKM must be TRUE (set in "script.R")
    data <- read.delim("droso_FPKMs.txt")
    expData <- "FPKM"
#     folderAnalysis <- paste0(folder, "/drosophila/") #can be changed according to your folder's path
  }else if(organism=="Caeno"){
    data <- read.delim("caeno_FPKMs.txt")
    expData <- "FPKM"
    FPKM <-TRUE
#     folderAnalysis <- paste0(folder, "/caeno/") #can be changed according to your folder's path
  }
  symbRNA <- setParamOrganism(organism)$symbRNA  
  speciesName <- setParamOrganism(organism)$speciesName
  tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA)) 
  data <- na.omit(data)
  
  if(!FPKM){
    if(expData == "Bgee"){
      data <- data[,c("Gene.ID", "Stage.name", "Chip.ID", "Log.of.normalized.signal.intensity")]
      colnames(data) <- c("Ensembl.Gene.ID", "Stage.name", "Chip.ID", "Exp.Value")  
      if(organism=="Mus" || (organism=="Drosophila" & !FPKM)){
        tab1 <- tab1[,c("Gene.ID", "Stage.name", "Chip.ID", "Log.of.normalized.signal.intensity")]
        colnames(tab1) <- c("Ensembl.Gene.ID", "Stage.name", "Chip.ID", "Exp.Value")  
        tab2 <- tab2[,c("Gene.ID", "Stage.name", "Chip.ID", "Log.of.normalized.signal.intensity")]
        colnames(tab2) <- c("Ensembl.Gene.ID", "Stage.name", "Chip.ID", "Exp.Value")  
      }
    }else if(expData =="BgeeRNA"){
      data <- data[,c("Gene.ID", "Stage.name", "Library.ID", "RPKM")]
      colnames(data) <- c("Ensembl.Gene.ID", "Stage.name", "Library.ID", "Exp.Value")  
    }
    
    if(organism == "Drosophila"){
      newLevels <- c("stage5","stage7","stage2","stage8","stage9",
                     "stage4","stage1","stage3","stage6")
    }else if(organism == "Danio"){
      newLevels <- c("stage3","stage4","stage5","stage2","stage10","stage14",
                     "stage11","stage12","stage13","stage9","stage8","stage7",
                     "stage6","stage1")
    }else if(organism == "Mus"){
      newLevels <- c("stage15","stage2","stage3","stage4","stage5","stage6",
                     "stage1","stage7","stage8","stage9","stage10","stage11",
                     "stage12","stage13","stage14")
    }else if(organism == "Xenopus"){
      newLevels <- c("stage1","stage2","stage3","stage9","stage8","stage10",
                     "stage11","stage12","stage13","stage14","stage4","stage5",
                     "stage6","stage7")
    }
    
    levels(data$Stage.name) <- as.factor(newLevels)
    
    ### Reshape the data
    
    if(expData == "Bgee"){
      stageNames <- levels(data$Stage.name)  
      data$Stage.ID <- paste(data$Stage.name, data$Chip.ID, sep=".")
      data <- data[,c("Ensembl.Gene.ID", "Stage.ID", "Exp.Value")]
    } else if(expData == "BgeeRNA"){
      stageNames <- levels(data$Stage.name)
      data$Stage.ID <- paste(data$Stage.name, data$Library.ID, sep=".")
      data <- data[,c("Ensembl.Gene.ID", "Stage.ID", "Exp.Value")]
    }

    data <- aggregate(Exp.Value~Ensembl.Gene.ID+Stage.ID, FUN="mean",data=data)
    data <- dcast(data, Ensembl.Gene.ID~Stage.ID, value.var="Exp.Value")

    row.names(data) <- data$Ensembl.Gene.ID
    data$Ensembl.Gene.ID <- NULL
    
    replicates <- lapply(strsplit(colnames(data), split=".", fixed=TRUE),function(x){x[1]})
    stages <- unique(unlist(replicates))
    nRep <- table(unlist(replicates))
    
    repliNames <- vector()
    
    for(i in stages){
      rep <- nRep[names(nRep)==i]
      repli <- letters[1:rep]
      stage <- rep(i, rep)
      repliN <- paste(stage, repli, sep=".")
      repliNames <- append(repliNames,repliN)
    }
  } else if(FPKM){
    data <- aggregate(.~ gene_id, FUN="mean",data=data)
    row.names(data) <- data$gene_id
    data$gene_id <- NULL
    repliNames <- paste0("stage", seq(1:length(colnames(data))))
  }
  
  colnames(data) <- repliNames
  
  data <- na.omit(data)
  
  if(organism == "Mus"){
    data <- data[,c(19:length(colnames(data)), 1:18)]
  }else if(organism == "Danio"){
    data <- data[,c(11:length(colnames(data)),1:10)]
  }else if(organism == "Xenopus"){
    data <- data[,c(11:length(colnames(data)), 1:10)]
  }

  #*****
  # Heatmap correlations between replicates or stages, before normalization
  #*****
  
  data <- as.matrix(data)
  covM <- cor(data, use="complete")
  
  melted_covM <- melt(covM)

  x11()
  ggplot(data=melted_covM, aes(x=Var1, y=Var2, fill=value))+
    geom_tile()+
    ggtitle(tit)+
    scale_fill_gradient2(low = "aliceblue", high ="darkblue", mid = "steelblue1", 
                         midpoint = 0.5, limit = c(0,1), name="Correlation") +
    theme(axis.text.x  = element_text(angle=90, hjust=0, vjust=0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  ggsave(filename=paste0(folderAnalysis,organism,"HeatmapCorr.png"))
  
  #*****
  # PCA before normalization 
  #*****
  condition <- "Without normalization"
  # Plot with colors by stages
  pca1 <- prcomp(t(data))
  
  x <- pca1$x[,1]
  y <- pca1$x[,2]
  stages <- lapply(strsplit(row.names(pca1$x), split=".", fixed=TRUE),function(x){x[1]})
  colStages <- as.numeric(gsub("stage", "", stages))
  labx <- paste0(colnames(summary(pca1)$importance)[1], " (", round(summary(pca1)$importance[2,1],4)*100, "%)")
  laby <- paste0(colnames(summary(pca1)$importance)[2], " (", round(summary(pca1)$importance[2,2],4)*100, "%)")
  
  x11()
  png(paste0(folderAnalysis, organism, "plotPCA.png"))
  plot(x,y, col=colStages, pch=(colStages%%25), xlab=labx,ylab=laby,main=tit,
       sub=condition, col.sub="red")
  legend("center", unique(unlist(stages)), col=unique(colStages), pch=unique(colStages%%25),
         cex=0.6, bty="n")
  dev.off()
  
  #*****
  # PCA after normalization 
  #*****
  tempData <- data
  condition <- "With QN normalization"
  #Normalization
  data <- as.matrix(data)
  if(expData == "Bgee"){
    data <- normalize.quantiles(data)
  }else if(expData == "FPKM" || expData == "BgeeRNA"){
    data <- log2(data)
    data <- normalize.quantiles(data)
    data[data==-Inf] <- log2(1.0001)
  }
  colnames(data) <- colnames(tempData)
  
  # Plot with colors by stages
  pca1 <- prcomp(t(data))
  
  x <- pca1$x[,1]
  y <- pca1$x[,2]
  stages <- lapply(strsplit(row.names(pca1$x), split=".", fixed=TRUE),function(x){x[1]})
  colStages <- as.numeric(gsub("stage", "", stages))

  labx <- paste0(colnames(summary(pca1)$importance)[1], " (", round(summary(pca1)$importance[2,1],4)*100, "%)")
  laby <- paste0(colnames(summary(pca1)$importance)[2], " (", round(summary(pca1)$importance[2,2],4)*100, "%)")
  
  x11()
  png(paste0(folderAnalysis, organism, "plotPCAafterQN.png"))
  plot(x,y, col=colStages, pch=(colStages%%25), xlab=labx,ylab=laby,main=tit,
       sub=condition, col.sub="red")
  legend("center", unique(unlist(stages)), col=unique(colStages), pch=unique(colStages%%25),
         cex=0.6, bty="n")
  dev.off()
  
  
  #*****
  # Heatmap correlations between replicates after normalization 
  #*****
  
  data <- as.matrix(data)
  covM <- cor(data, use="complete")
  
  melted_covM <- melt(covM)
  x11()
    ggplot(data=melted_covM, aes(x=as.factor(Var1), y=as.factor(Var2), fill=value))+
    geom_tile()+
    ggtitle(tit)+
    scale_fill_gradient2(low = "aliceblue", high ="darkblue", mid = "steelblue1", 
                         midpoint = 0.5, limit = c(0,1), name="Correlation") +
    theme(axis.text.x  = element_text(angle=90, hjust=0, vjust=0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  # ggsave(filename=paste0(folder,organism,symbRNA, "HeatmapCorrAfterQN.png"))
  ggsave(filename=paste0(folderAnalysis,organism, "HeatmapCorrAfterQN.png"))
  dev.off()
}

#--------------------------------------------------i--------------------------------------
# PCA, heatmap and boxplot analysis when 2 datasets
#-----------------------------------------------------------------------------------------
# function used for fly and mouse microarray datasets
# PCA plot and heatmaps for : without, ComBat, ComBat+QN, QN (on each set separately)+ComBat normalization
# takes the name of the organism in argument

plotHeatmapPCAtwoSets <- function(organism) {  
  # Set basic parameters (setParamOrganism = class defined elsewhere in this script)
  FPKM <- FALSE
  symbRNA <- setParamOrganism(organism)$symbRNA  
  speciesName <- setParamOrganism(organism)$speciesName
  expData <- setParamOrganism(organism)$expData
  stagesName <- setParamOrganism(organism)$stagesName
  tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA))
  # Import data
  #folder : can be changed according to your folder's path 
  #         but already defined if script.R run from the beginning
  if(organism=="Mus"){
    tab1 <- read.delim("musEMEXP51.tsv") #  (Affymetrix data1 from Bgee)
    tab2 <- read.delim("musEMTAB368.tsv") #  (Affymetrix data2 from Bgee)
    folder <- paste0(getwd(),"/mouse/")  #can be changed according to your folders' path
    stagesName <- sub(" \\(mouse\\)", "", stagesName)
  }else if(organism=="Drosophila"){
    tab1 <- read.delim("drosoEMTAB379.tsv")
    tab2 <- read.delim("drosoGSE3955.tsv")
    folder <- paste0(getwd(),"/drosophila/") #can be changed according to your folders' path
    stagesName <- sub(" \\(Drosophila\\)", "", stagesName)
  }
  data <- rbind(tab1, tab2)
  data <- na.omit(data)
  temp0 <- data 
  
  # Reshape the data
  if(expData == "Bgee"){
    stageNames <- levels(data$Stage.name)  
    data <- data[,c("Gene.ID", "Stage.name",  "Chip.ID", "Log.of.normalized.signal.intensity")]
    colnames(data) <- c("Ensembl.Gene.ID", "Stage.name",  "Chip.ID", "Log.of.normalized.signal.intensity")
  } 
  
  if(organism == "Drosophila"){
    newLevels <- c("stage5","stage7","stage2","stage8","stage9",
                   "stage4","stage1","stage3","stage6")
  }else if(organism == "Mus"){
    newLevels <- c("stage15","stage2","stage3","stage4","stage5","stage6",
                   "stage1","stage7","stage8","stage9","stage10","stage11",
                   "stage12","stage13","stage14")
  }
  
  levels1 <- newLevels[1:length(levels(tab1$Stage.name))]
  levels2 <- newLevels[(length(levels(tab1$Stage.name))+1):length(newLevels)]
  
  levels(data$Stage.name) <- as.factor(newLevels)
  data$Stage.ID <- paste(data$Stage.name, data$Chip.ID, sep=".")
  data <- data[,c("Ensembl.Gene.ID", "Stage.ID", "Log.of.normalized.signal.intensity")]
  data <- aggregate(Log.of.normalized.signal.intensity~Ensembl.Gene.ID+Stage.ID, FUN="mean",data=data)
  data <- dcast(data, Ensembl.Gene.ID~Stage.ID, value.var="Log.of.normalized.signal.intensity")
  
  colnames(data)<-sub("(stage[0-9]+)(\\.).*", "\\1",colnames(data))
  
  for(i in unique(colnames(data))){
    k = 1
    for(j in 1:length(colnames(data))){
      if(i==colnames(data)[j]){
        colnames(data)[j] <- paste(i, letters[k], sep=".")
        k = k+1
      }
    }
  }

  row.names(data) <- data$Ensembl.Gene.ID
  data$Ensembl.Gene.ID <- NULL
  
  temp <-sub("stage([0-9]+)(\\.).*", "\\1",colnames(data)[-1])
  temp <- as.character(sort(as.numeric(temp)))
  for(i in temp){
    k = 1
    for(j in 1:length(temp)){
      if(i==temp[j]){
        temp[j] <- paste(i, letters[k], sep=".")
        k = k+1
      }
    }
  }
  temp <- paste0("stage", temp)
  data <- data[,c(temp)]
  data <- as.matrix(data)
  
  #save matrix with expression values (raw values)
  write.table(data, file=paste0(folder, organism, "MatrixExpression_Extended.txt"))
  
  #*****
  # Raw data (log) : heatmap, PCA and boxplot 
  #*****
  data1 <- read.table(file=paste0(folder, organism, "MatrixExpression_Extended.txt"))
  data1 <- na.omit(data1)
  condition <- "Without normalization" #needed for the figure
  
  #### Heatmap
  data <- as.matrix(data1)
  covM <- cor(data, use="complete")
  melted_covM <- reshape2::melt(covM)
  x11()
  ggplot(data=melted_covM, aes(x=Var1, y=Var2, fill=value))+
    geom_tile()+
    ggtitle(tit)+
    scale_fill_gradient2(low = "aliceblue", high ="darkblue", mid = "steelblue1", 
                         midpoint = 0.5, limit = c(0,1), name="Correlation") +
    theme(axis.text.x  = element_text(angle=90, hjust=0, vjust=0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  ggsave(filename=paste0(folder,organism, "HeatmapCorrRawData.png"))
  
  #### PCA
  # Plot with colors by stages
  pca1 <- prcomp(t(data1))
  x <- pca1$x[,1]
  y <- pca1$x[,2]
  stages <- lapply(strsplit(row.names(pca1$x), split=".", fixed=TRUE),function(x){x[1]})
  colStages <- gsub("stage", "", stages)
  
  x11()
  png(paste0(folder, organism, "PCArawDataStages.png"))
  labx <- paste0(colnames(summary(pca1)$importance)[1], " (", round(summary(pca1)$importance[2,1],4)*100, "%)")
  laby <- paste0(colnames(summary(pca1)$importance)[2], " (", round(summary(pca1)$importance[2,2],4)*100, "%)")
  plot(x,y, col=colStages, pch=as.numeric(colStages), xlab=labx,ylab=laby, main=tit, 
       sub=condition, col.sub="red")
  legend("center", unique(unlist(stages)), col=unique(colStages), pch=unique(as.numeric(colStages)),
         cex=1, bty="n")
  dev.off()
  
  #Plot with colors by set (for Mus and Drosophila - Bgee data)
  colSet <- lapply(stages, function(x){
    if(x %in% levels1){
      y <- "1"
    } else if(x %in% levels2){
      y <- "2"
    } })
  
  colSet <- unlist(colSet)
  
  x11()
  png(paste0(folder, organism, "PCArawDataSet.png"))
  labx <- paste0(colnames(summary(pca1)$importance)[1], " (", round(summary(pca1)$importance[2,1],4)*100, "%)")
  laby <- paste0(colnames(summary(pca1)$importance)[2], " (", round(summary(pca1)$importance[2,2],4)*100, "%)")
  plot(x,y, col=colSet, pch=as.numeric(colSet), xlab=labx,ylab=laby, main=tit,
       sub=condition, col.sub="red")
  legend("center", c("set1", "set2"), col=unique(colSet), pch=unique(as.numeric(colSet)),
         cex=1, bty="n")
  dev.off()
  
  #### Boxplot
  data <- data1
  stages <- unique(sub("(stage[0-9]+)(\\.).*", "\\1",colnames(data)))
  stagesT <- paste0(stages, ".") #e.g. if stage1 and stage11
  aggData <- sapply(stagesT, function(x) rowMeans(data [, grep(x, names(data))]))
  colnames(aggData) <- stages
  
  x11()
  png(paste0(folder, organism, "BoxplotRawData.png"))
  boxplot(aggData, main=paste0(condition), ylab="Log of normalized signal intensity", las=2,
          names=stagesName, cex.axis=0.3, cex.main=0.8, cex.sub=0.4, sub=speciesName)
  dev.off()
  
  #*****
  # ComBat normalization : heatmap, PCA and boxplot 
  #*****
  data <- read.table(file=paste0(folder, organism, "MatrixExpression_Extended.txt"))
  expData <- na.omit(as.matrix(data))
  condition <- "With ComBat normalization"
    
  # Design matrix : no variable of interest, adjustment variable is set1/set2
  func <- function(x){
    devS <-sub("(stage[0-9]+)(\\.).*", "\\1",x)
    ifelse(devS %in% levels1, x <- 1, x <- 2)
  }
  set <-unlist(lapply(colnames(expData), func))
  nStages <- c(1:ncol(expData))
  design <- data.frame(columns=nStages, condition=set)
  design <- as.matrix(design) 
  write.table(design, file=paste0(folder, organism, "MatrixDesign_Extended.txt"))
  adjData <- read.table(paste0(folder, organism, "MatrixDesign_Extended.txt"))
  
  # ComBat normalization : 
  comBatObj <- ComBat(expData, adjData[,2])
  write.table(comBatObj, paste0(folder, organism, "ComBatExpression_Extended.txt"))
  data1 <- read.table(file=paste0(folder, organism, "ComBatExpression_Extended.txt"))
  data1 <- na.omit(data1)
  
  #### Heatmap
  data <- as.matrix(data1)
  covM <- cor(data, use="complete")
  melted_covM <- melt(covM)
  x11()
  ggplot(data=melted_covM, aes(x=Var1, y=Var2, fill=value))+
    geom_tile()+
    ggtitle(tit)+
    scale_fill_gradient2(low = "aliceblue", high ="darkblue", mid = "steelblue1", 
                         midpoint = 0.5, limit = c(0,1), name="Correlation") +
    theme(axis.text.x  = element_text(angle=90, hjust=0, vjust=0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  ggsave(filename=paste0(folder,organism, "HeatmapCorrComBat.png"))
  
  #### PCA
  # Plot with colors by stages
  pca1 <- prcomp(t(data1))
  x <- pca1$x[,1]
  y <- pca1$x[,2]
  stages <- lapply(strsplit(row.names(pca1$x), split=".", fixed=TRUE),function(x){x[1]})
  colStages <- gsub("stage", "", stages)
  
  x11()
  png(paste0(folder, organism, "PCAComBatStages.png"))
  labx <- paste0(colnames(summary(pca1)$importance)[1], " (", round(summary(pca1)$importance[2,1],4)*100, "%)")
  laby <- paste0(colnames(summary(pca1)$importance)[2], " (", round(summary(pca1)$importance[2,2],4)*100, "%)")
  plot(x,y, col=colStages, pch=as.numeric(colStages), xlab=labx,ylab=laby,main=tit, 
       sub=condition, col.sub="red")
  legend("center", unique(unlist(stages)), col=unique(colStages), pch=unique(as.numeric(colStages)),
         cex=1, bty="n")
  dev.off()
  
  #Plot with colors by set (for Mus and Drosophila - Bgee data)
  colSet <- lapply(stages, function(x){
    if(x %in% levels1){
      y <- "1"
    } else if(x %in% levels2){
      y <- "2"
    } })
  
  colSet <- unlist(colSet)
  
  x11()
  png(paste0(folder, organism, "PCAComBatSet.png"))
  labx <- paste0(colnames(summary(pca1)$importance)[1], " (", round(summary(pca1)$importance[2,1],4)*100, "%)")
  laby <- paste0(colnames(summary(pca1)$importance)[2], " (", round(summary(pca1)$importance[2,2],4)*100, "%)")
  plot(x,y, col=colSet, pch=as.numeric(colSet), xlab=labx,ylab=laby,main=tit,
       sub=condition, col.sub="red")
  legend("center", c("set1", "set2"), col=unique(colSet), pch=unique(as.numeric(colSet)),
         cex=1, bty="n")
  dev.off()
  
  #### Boxplot
  data <- data1
  stages <- unique(sub("(stage[0-9]+)(\\.).*", "\\1",colnames(data)))
  stagesT <- paste0(stages, ".") #e.g. if stage1 and stage11
  aggData <- sapply(stagesT, function(x) rowMeans(data [, grep(x, names(data))]))
  colnames(aggData) <- stages
  
  x11()
  png(paste0(folder, organism, "BoxplotComBat.png"))
  boxplot(aggData, main=paste0(condition), ylab="Log of normalized signal intensity", las=2,
          names=stagesName, cex.axis=0.3, cex.main=0.8, cex.sub=0.4, sub=speciesName)
  dev.off()
  
  #*****
  # ComBat + QN normalization : heatmap, PCA and boxplot 
  #*****
  data1 <- read.table(file=paste0(folder, organism, "ComBatExpression_Extended.txt"))
  data1 <- na.omit(data1)
  condition <- "With ComBat + QN normalization"
  
  # Quantile normalization
  expData <- as.matrix(data1)
  stages <- colnames(expData)
  genes <- row.names(expData)
  expData <- normalize.quantiles(expData)
  colnames(expData) <- stages
  row.names(expData) <- genes
  write.table(expData, paste0(folder, organism, "ComBatQN_Extended.txt"))
  data1 <- read.table(paste0(folder, organism, "ComBatQN_Extended.txt"))
  
  #### Heatmap
  data <- as.matrix(data1)
  covM <- cor(data, use="complete")
  melted_covM <- melt(covM)
  x11()
  ggplot(data=melted_covM, aes(x=Var1, y=Var2, fill=value))+
    geom_tile()+
    ggtitle(tit)+
    scale_fill_gradient2(low = "aliceblue", high ="darkblue", mid = "steelblue1", 
                         midpoint = 0.5, limit = c(0,1), name="Correlation") +
    theme(axis.text.x  = element_text(angle=90, hjust=0, vjust=0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  ggsave(filename=paste0(folder,organism, "HeatmapCorrComBatQN.png"))
  
  #### PCA
  # Plot with colors by stages
  pca1 <- prcomp(t(data1))
  x <- pca1$x[,1]
  y <- pca1$x[,2]
  stages <- lapply(strsplit(row.names(pca1$x), split=".", fixed=TRUE),function(x){x[1]})
  colStages <- gsub("stage", "", stages)
  
  x11()
  png(paste0(folder, organism, "PCAComBatQNstages.png"))
  labx <- paste0(colnames(summary(pca1)$importance)[1], " (", round(summary(pca1)$importance[2,1],4)*100, "%)")
  laby <- paste0(colnames(summary(pca1)$importance)[2], " (", round(summary(pca1)$importance[2,2],4)*100, "%)")
  plot(x,y, col=colStages, pch=as.numeric(colStages), xlab=labx, ylab=laby, main=tit,
       sub=condition, col.sub="red")
  legend("center", unique(unlist(stages)), col=unique(colStages), pch=unique(as.numeric(colStages)),
         cex=1, bty="n")
  dev.off()
  
  #Plot with colors by set (for Mus and Drosophila - Bgee data)
  colSet <- lapply(stages, function(x){
    if(x %in% levels1){
      y <- "1"
    } else if(x %in% levels2){
      y <- "2"
    } })
  
  colSet <- unlist(colSet)
  
  x11()
  png(paste0(folder, organism, "PCAComBatQNset.png"))
  labx <- paste0(colnames(summary(pca1)$importance)[1], " (", round(summary(pca1)$importance[2,1],4)*100, "%)")
  laby <- paste0(colnames(summary(pca1)$importance)[2], " (", round(summary(pca1)$importance[2,2],4)*100, "%)")
  plot(x,y, col=colSet, pch=as.numeric(colSet), xlab=labx,ylab=laby,main=tit,
       sub=condition, col.sub="red")
  legend("center", c("set1", "set2"), col=unique(colSet), pch=unique(as.numeric(colSet)),
         cex=1, bty="n")
  dev.off()
  
  #### Boxplot
  data <- data1
  stages <- unique(sub("(stage[0-9]+)(\\.).*", "\\1",colnames(data)))
  stagesT <- paste0(stages, ".") #e.g. if stage1 and stage11
  aggData <- sapply(stagesT, function(x) rowMeans(data [, grep(x, names(data))]))
  colnames(aggData) <- stages
  
  x11()
  png(paste0(folder, organism, "BoxplotComBatQN.png"))
  boxplot(aggData, main=paste0(condition), ylab="Log of normalized signal intensity", las=2,
          names=stagesName, cex.axis=0.3, cex.main=0.8, cex.sub=0.4, sub=speciesName)
  dev.off()
  
  #table to use for all analyses :
  colnames(aggData) <- setParamOrganism(organism)$stagesName
  write.table(aggData, file=paste0(folder,organism, "ComBatExpression.txt"))
  
  ##*****
  # QN (by set) + ComBat normalization : heatmap, PCA and boxplot 
  ##*****
  # Prepare set1
  set1 <- tab1[,c("Gene.ID", "Stage.name",  "Chip.ID", "Log.of.normalized.signal.intensity")]
  colnames(set1) <- c("Ensembl.Gene.ID", "Stage.name",  "Chip.ID", "Log.of.normalized.signal.intensity")
  levels(set1$Stage.name) <- as.factor(levels1)
  set1$Stage.ID <- paste(set1$Stage.name, set1$Chip.ID, sep=".")
  set1 <- set1[,c("Ensembl.Gene.ID", "Stage.ID", "Log.of.normalized.signal.intensity")]
  set1 <- aggregate(Log.of.normalized.signal.intensity~Ensembl.Gene.ID+Stage.ID, FUN="mean",data=set1)
  set1 <- dcast(set1, Ensembl.Gene.ID~Stage.ID, value.var="Log.of.normalized.signal.intensity")
  
  t <- colnames(set1)
  colnames(set1)<-sub("(stage[0-9]+)(\\.).*", "\\1",colnames(set1))
  for(i in unique(colnames(set1))){
    k = 1
    for(j in 1:length(colnames(set1))){
      if(i==colnames(set1)[j]){
        colnames(set1)[j] <- paste(i, letters[k], sep=".")
        k = k+1
      }
    }
  }
  t <- colnames(set1)
  row.names(set1) <- set1$Ensembl.Gene.ID
  set1$Ensembl.Gene.ID <- NULL
  
  # Build design matrix (for the moment = same as before)
  set <-unlist(lapply(colnames(expData), func))
  nStages <- c(1:ncol(expData))
  design <- data.frame(columns=nStages, condition=set)
  design <- as.matrix(design) 
  write.table(design, file=paste0(folder, organism, "MatrixDesign2sets_Extended.txt"))
  adjData <- read.table(paste0(folder, organism, "MatrixDesign2sets_Extended.txt"))
  
  temp <-sub("stage([0-9]+)(\\.).*", "\\1",colnames(set1)[-1])
  temp <- as.character(sort(as.numeric(temp)))
  for(i in temp){
    k = 1
    for(j in 1:length(temp)){
      if(i==temp[j]){
        temp[j] <- paste(i, letters[k], sep=".")
        k = k+1
      }
    }
  }
  temp <- paste0("stage", temp)
  set1 <- set1[,c(temp)]
  set1 <- as.matrix(set1)
  
  
  #Prepare set2
  set2 <- tab2[,c("Gene.ID", "Stage.name",  "Chip.ID", "Log.of.normalized.signal.intensity")]
  colnames(set2) <- c("Ensembl.Gene.ID", "Stage.name",  "Chip.ID", "Log.of.normalized.signal.intensity")
  levels(set2$Stage.name) <- as.factor(levels2)
  set2$Stage.ID <- paste(set2$Stage.name, set2$Chip.ID, sep=".")
  set2 <- set2[,c("Ensembl.Gene.ID", "Stage.ID", "Log.of.normalized.signal.intensity")]
  set2 <- aggregate(Log.of.normalized.signal.intensity~Ensembl.Gene.ID+Stage.ID, FUN="mean",data=set2)
  set2 <- dcast(set2, Ensembl.Gene.ID~Stage.ID, value.var="Log.of.normalized.signal.intensity")
  
  t <- colnames(set2)
  colnames(set2)<-sub("(stage[0-9]+)(\\.).*", "\\1",colnames(set2))
  for(i in unique(colnames(set2))){
    k = 1
    for(j in 1:length(colnames(set2))){
      if(i==colnames(set2)[j]){
        colnames(set2)[j] <- paste(i, letters[k], sep=".")
        k = k+1
      }
    }
  }
  t <- colnames(set2)
  row.names(set2) <- set2$Ensembl.Gene.ID
  set2$Ensembl.Gene.ID <- NULL
  
  temp <-sub("stage([0-9]+)(\\.).*", "\\1",colnames(set2)[-1])
  temp <- as.character(sort(as.numeric(temp)))
  for(i in temp){
    k = 1
    for(j in 1:length(temp)){
      if(i==temp[j]){
        temp[j] <- paste(i, letters[k], sep=".")
        k = k+1
      }
    }
  }
  temp <- paste0("stage", temp)
  set2 <- set2[,c(temp)]
  set2 <- as.matrix(set2)
  
  # Quantile normalization set1
  expData1 <- set1
  stages <- colnames(expData1)
  genes <- row.names(expData1)
  expData1 <- normalize.quantiles(expData1)
  colnames(expData1) <- stages
  row.names(expData1) <- genes
  
  # Quantile normalization set2
  expData2 <- set2
  stages <- colnames(expData2)
  genes <- row.names(expData2)
  expData2 <- normalize.quantiles(expData2)
  colnames(expData2) <- stages
  row.names(expData2) <- genes
  
  # Bind the two datasets
  expData <- merge(expData1, expData2, by="row.names")
  
  # Sort the columns
  temp <-sub("stage([0-9]+)(\\.).*", "\\1",colnames(expData)[-1])
  temp <- as.character(sort(as.numeric(temp)))
  for(i in temp){
    k = 1
    for(j in 1:length(temp)){
      if(i==temp[j]){
        temp[j] <- paste(i, letters[k], sep=".")
        k = k+1
      }
    }
  }
  temp <- paste0("stage", temp)
  expData <- expData[,c(temp)]
  expData <- as.matrix(expData)
  write.table(expData, paste0(folder, organism, "QNmerge2sets_Extended.txt"))
  
  # ComBat normalization : 
  comBatObj <- ComBat(expData, adjData[,2])
  write.table(comBatObj, paste0(folder, organism, "QNComBatExpression_Extended.txt"))
  data1 <- read.table(file=paste0(folder, organism, "QNComBatExpression_Extended.txt"))
  data1 <- na.omit(data1)
  condition <- "With QN (by set) + ComBat normalization"
    
  #### Heatmap
  data <- as.matrix(data1)
  covM <- cor(data, use="complete")
  melted_covM <- melt(covM)
  x11()
  ggplot(data=melted_covM, aes(x=Var1, y=Var2, fill=value))+
    geom_tile()+
    ggtitle(tit)+
    scale_fill_gradient2(low = "aliceblue", high ="darkblue", mid = "steelblue1", 
                         midpoint = 0.5, limit = c(0,1), name="Correlation") +
    theme(axis.text.x  = element_text(angle=90, hjust=0, vjust=0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  ggsave(filename=paste0(folder,organism, "HeatmapCorrQNComBat.png"))
  
  #### PCA
  # Plot with colors by stages
  pca1 <- prcomp(t(data1))
  x <- pca1$x[,1]
  y <- pca1$x[,2]
  stages <- lapply(strsplit(row.names(pca1$x), split=".", fixed=TRUE),function(x){x[1]})
  colStages <- gsub("stage", "", stages)
  
  x11()
  png(paste0(folder, organism, "PCAqnComBatStages.png"))
  labx <- paste0(colnames(summary(pca1)$importance)[1], " (", round(summary(pca1)$importance[2,1],4)*100, "%)")
  laby <- paste0(colnames(summary(pca1)$importance)[2], " (", round(summary(pca1)$importance[2,2],4)*100, "%)")
  plot(x,y, col=colStages, pch=as.numeric(colStages), xlab=labx, ylab=laby, main=tit,
       sub=condition, col.sub="red")
  legend("center", unique(unlist(stages)), col=unique(colStages), pch=unique(as.numeric(colStages)),
         cex=0.6, bty="n")
  dev.off()
  
  #Plot with colors by set (for Mus and Drosophila - Bgee data)
  colSet <- lapply(stages, function(x){
    if(x %in% levels1){
      y <- "1"
    } else if(x %in% levels2){
      y <- "2"
    } })
  
  colSet <- unlist(colSet)
  
  x11()
  png(paste0(folder, organism, "PCAqnComBatSet.png"))
  labx <- paste0(colnames(summary(pca1)$importance)[1], " (", round(summary(pca1)$importance[2,1],4)*100, "%)")
  laby <- paste0(colnames(summary(pca1)$importance)[2], " (", round(summary(pca1)$importance[2,2],4)*100, "%)")
  plot(x,y, col=colSet, pch=as.numeric(colSet), xlab=labx, ylab=laby, main=tit,
       sub=condition, col.sub="red")
  legend("center", c("set1", "set2"), col=unique(colSet), pch=unique(as.numeric(colSet)),
         cex=0.8, bty="n")
  dev.off()
  
  #### Boxplot
  data <- data1
  stages <- unique(sub("(stage[0-9]+)(\\.).*", "\\1",colnames(data)))
  stagesT <- paste0(stages, ".") #e.g. if stage1 and stage11
  aggData <- sapply(stagesT, function(x) rowMeans(data [, grep(x, names(data))]))
  colnames(aggData) <- stages
  
  x11()
  png(paste0(folder, organism, "BoxplotQNComBat.png"))
  boxplot(aggData, main=paste0(condition), ylab="Log of normalized signal intensity", las=2,
          names=stagesName, cex.axis=0.3, cex.main=0.8, cex.sub=0.4, sub=speciesName)
  dev.off()
}  
