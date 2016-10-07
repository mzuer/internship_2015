###############################################################
# This file can be used to plot individuals figures           #
# but you may need to change folders' path and files' name    #
###############################################################
# Functions and scripts to plot figures without running all the script ("script.R"), i.a.:
# for figures with continuous/discrete axis (with all stages or until birth), 
# for connectivity, for correlations between parameters, for correlations concerning broad/specific genes
# etc. Names of the files and folders should be changed manually.
rm(list=ls())
require(ggplot2)
library(gridExtra)
FPKM <- TRUE #FALSE

#choose the organisms to plot the results 
#and the folders where there are tables with the results to plot
# !work here if all files end with "PartialspearmanCorStagesBgeeOriginal.txt")
organisms <- c("Mus", "Danio", "Xenopus", "Drosophila", "Caeno")#, "Danio", "Xenopus", "Drosophila", "Caeno")
# organisms <- c("Mus", "Danio", "Xenopus", "Drosophila", "Caeno")
# organisms <- c("Mus",  "Drosophila", "Caeno")  #for phyletic age
# organisms <- c("Danio")
if(!("Caeno" %in% organisms)){
  folder <- c("mouse", "zebrafish", "xenopus", "drosophila")  
#   folder <- c( "zebrafish")
}else{
  folder <- c("mouse", "zebrafish", "xenopus", "drosophila", "caeno")
#   folder <- c("mouse", "drosophila", "caeno")  #for phyetic age
}
folderAnalysis <- paste("~/Bureau/stage15/")  

#the variables that will store the dataframe (e.g. MusCorr)
name <- sapply(organisms, function(x) {x<-paste0(x, "Corr")}) 

#store the dataframe with the results in a variable for each organisms
for (i in organisms){
  k <- match(i, organisms)
  a <- read.table(paste(folderAnalysis, folder[k],"/analysis/", i, 
                                     "PartialspearmanCorStagesBgeeOriginal.txt",sep=""),
                          header=T)
  assign(paste(name[[k]]), a)
}

#####
# Plot according to stages: discrete axis -----------------------------------------------------------------
#####
ggplotCorrelations <- function(x){

    plotList <- list()
  
    if((x=="P.1"||x == "LRT") & "Caeno" %in% organisms){
      b <- length(name)-1
    }else{
      b <- length(name)
    }
    #define color for the x-axis according to stages of development "equivalence"

    for(i in 1:b){
        table <- paste0(name[i],x)  
        #name the table e.g. MusCorrOmega
        assign(table, get(name[i])[which(regexpr("stage", get(name[i])$x1)>0 & regexpr(x, get(name[i])$x2)>0),]) 
        #assign df with only stages for "x" parameter to the name e.g. MusCorrOmega
        
        data <- get(table)
        data$significant <- as.factor(data$significant) #need as factor for the colour below
        
        data$x1 <- factor(as.character(data$x1),levels=data$x1)
        
        #Color as follows :
        # green -> zygote, cleavage
        # orange -> blastula, gastrula
        # blue -> neurula
        # red -> organogenesis, segmentation
        # black -> post-embryonic
        if(regexpr("Mus", name[i])>0){
#           colLabel <- c(rep("green",3), rep("orange",3), rep("blue", 1), rep("red",2), rep("black", 1))
          colLabel <- c(rep("green",3), rep("orange",3), rep("blue", 1), rep("red",2))
          speciesName <- "M. musculus"
        }else if(regexpr("Danio", name[i])>0){
#           colLabel <- c(rep("green",1), rep("orange", 2), rep("blue", 2), rep("red",4), rep("black", 5))
          colLabel <- c(rep("green",1), rep("orange", 2), rep("blue", 2), rep("red",4))
          speciesName <- "D. rario"
        }else if(regexpr("Drosophila", name[i])>0){
#           colLabel <- c(rep("green",1), rep("orange", 2), rep("blue",2), rep("red", 2),rep("black", 2))
          colLabel <- c(rep("green",1), rep("orange", 2), rep("blue",2), rep("red", 2))
          speciesName <- "D. melanogaster"
        }else if(regexpr("Caeno", name[i])>0 & !FPKM){
#           colLabel <- c(rep("orange", 1), rep("blue",1), rep("red", 1),rep("black", 2))
          colLabel <- c(rep("orange", 1), rep("blue",1), rep("red", 1))
          speciesName <- "C. elegans"
        }else if(regexpr("Xenopus", name[i])>0){
          colLabel <- c(rep("green",5), rep("orange", 4), rep("blue",2), rep("red", 3))
          speciesName <- "X. tropicalis"
#         }else if(regexpr("Drosophila", name[i])>0 & FPKM){
#           colLabel <- c(rep("green",1), rep("orange",1), rep("blue", 1), rep("red",9), rep("black", 18))
        }else if (regexpr("Caeno", name[i])>0 & FPKM){
          colLabel <- c(rep("green",2), rep("orange",8), rep("blue", 4), rep("red", 10))
#           colLabel <- c(rep("green",2), rep("orange",8), rep("blue", 4), rep("red", 10), rep("black", 11))
          speciesName <- "C. elegans"
        }
        
        data <- data[c(1:length(colLabel)),]

        p<- ggplot(data, aes(x=x1, y=corValue,group=1,colour=significant)) +
          geom_point(size=4) + geom_line(colour="black")+
          scale_colour_manual(values=c("0"="chocolate","1"="green"),guide="none")+
          scale_y_continuous("Spearmans' (partial) correlation coefficient")+
          scale_x_discrete(paste0(x," - ",organisms[i]))+
          theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
                axis.text.y = element_text(colour="black"),
                axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
                axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
        theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
        
        plotList[[i]] <- p
          
     ggsave(filename=paste0(folderAnalysis,x,organisms[i],".png"))
    }

if(x!="P.1" & !("Caeno" %in% organisms)){
  pdf(paste0(folderAnalysis,x,"WithoutCaeno.pdf"), onefile=TRUE)
  for(i in 1:b){
    print(plotList[[i]])
  }
  dev.off()
} else{
  pdf(paste0(folderAnalysis,x,".pdf"), onefile=TRUE)
  for(i in 1:b){
    print(plotList[[i]])
  }
  dev.off()
 }
}
# ggplotCorrelations("Tau")
ggplotCorrelations("LRT")
ggplotCorrelations("Omega")
ggplotCorrelations("P.1")
ggplotCorrelations("Paralogs.Number")
ggplotCorrelations("X..GC.content")
ggplotCorrelations("CDS.Length")
ggplotCorrelations("Phyletic.Age")

####
# Plot correlations with continuous x-axis ---------------------------------------------------------------------
####

plotCorrContX <- function(x){
  
  plotList <- list()
  
  if(x=="P.1" & "Caeno" %in% organisms){
    b <- length(name)-1
  }else{
    b <- length(name)
  }
  #define color for the x-axis according to stages of development "equivalence"
  
  for(i in 1:b){
    table <- paste0(name[i],x)  
    #name the table e.g. MusCorrOmega
    assign(table, get(name[i])[which(regexpr("stage", get(name[i])$x1)>0 & regexpr(x, get(name[i])$x2)>0),]) 
    #assign df with only stages for "x" parameter to the name e.g. MusCorrOmega
    
    data <- get(table)
    data$significant <- as.factor(data$significant) #need as factor for the colour below
    
    data$x1 <- factor(as.character(data$x1),levels=data$x1)
    
    stagesName <- data$x1
    
    if(regexpr("Mus", name[i])>0){
      colLabel <- c(rep("green",3), rep("orange",3), rep("blue", 1), rep("red",2), rep("black", 1))
      speciesName <- "M. musculus"
    }else if(regexpr("Danio", name[i])>0){
      colLabel <- c(rep("green",1), rep("orange", 2), rep("blue", 2), rep("red",4), rep("black", 5))
      speciesName <- "D. rario"
    }else if(regexpr("Drosophila", name[i])>0){
      colLabel <- c(rep("green",1), rep("orange", 2), rep("blue",2), rep("red", 2),rep("black", 2))
      speciesName <- "D. melanogaster"
    }else if(regexpr("Caeno", name[i])>0 & !FPKM){
      colLabel <- c(rep("orange", 1), rep("blue",1), rep("red", 1),rep("black", 2))
      speciesName <- "C. elegans"
    }else if(regexpr("Xenopus", name[i])>0){
      colLabel <- c(rep("green",5), rep("orange", 4), rep("blue",2), rep("red", 3))
      speciesName <- "X. tropicalis"
      #         }else if(regexpr("Drosophila", name[i])>0 & FPKM){
      #           colLabel <- c(rep("green",1), rep("orange",1), rep("blue", 1), rep("red",9), rep("black", 18))
    }else if (regexpr("Caeno", name[i])>0 & FPKM){
      colLabel <- c(rep("green",2), rep("orange",8), rep("blue", 4), rep("red", 10), rep("black", 11))
      speciesName <- "C. elegans"
    }
    
    if(FPKM){
      data <- data[c(1:length(colLabel)),]
    }

    if(regexpr("Mus", name[i])>0){  #scale : dpc
      #stages : 
      data$x1 <- c(0, 1, 2, 3, 4, 4.5, 7.5, 8.5, 9.5, 10.5, 12, 14, 16, 18, 42)
    }else if(regexpr("Danio", name[i])>0){  #scale : hours of dvpmt at 28.5°C
      data$x1 <- c(0, 6, 8, 9, 10, 11.66, 16, 24, 30, 48, 96, 120, 336, 720)
    }else if(regexpr("Drosophila", name[i])>0){ #scale : minutes after fertilization
      data$x1 <- c(42.5, 110, 155, 275, 290, 450, 1200, 2190, 3630)
    }else if(regexpr("Xenopus", name[i])>0){  #hours post-fertilization
      data$x1 <- c(1.5, 2, 2.25, 2.75, 3, 5, 7, 8, 9, 17.5, 20.75, 32.5, 44.5, 66)
      #     }else if(regexpr("Drosophila", name[i])>0 & FPKM){ #scale : hours
      #       data$x1 <- c(seq(1, 23,2), 36.5, 60.5, 96)
    }else if(regexpr("Caeno", name[i])>0 & FPKM){ #scale : minutes
      data$x1 <- c(seq(0,240,30),seq(300, 720,30), 850)
    }
    
    data$x1 <- as.numeric(as.character(data$x1))
    
    p <- ggplot(data, aes(x=x1, y=corValue,group=1,colour=significant)) +
      geom_point(size=4) + geom_line(colour="black")+
      scale_colour_manual(values=c("0"="chocolate","1"="green"),guide="none")+
      scale_y_continuous("Spearmans' (partial) correlation coefficient")+
      scale_x_continuous(paste0(x," - ",organisms[i]), breaks=data$x1, labels=stagesName)+
      theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
            axis.text.y = element_text(colour="black"),
            axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
            axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
      theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
    
    plotList[[i]] <- p
    
    ggsave(filename=paste0(folderAnalysis,x,organisms[i],"Continuous.png"))
  }
  
  if(x!="P.1" & !("Caeno" %in% organisms)){
    pdf(paste0(folderAnalysis,x,"WithoutCaenoContinuous.pdf"), onefile=TRUE)
    for(i in 1:b){
      print(plotList[[i]])
    }
    dev.off()
  } else{
    pdf(paste0(folderAnalysis,x,"Continuous.pdf"), onefile=TRUE)
    for(i in 1:b){
      print(plotList[[i]])
    }
    dev.off()
  }
}
# organisms <- organisms[which(organisms!="Caeno")]
plotCorrContX("Omega")


#####
# Plot all parameters except stages ---------------------------------------------------------------------
#####
# finally, the next section was used for the figure with correlations between parameters 
# (this one is similar as the following, but without mean/median expression)

# if("Danio" %in% organisms || "Xenopus" %in% organisms){
#   phyleticAge <- FALSE  #all do not have phyletic age
# } else{
#   phyleticAge <- TRUE
# }
# no PA available for all organisms !

if("Caeno" %in% organisms){
  selectome <- FALSE
} else{
  selectome <- TRUE
}

param <- sapply(organisms, function(x) {x<-paste0(x, "Param")}) #the name of the df for each organism

for (i in 1:length(organisms)){
     if(selectome){
          assign(param[i],get(name[i])[which(regexpr("stage", get(name[i])$x1)<0 & 
                                         regexpr("stage", get(name[i])$x2)<0),])
      }else if(!selectome){
        assign(param[i],get(name[i])[which(regexpr("stage", get(name[i])$x1)<0 & 
                                           regexpr("stage", get(name[i])$x2)<0 &
                                           regexpr("LRT", get(name[i])$x1)<0 &
                                           regexpr("LRT", get(name[i])$x2)<0 &
                                           regexpr("P.1", get(name[i])$x1)<0 &
                                           regexpr("P.1", get(name[i])$x2)<0),])
    }
        
  temp <- get(param[i])
  rownames(temp) <- paste(temp$x1,temp$x2,sep="*")
  temp$x1 <- NULL
  temp$x2 <- NULL
  temp <- temp[,c("corValue","significant")]
  colnames(temp) <- sapply(colnames(temp),function(x){x<-paste0(x,".",organisms[i])})
  assign(param[i],temp)
  
    print(nrow(get(param[i]))) #all should have the same number
}

data <- data.frame(row.names=rownames(get(param[i])))

for (i in 1:length(organisms)){  
  data <- cbind(data,get(param[i]))
}

dataP <- data[,regexpr("significant",colnames(data))>0]
dataC <- data[,regexpr("cor",colnames(data))>0]

if("Caeno" %in% organisms){
  pdf(paste0(folderAnalysis,"ParamCorrelationWithoutLRTandP1.pdf"), onefile=TRUE)
}else{
  pdf(paste0(folderAnalysis,"ParamCorrelation.pdf"), onefile=TRUE)
}

start <- 0.1
step <- 0.5
lab <- rownames(dataC)

lab2 <- gsub("\\*", "\\*\n", lab)


xgrad <- seq(start, (nrow(dataC)-1)*step+start, step)
yminmax <- max(abs(range(dataC)))
ygrad <- seq(-round(yminmax,1),round(yminmax,1), 0.2)

plot(x=start,y=dataC[1,1], ylab="Spearman's (partial) correlation coefficient", xlab="",
     pch=0, xlim=c(0,((nrow(dataC)-1)*step+4)),col=ifelse(dataP[1,1]==0, "coral", "chartreuse4"), lwd=1.8,
     cex=1.5,ylim=(c(-yminmax,yminmax)+c(-0.015,0.015)),xaxt="none", yaxt="none") 

axis(1, col.axis="black", las=2, labels=lab2, at=xgrad, cex.axis=0.5)  
axis(2, col.axis="black", las=2, at=ygrad)  

#nrow-1 * 0.5 => 0.5 inbetween, 0.1 left and more right for the legend
#col : chartreuse4 if signif, blue if not signif

segments(0,0,max(xgrad),0,col="grey")
abline(v=xgrad, lty=2, col="grey")

pos_x <- 0.1 #increment +0.5 at each row
for(i in 1:nrow(dataC)){
  b <- ifelse(i==1, 2, 1)  #for the first row, we dont'need to plot the first point
  
  for(j in b:ncol(dataC)){
    points(x=pos_x, y=dataC[i,j], pch=0+(j-1), col=ifelse(dataP[i,j]==0, "coral", "chartreuse4"), 
           cex=1.5,lwd=1.5)
  }
  pos_x <- pos_x +0.5
}

legend("topright", legend = c(organisms),pch=0:(length(organisms)-1),pt.cex=0.8, cex=0.8, bty="n")
legend("bottomright", legend = c("significant", "not significant"),pch=c(0,0),
       col=c("chartreuse4","coral"),pt.cex=0.8, cex=0.8, bty="n")


dev.off()

#------------------------------------------------------------------------------------------------
#Plot CorBgeeOriginal
#------------------------------------------------------------------------------------------------
# ! needs to rename folder/file in order that it works
# use to plot the figure with the values for correlation between all parameters and mean/median expression 
# values for all organisms are represented (different symbol for each organism)

organisms <- c("Mus", "Danio",  "DrosophilaM", "DrosophilaR", "Xenopus", "Caeno")
folder <- c("mouse", "zebrafish", "drosophila", "drosophila","xenopus", "caeno")
folderAnalysis <- paste("~/Bureau/stage15/")  
selectome <-FALSE
param <- sapply(organisms, function(x) {x<-paste0(x, "CorrParam")}) 

#store the dataframe with the results in a variable for each organisms
for (i in organisms){
  k <- match(i, organisms)
  a <- read.table(paste(folderAnalysis, folder[k],"/analysis/", i, 
                        "PartialspearmanCorBgeeOriginal.txt",sep=""),
                  header=T)
  assign(paste(param[[k]]), a)

  if(!selectome){
    assign(param[i],get(param[i])[which(  regexpr("LRT", get(param[i])$x1)<0 &
                                         regexpr("LRT", get(param[i])$x2)<0 &
                                         regexpr("P.1", get(param[i])$x1)<0 &
                                         regexpr("P.1", get(param[i])$x2)<0),])
  }
temp <- get(param[i])
rownames(temp) <- paste(temp$x1,temp$x2,sep=" vs. ")
temp$x1 <- NULL
temp$x2 <- NULL
temp <- temp[,c("corValue","significant")]
colnames(temp) <- sapply(colnames(temp),function(x){x<-paste0(x,".",organisms[i])})
assign(param[i],temp)

print(nrow(get(param[i]))) #all should have the same number
}

data <- data.frame(row.names=rownames(get(param[i])))

for (i in 1:length(organisms)){  
  data <- cbind(data,get(param[i]))
}

dataP <- data[,regexpr("significant",colnames(data))>0]
dataC <- data[,regexpr("cor",colnames(data))>0]

if("Caeno" %in% organisms){
#   pdf(paste0(folderAnalysis,"ParamCorrelationWithoutLRTandP1.pdf"), onefile=TRUE)
  svg(paste0(folderAnalysis,"ParamCorrelationWithoutLRTandP1.svg"))
}else{
#   pdf(paste0(folderAnalysis,"ParamCorrelation.pdf"), onefile=TRUE)
  png(paste0(folderAnalysis,"ParamCorrelation.png"))
}

start <- 0.1
step <- 0.5
lab <- rownames(dataC)

xgrad <- seq(start, (nrow(dataC)-1)*step+start, step)
yminmax <- max(abs(range(dataC)))
ygrad <- seq(-round(yminmax,1),round(yminmax,1), 0.2)

op <- par(mar = c(10,4,4,2) + 0.1)
plot(x=start,y=dataC[1,1], ylab="Spearman's (partial) correlation coefficient", xlab="",
     pch=0, xlim=c(0,((nrow(dataC)-1)*step+4)),col=ifelse(dataP[1,1]==0, "coral", "chartreuse4"), lwd=1.8,
     cex=1.5,ylim=(c(-yminmax,yminmax)+c(-0.015,0.015)),xaxt="none", yaxt="none") 

axis(1, col.axis="black", las=2, labels=lab, at=xgrad, cex.axis=0.5)  
axis(2, col.axis="black", las=2, at=ygrad)  

#nrow-1 * 0.5 => 0.5 inbetween, 0.1 left and more right for the legend
#col : chartreuse4 if signif, blue if not signif

segments(0,0,max(xgrad),0,col="grey")
abline(v=xgrad, lty=2, col="grey")

pos_x <- 0.1 #increment +0.5 at each row
for(i in 1:nrow(dataC)){
  b <- ifelse(i==1, 2, 1)  #for the first row, we dont'need to plot the first point
  
  for(j in b:ncol(dataC)){
    points(x=pos_x, y=dataC[i,j], pch=0+(j-1), col=ifelse(dataP[i,j]==0, "coral", "chartreuse4"), 
           cex=1.5,lwd=1.5)
  }
  pos_x <- pos_x +step
}
FPKM <- FALSE
organisms <- c("Mus", "Danio", "Drosophila", "Drosophila", "Xenopus", "Caeno")
speciesName <- unlist(lapply(organisms, function(x){setParamOrganism(x)$speciesName}))
symbRNA <- unlist(lapply(organisms, function(x){setParamOrganism(x)$symbRNA}))
FPKM <- TRUE;symbRNA[4] <- setParamOrganism(organisms[4])$symbRNA
labNames <- paste0(speciesName, " (", symbRNA, ")")
legend("topright", legend = c(labNames),pch=0:(length(organisms)-1),
       cex=0.5, pt.cex=0.5,  bty="n", text.font=3)
legend("bottomright", legend = c("significant", "not significant"),pch=c(0,0),
       col=c("chartreuse4","coral"),pt.cex=0.8, cex=0.5, bty="n")


dev.off()

par(op)

#------------------------------------------------------------------------------------------------
#Plot individual figure with discrete axis
#------------------------------------------------------------------------------------------------
plotDiscreteOne <- function(x, organism){
  
  data<-read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""),header=T)
  data <- data[which(regexpr("stage", data$x1)>0 & regexpr(x, data$x2)>0),] 
  data$significant <- as.factor(data$significant) #need as factor for the colour below
  
  colLabel <- graphParamAll(organism)$colLabel
  speciesName <- graphParamAll(organism)$speciesName
  symbRNA <- graphParamAll(organism)$symbRNA
  
  data <- data[c(1:length(colLabel)),]
  
  data$x1 <- factor(as.character(data$x1),levels=data$x1)
  
  tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA))
  
  x11()
  ggplot(data, aes(x=x1, y=corValue,group=1,colour=significant)) +
    geom_point(size=4) + geom_line(colour="black")+
    scale_colour_manual(values=c("0"="chocolate","1"="green"),guide="none")+
    scale_y_continuous("Spearmans' (partial) correlation coefficient")+
    ggtitle(tit)+
    scale_x_discrete(paste0(x))+
    theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
          axis.text.y = element_text(colour="black"),
          axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
    theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
  
  if(organism=="Drosophila" & FPKM){
    ggsave(filename=paste0(folderAnalysis,x,organism,"DiscreteFPKM.png"))
  } else{
    ggsave(filename=paste0(folderAnalysis,x,organism,"Discrete.png"))
  }
}
organism <- "Mus"
expDataSource <- "Bgee"
FPKM <- FALSE
folderAnalysis <- "~/Bureau/stage15/mouse/"

parameters <- c("CDS.Length", "Omega", "P.1", "LRT", "Paralogs.Number", "Phyletic.Age")
# parameters <- c("CDS.Length", "Omega", "Paralogs.Number", "Phyletic.Age") #for Caeno

for(i in parameters)
  plotDiscreteOne(i, organism)  #x, organism

#------------------------------------------------------------------------------------------------
#Plot individual figure with continuous x-axis: Omega with only broad/specific genes
#------------------------------------------------------------------------------------------------
source("scriptFunctions.R")
organism <-"Drosophila"     # "Danio" 
folder <- "drosophila/"   #"mouse/" # "zebrafish/"
folderAnalysis <- paste0("~/Bureau/stage15/",folder)  
part <- "Partial"
corMethod <- "spearman"
expDataSource <- "Bgee"
x <- "Omega"
FPKM <- FALSE
type <- "Broad"

data<-read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, tolower(type),".txt",sep=""),header=T)

data <- data[which(regexpr("stage", data$x1)>0 & regexpr(x, data$x2)>0),] 

data$significant <- as.factor(data$significant) #need as factor for the colour below

data$x1 <- factor(as.character(data$x1),levels=data$x1)

stagesName <- data$x1
    
colLabel <- graphParamConti(organism)$colLabel
speciesName <- graphParamConti(organism)$speciesName
stageNum <- graphParamConti(organism)$stageNum
symbRNA <- graphParamConti(organism)$symbRNA

data <- data[c(1:length(colLabel)),]
stagesName <- stagesName[c(1:length(colLabel))]
data$x1 <- stageNum

data$x1 <- as.numeric(as.character(data$x1))
tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA))
x11()
p <- ggplot(data, aes(x=x1, y=corValue,group=1,colour=significant)) +  geom_point(size=4) + geom_line(colour="black")+
  scale_colour_manual(values=c("0"="chocolate","1"="green"),guide="none")+
  scale_y_continuous("Spearmans' (partial) correlation coefficient")+
  scale_x_continuous(type, breaks=data$x1, labels=stagesName)+
  ggtitle(tit)+
  theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
  theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())

  ggsave(filename=paste0(folderAnalysis,x,organism,type,"Conti.png"))

#------------------------------------------------------------------------------------------------
#Plot connectivity with continuous axis
#------------------------------------------------------------------------------------------------
source("scriptFunctions.R")
folder <- "drosophila/"   #"mouse/" # "zebrafish/"
part <- "Partial"
corMethod <- "spearman"
x <- "Connectivity"
expDataSource <- "Bgee"

plotContinuousOne <- function(x, organism){

  data<-read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "OriginalWithConnect.txt",sep=""),header=T)
  
  data <- data[which(regexpr("stage", data$x1)>0 & regexpr(x, data$x2)>0),] 
  
  data$significant <- as.factor(data$significant) #need as factor for the colour below
  
  data$x1 <- factor(as.character(data$x1),levels=data$x1)
  
  stagesName <- data$x1
  
  colLabel <- graphParamConti(organism)$colLabel
  speciesName <- graphParamConti(organism)$speciesName
  stageNum <- graphParamConti(organism)$stageNum
  symbRNA <- graphParamConti(organism)$symbRNA
  
  data <- data[c(1:length(colLabel)),]
  stagesName <- stagesName[c(1:length(colLabel))]
  data$x1 <- stageNum
  data$x1 <- as.numeric(as.character(data$x1))
  
  tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA))
  x11()
  p <- ggplot(data, aes(x=x1, y=corValue,group=1,colour=significant)) +
    #   geom_point(aes(colour=significant),size=4) + geom_line(colour="black")+
    geom_point(size=4) + geom_line(colour="black")+
    #   scale_colour_gradient(low="white", high="chartreuse4",guide="none")+
    #           scale_colour_manual(breaks=c("0","1"), values=c("chocolate","green"),guide="none")+
    scale_colour_manual(values=c("0"="chocolate","1"="green"),guide="none")+
    scale_y_continuous("Spearmans' (partial) correlation coefficient")+
    scale_x_continuous(x, breaks=data$x1, labels=stagesName)+
    ggtitle(tit)+
    theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
          axis.text.y = element_text(colour="black"),
          axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
    theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
  
  if(organism=="Drosohpila" & FPKM){
    ggsave(filename=paste0(folderAnalysis,x,organism,"ConnectFPKM.png"))
  }else{
    ggsave(filename=paste0(folderAnalysis,x,organism,"Connect.png"))
  }
}
#To change before calling the function :
organism <-"Mus" #Mus"# "Drosophila"
expDataSource <- "Bgee"
FPKM <- FALSE
folderAnalysis <- "~/Bureau/stage15/mouse/"

parameters <- c("CDS.Length", "Omega", "P.1", "LRT", "Paralogs.Number", "Phyletic.Age")
# parameters <- c("CDS.Length", "Omega", "Paralogs.Number", "Phyletic.Age") #for Caeno
i="Connectivity"
# for(i in parameters)
  plotContinuousOne(i, organism)  #x, organism

#------------------------------------------------------------------------------------------------
#Plot individual figure with continuous axis
#------------------------------------------------------------------------------------------------
source("scriptFunctions.R")
organism <-"Drosophila"     # "Danio" 
folder <- "drosophila/analysis_ComBat_01.09/"   #"mouse/" # "zebrafish/"
folderAnalysis <- paste0("~/Bureau/stage15/",folder)  
part <- "Partial"
corMethod <- "spearman"
x <- "Omega"
expDataSource <- "Bgee"
FPKM <- FALSE

plotContinuousOne <- function(x, organism){
  
  data<-read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""),header=T)
  
  data <- data[which(regexpr("stage", data$x1)>0 & regexpr(x, data$x2)>0),] 
  
  data$significant <- as.factor(data$significant) #need as factor for the colour below
  
  data$x1 <- factor(as.character(data$x1),levels=data$x1)
  
  stagesName <- data$x1
  
  colLabel <- graphParamConti(organism)$colLabel
  speciesName <- graphParamConti(organism)$speciesName
  stageNum <- graphParamConti(organism)$stageNum
  symbRNA <- graphParamConti(organism)$symbRNA
  
  data <- data[c(1:length(colLabel)),]
  stagesName <- stagesName[c(1:length(colLabel))]
  data$x1 <- stageNum
  data$x1 <- as.numeric(as.character(data$x1))
  
  tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA))
  x11()
  p <- ggplot(data, aes(x=x1, y=corValue,group=1,colour=significant)) +
    geom_point(size=4) + geom_line(colour="black")+
    scale_colour_manual(values=c("0"="chocolate","1"="green"),guide="none")+
    scale_y_continuous("Spearmans' (partial) correlation coefficient")+
    scale_x_continuous(x, breaks=data$x1, labels=stagesName)+
    ggtitle(tit)+
    theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
          axis.text.y = element_text(colour="black"),
          axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
    theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
  
  if(organism=="Drosohpila" & FPKM){
    ggsave(filename=paste0(folderAnalysis,x,organism,"ContinuousFPKM.png"))
  }else{
    ggsave(filename=paste0(folderAnalysis,x,organism,"Continuous.png"))
  }
}
#To change before calling the function :
organism <-"Drosophila" #Mus"# "Drosophila"
expDataSource <- "Bgee"
FPKM <- FALSE
folderAnalysis <- "~/Bureau/stage15/drosophila/analysis_ComBat_01.09/"

parameters <- c("CDS.Length", "Omega", "P.1", "LRT", "Paralogs.Number", "Phyletic.Age")
# parameters <- c("CDS.Length", "Omega", "Paralogs.Number", "Phyletic.Age") #for Caeno
for(i in parameters)
  plotContinuousOne(i, organism)  #x, organism


#------------------------------------------------------------------------------------------------
#plot with >1 parameter
#------------------------------------------------------------------------------------------------
#use to plot variation of the correlation coefficients for multiple parameters
#but bad results because the values are very different (it would be better to plot "relative" values)

organism <-"Xenopus"     # "Danio" 
folder <- "xenopus/analysis_14.08_PA/"   #"mouse/" # "zebrafish/"
folderAnalysis <- paste0("~/Bureau/stage15/",folder)  
part <- "Partial"
corMethod <- "spearman"
expDataSource <- "BgeeRNA"
x <- "Paralogs.Number"
y <- "Phyletic.Age"
z <- "Omega"
FPKM <- F

data<-read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""),header=T)

data1 <- data[which(regexpr("stage", data$x1)>0 & regexpr(x, data$x2)>0),] 
data1 <- data1[,c("x1", "x2", "corValue", "significant")]

data2 <- data[which(regexpr("stage", data$x1)>0 & regexpr(y, data$x2)>0),] 
data2 <- data2[,c("x1", "x2", "corValue", "significant")]

data3 <- data[which(regexpr("stage", data$x1)>0 & regexpr(z, data$x2)>0),] 
data3 <- data3[,c("x1", "x2", "corValue", "significant")]

data <- rbind(data1, data2, data3)

temp <- data

data$significant <- as.factor(data$significant) #need as factor for the colour below

if(organism=="Mus"){
  colLabel <- c(rep("green",3), rep("orange",3), rep("blue", 1), rep("red",5), rep("black", 3))
  speciesName <- "M. musculus"
}else if(organism=="Danio"){
  colLabel <- c(rep("green",1), rep("orange", 2), rep("blue", 2), rep("red",4), rep("black", 5))
  speciesName <- "D. rerio"
}else if(organism=="Drosophila" & !FPKM){
  colLabel <- c(rep("green",1), rep("orange", 2), rep("blue",2), rep("red", 2),rep("black", 2))
  speciesName <- "D. melanogaster"
}else if(organism=="Caeno" & !FPKM){
  colLabel <- c(rep("orange", 1), rep("blue",1), rep("red", 1),rep("black", 2))
  speciesName <- "C. elegans"
}else if(organism=="Xenopus"){
  colLabel <- c(rep("green",5), rep("orange", 4), rep("blue",2), rep("red", 3))
  speciesName <- "X. tropicalis"
}else if(organism=="Drosophila" & FPKM){
  colLabel <- c(rep("green",1), rep("orange",1), rep("blue", 1), rep("red",9), rep("black", 15))
  speciesName <- "D. melanogaster"
} else if(organism=="Caeno" & FPKM){
  colLabel <- c(rep("green",2), rep("orange",8), rep("blue", 4), rep("red",10), rep("black", 11))
  speciesName <- "C. elegans"
}

data$x1 <- factor(as.character(data$x1),levels=data$x1)
data$x2 <- factor(as.character(data$x2),levels=data$x2)

colnames(data) <- c("x1", "Parameter", "corValue", "Significant")

x11()
ggplot(data, aes(x=x1, y=corValue,group=Parameter, shape=Significant), guide="none") +
  geom_point(size=3)+
  geom_line(aes(colour=Parameter))+
  scale_y_continuous("Spearmans' (partial) correlation coefficient")+
  scale_x_discrete(paste0(organism))+
  theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
  theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())

ggsave(filename=paste0(folderAnalysis, organism,"MultiParamPlot.png"))
