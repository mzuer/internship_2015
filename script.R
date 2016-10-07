###############################################################
# Script for analysis Mus/Danio/Drosophila/Xenopus/Caeno data #
# Microarray (Affymetrix) & RNA-seq data                      #
###############################################################

# ---------------------------------------------------------------------------------------
# PREAMBLE
# ---------------------------------------------------------------------------------------

### Clear workspace 
rm(list = ls())

#####*****Choose working directory*****#####
#setwd("")
setwd("~/Bureau/stage15")

#####*****Used libraries*****#####
library("corpcor")          
library("plyr")             
library("reshape")          
library("gplots")           
library("lattice")          
library("latticeExtra")     
library("xtable")           
library("circlize")
library("RCytoscape")       
library("preprocessCore") 
library("reshape2") #to reshape dataframes
library("ggplot2")
library("stringr")
library("sva")
source("scriptFunctions.R") #some functions are define in another file

### Test Cytoscape connection
cy <- CytoscapeConnection()
pluginVersion(cy)

#####*****Choose general parameters*****#####
##### Color used for graphes
my.col <- colorRampPalette(c("#FFFFFF", "black", "blue", "#FA8072","#00A2FF", "#00CC00", "#E0E0E0"))(7) 
#1: Backgroundcolor for all graphs, #2: Foregroundcolor for all graphs (E6E6E6), #3: Fill for histograms,
#4: Red, for boxplots, #5: Blue, for boxplots, #6: Green, for boxplots, #7: Light gray

### Organism and data source :
organism <- "Drosophila" #"Danio" # Mus"#"Drosophila" #"Mus"   #Danio  #Drosophila

### modEncode data ? (for Drosophila and Caeno: if TRUE, values from Li et al. 2014)
if(organism == "Caeno" || organism == "Drosophila"){
  FPKM <- FALSE #TRUE or FALSE: choice for Caeno and Drosophila only
}else{
  FPKM <- FALSE  #FALSE only for other organisms
} 

### Correlation method
corMethod <- "spearman"
partial <- TRUE

### Different parameters depending of the organism
#(could also be set in a "classe" defined in "scriptFunctions.R" - not done yet)
if (organism == "Mus") {
	#folder <- paste("~/Folder with data for analysis")
  folder <- paste("~/Bureau/stage15/mouse/")
  #folderAnalysis <- paste("~/Folder with data for analysis")
  folderAnalysis <- paste("~/Bureau/stage15/mouse/")  
	organismName <- "Mouse"
	dataSource <- "EnsV80"
  expDataSource <- "Bgee"
  twoDataSets <- TRUE #need to be TRUE to perform ComBat normalization (2 datasets)
} else if(organism=="Danio"){
  #folder <- paste("~/Folder with data for analysis")
  folder <- paste("~/Bureau/stage15/zebrafish/")
  #folderAnalysis <- paste("~/Folder with data for analysis")
  folderAnalysis <- paste("~/Bureau/stage15/zebrafish/")  
  organismName <- "Zebrafish"
  dataSource <- "EnsV80"
  expDataSource <- "Bgee"
  twoDataSets <- FALSE
} else if(organism=="Drosophila"){
  #folder <- paste("~/Folder with data for analysis")
  folder <- paste("~/Bureau/stage15/drosophila/")
  #folderAnalysis <- paste("~/Folder with data for analysis")
  folderAnalysis <- paste("~/Bureau/stage15/drosophila/")  
  organismName <- "Drosophila"
  dataSource <- "EnsV80"
  expDataSource <- "Bgee"
  twoDataSets <- ifelse(FPKM, FALSE, TRUE)  #or: twoDataSets <- !FPKM (TRUE only if Bgee data: 2 datasets)
} else if (organism == "Xenopus") {
  #folder <- paste("~/Folder with data for analysis")
  folder <- paste("~/Bureau/stage15/xenopus/analysis/")
  #folderAnalysis <- paste("~/Folder with data for analysis")
  folderAnalysis <- paste("~/Bureau/stage15/xenopus/analysis/")  
  organismName <- "Xenopus"
  dataSource <- "EnsV80"
  expDataSource <- "BgeeRNA"
  twoDataSets <- FALSE
} else if (organism == "Caeno") {
  #folder <- paste("~/Folder with data for analysis")
  folder <- paste("~/Bureau/stage15/caeno/analysis/")
  #folderAnalysis <- paste("~/Folder with data for analysis")
  folderAnalysis <- paste("~/Bureau/stage15/caeno/analysis/")  
  organismName <- "Caeno"
  dataSource <- "EnsV80"
  expDataSource <- "Bgee"
  twoDataSets <- FALSE
} 

### For having values >1 (used later for normalization)
if (expDataSource == "Bgee" ) {
  expNorm <- 1
}else if (expDataSource == "BgeeRNA" ) {
  expNorm <- 10000  
}

### Name of the column is different depending RNA-seq or microarray data from Bgee
if (expDataSource == "Bgee" ) {
  values <- "Log.of.normalized.signal.intensity"
}else if (expDataSource == "BgeeRNA" ) {
  values <- "RPKM"
}

### Should phyletic age be included in the analysis ?
# from OGEE: no data for frog and zebrafish (1: mammalia, 6: cellular organism)
# from Ensembl (perl API): data for all (in MYA), so now always TRUE
if(organism=="Mus" || organism=="Drosophila" || organism=="Caeno"){
  phylAge <- TRUE #FALSE => include phyletic age in the analysis ? 
} else if(organism=="Danio" || organism == "Xenopus"){ 
  phylAge <- TRUE   # not the choice for Danio & Xenopus for OGEE data
}

### Should essentiality be included in the analysis ?
if(organism=="Mus" || organism=="Drosophila" || organism=="Caeno"){
  essentiality <- TRUE #FALSE => include phyletic age in the analysis ? For Mus, Drosophila only
} else if(organism=="Danio" || organism == "Xenopus"){ 
  essentiality <- FALSE   # not the choice for Danio & Xenopus for the moment
}

### Symbol (used later): R for RNA-seq (BgeeRNA or Li et al. 2014 data), M for microarray (Bgee)
if(expDataSource == "BgeeRNA" || FPKM){
  symbRNA <- "R"
}else{
  symbRNA <- "M"
}

# ---------------------------------------------------------------------------------------
# PREPARE THE INPUT
# ---------------------------------------------------------------------------------------

#####*****Input data*****#####
### Table with the age of the taxons (used for computing phyletic age)
ageTable <- read.table("ageTable.txt", sep=" ", header=T)

### Filters in Ensembl Biomart: 
#in Genes "protein_coding"
#Save: in CSV and Unique results only
#Data downloaded without the Transcript.ID

if(organism=="Mus"){
  ### Gene information
  #Ensembl.Gene.ID, Associated.Gene.Name, X..GC.content
  #orgGenes <- read.delim(paste("~/Gene strucuture data"))
  orgGenes <- read.delim(paste("MusGenesEnsV80.txt"),sep="\t") 

  ### Paralog information
  #Ensembl.Gene.ID, "Organism".Paralog.Ensembl.Gene.ID, Homology.Type
  #orgParalogs <- read.table(paste("~/Paralog information"), sep="\t", header=TRUE)
  orgParalogs <- read.table(paste("MusParalogsEnsV80.txt"), sep="\t", header=TRUE)
  #check import : colnames(orgParalogs)
  
  ### PPI information (protein-protein interaction) 
  #sciName, kingdom, taxID, protein, locus, scorecutoff, connectivity, toppercentile
  #The number of direct neighbors of genes in protein-protein network
  #orgConnectivity <- read.table(paste("~/PPI information"), sep="\t", header=TRUE)
  orgConnectivity <- read.table(paste("connectivityMus.txt"), sep="\t", header=TRUE)
  
  ### Protein gene connection
  #Ensembl.Gene.ID, Ensembl.Transcript.ID, Ensembl.Protein.ID
  # => need transcript ID for merging with connectivity
  #orgProtein <- read.table(paste("~/Protein information"), sep="\t", header=TRUE)
  orgProtein <- read.table(paste("MusProteinEnsV80.txt"), sep="\t", header=TRUE)
  
  if(phylAge){
        ### Phyletic ages of genes from OGEE
        #sciName, kingdom, locus, phyleticage, phyleticageTaxID, taxID, sciName.1, scorecutoff
         #orgPhyleticage <- read.table(paste("~/Phyletic age information"), sep="\t", header=TRUE)
         #orgPhyleticage <- read.table(paste("phyleticageMus.txt"), sep="\t", header=TRUE)

        # Phyletic age retrieved from Ensembl
        orgPhyleticage <- read.table("ageMus.txt", sep="\t")
  }
  
  if(essentiality){
  ### Essentiality of genes
    #omim_pheno_accession, omim_gene_descriptio, omim_morbid_description, mgi_accession, ens_mouse_id, 
    #ens_human_id, orthology_type, human_omim_desc_essentiality, mouse_ontology_essentiality
    #orgEssentiality <- read.table(paste("~/Essentiality information"), sep="\t", header=TRUE)
    orgEssentiality <- read.table(paste("Omim_essentiality_data_clinical_all_pheno.txt"), sep="\t", header=TRUE)
  }
  
  ### Omega from Selectome
    #Ensembl.Gene.ID, Ensembl.Transcript.ID, Omega.0, Omega.2, LRT, P.value, Q.value, 
    #Positive.Selection, P.0, P.1, P2 
    #orgSelectome <- read.table(paste("~/Evolutionary rate information"), sep="\t", header=TRUE)
    orgSelectome <- read.table(paste("Mouse_selectome_Euteleostomi.txt"), sep="\t", header=TRUE)  
    orgSelectome <- orgSelectome[,c("gene","omega0","lrt","selected","p0","p1")]
    colnames(orgSelectome) <- c("Ensembl.Gene.ID","Omega.0","LRT","Positive.Selection","P.0","P.1")
     
  ### Structure  
    #Ensembl.Gene.ID, Ensembl.Transcript.ID, CDS.Length, Exon.Rank.in.Transcript, Exon.Chr.Start..bp., 
    #Exon.Chr.End..bp.
    #Transcript ID needed ! If not transcript ID, I have many rows for 1 gene ID !!!
    #orgStructure <- read.table(paste("~/Gene structure information"), sep=",", header=TRUE)
    orgStructure <- read.table(paste("MusStructureEnsV80.txt"), sep="\t", header=TRUE)
  
  ### Gene expression 
    #orgExpression <- read.table(paste("~/Expression information"), sep="\t", header=TRUE)
    #orgExpression loaded later (ComBat)
    #tab1 <- read.delim("musEMEXP51.tsv") 
    #tab2 <- read.delim("musEMTAB368.tsv")
    #orgExpression <- rbind(tab1, tab2)

} else if(organism=="Danio"){
  ### Gene information
  #Ensembl.Gene.ID, Associated.Gene.Name, X..GC.content
  #orgGenes <- read.delim(paste("~/Gene strucuture data"))
  orgGenes <- read.delim(paste("DanioGenesEnsV80.txt"),sep="\t") 
  
  ### Paralog information
  #Ensembl.Gene.ID, "Organism".Paralog.Ensembl.Gene.ID, Homology.Type
  #orgParalogs <- read.table(paste("~/Paralog information"), sep="\t", header=TRUE)
  orgParalogs <- read.table(paste("DanioParalogsEnsV80.txt"), sep="\t", header=TRUE)
  
  ### PPI information (protein-protein interaction) 
  #sciName, kingdom, taxID, protein, locus, scorecutoff, connectivity, toppercentile
  #The number of direct neighbors of genes in protein-protein network
  #orgConnectivity <- read.table(paste("~/PPI information"), sep="\t", header=TRUE)
  orgConnectivity <- read.table(paste("connectivity.txt"), sep="\t", header=TRUE)
  orgConnectivity <- subset(orgConnectivity, regexpr("Danio rerio", sciName)>0)
  
  ### Protein gene connection
  #Ensembl.Gene.ID, Ensembl.Protein.ID
  # => need transcript ID for merging with connectivity
  #orgProtein <- read.table(paste("~/Protein information"), sep="\t", header=TRUE)
  orgProtein <- read.table(paste("DanioProteinEnsV80.txt"), sep="\t", header=TRUE)
  
 if(phylAge){            
  ### Phyletic ages of genes from OGEE
              #sciName, kingdom, locus, phyleticage, phyleticageTaxID, taxID, sciName.1, scorecutoff
              #orgPhyleticage <- read.table(paste("~/Phyletic age information"), sep="\t", header=TRUE)
              #orgPhyleticage <- read.table(paste("phyleticage.txt"), sep="\t", header=TRUE)
              #orgPhyleticage <- subset(orgPhyleticage, regexpr("Danio rerio", sciName)>0)
              # from OGEE: 0 data
            
  ### Phyletic age from Ensembl
            orgPhyleticage <- read.table("ageDanio.txt", sep="\t")
 }
  ### Essentiality of genes
if(essentiality){
  # orgEssentiality <- read.table(paste("~/Essentiality information"), sep="\t", header=TRUE)
  # from ogee, only the experimental data (taken) [text mining : only 1 gene...]
  # from the downloaded file, remove the "#" at the beginning of the header row
  # taxID, sciName, locus, essential, pubmedID, dataType, dataSource
  orgEssentiality<- read.table(paste("Ogee_essentiality_expDanio.txt"), sep="\t", header=TRUE)
  #only 306
}
  ### Omega from Selectome
  #Ensembl.Gene.ID, Ensembl.Transcript.ID, Omega.0, Omega.2, LRT, P.value, Q.value, 
  #Positive.Selection, P.0, P.1, P2 
  #orgSelectome <- read.table(paste("~/Evolutionary rate information"), sep="\t", header=TRUE)
  orgSelectome <- read.table(paste("Danio_rerio_selectome_Euteleostomi.txt"), sep="\t", header=TRUE)
  orgSelectome <- orgSelectome[,c("gene","omega0","lrt","selected","p0","p1")]
  colnames(orgSelectome) <- c("Ensembl.Gene.ID","Omega.0","LRT","Positive.Selection","P.0","P.1")
  
  
  ### Structure  
  #Ensembl.Gene.ID, Ensembl.Transcript.ID, CDS.Length, Exon.Rank.in.Transcript, Exon.Chr.Start..bp., 
  #Exon.Chr.End..bp.
  #Transcript ID needed ! If not transcript ID, I have many rows for 1 gene ID !!!
  #orgStructure <- read.table(paste("~/Gene structure information"), sep=",", header=TRUE)
  orgStructure <- read.table(paste("DanioStructureEnsV80.txt"), sep="\t", header=TRUE)

  ### Gene expression 
  #orgExpression <- read.table(paste("~/Expression information"), sep="\t", header=TRUE)
  orgExpression <- read.delim("danioETABM33.tsv") #  (Affymetrix data from Bgee)

  } else if(organism=="Drosophila"){
    
  ### Gene information
  #Gene.stable.ID, Gene.Name, X..GC.content
  #orgGenes <- read.delim(paste("~/Gene strucuture data"))
  orgGenes <- read.delim(paste("DrosoGenes.txt"),sep="\t") 
  orgGenes <- orgGenes[,c("Gene.stable.ID", "Gene.name", "X..GC.content")]
  colnames(orgGenes) <- c("Ensembl.Gene.ID", "Associated.Gene.Name", "X..GC.content")
  
  ### Paralog information
  #Gene.stable.ID, Paralog.gene.stable.ID, Homology.type
  #orgParalogs <- read.table(paste("~/Paralog information"), sep="\t", header=TRUE)
  orgParalogs <- read.table(paste("DrosoParalogs.txt"), sep="\t", header=TRUE)
  orgParalogs <- orgParalogs[,c("Gene.stable.ID", "Paralog.gene.stable.ID", "Homology.type")]
  colnames(orgParalogs) <- c("Ensembl.Gene.ID", "Paralog.Ensembl.Gene.ID","Homology.Type")
  
  ### PPI information (protein-protein interaction) 
  #sciName, kingdom, taxID, protein, locus, scorecutoff, connectivity, toppercentile
  #The number of direct neighbors of genes in protein-protein network
  #orgConnectivity <- read.table(paste("~/PPI information"), sep="\t", header=TRUE)
  orgConnectivity <- read.table(paste("connectivity.txt"), sep="\t", header=TRUE)
  orgConnectivity <- subset(orgConnectivity, regexpr("Drosophila melanogaster", sciName)>0) 
  
  ### Protein gene connection
  #Gene.stable.ID, Protein.stable.ID
  # => need transcript ID for merging with connectivity
  #orgProtein <- read.table(paste("~/Protein information"), sep="\t", header=TRUE)
  orgProtein <- read.table(paste("DrosoProtein.txt"), sep="\t", header=TRUE)
  #check if this is the good dataset : colnames(orgProtein)
  orgProtein <- orgProtein[,c("Gene.stable.ID","Protein.stable.ID")]
  colnames(orgProtein) <- c("Ensembl.Gene.ID", "Ensembl.Protein.ID")
  
  if(phylAge){
    ### Phyletic ages of genes from OGEE
    #sciName, kingdom, locus, phyleticage, phyleticageTaxID, taxID, sciName.1, scorecutoff
    #orgPhyleticage <- read.table(paste("~/Phyletic age information"), sep="\t", header=TRUE)
    #orgPhyleticage <- read.table(paste("phyleticage.txt"), sep="\t", header=TRUE)
    #orgPhyleticage <- subset(orgPhyleticage, regexpr("Drosophila melanogaster", sciName)>0) 
    
    # Phyletic age from Ensembl
    orgPhyleticage <- read.table("ageDroso.txt", sep="\t")
  }
  
  ### Essentiality of genes
  if(essentiality){
  # orgEssentiality <- read.table(paste("~/Essentiality information"), sep="\t", header=TRUE)
  # from ogee, only the experimental data (taken) [text mining : only 1 gene...]
  # from the downloaded file, remove the "#" at the beginning of the header row (before read.table in R)
  # taxID, sciName, locus, essential, pubmedID, dataType, dataSource
  orgEssentiality<- read.table(paste("Ogee_essentiality_expDroso.txt"), sep="\t", header=TRUE) 
}

  ### Omega from Selectome
  #Ensembl.Gene.ID, Ensembl.Transcript.ID, Omega.0, Omega.2, LRT, P.value, Q.value, 
  #Positive.Selection, P.0, P.1, P2 
  #orgSelectome <- read.table(paste("~/Evolutionary rate information"), sep="\t", header=TRUE)
  orgSelectome <- read.table(paste("Melanogaster_selectome_Drosophila.txt"), sep="\t", header=TRUE)
  orgSelectome <- orgSelectome[,c("gene","omega0","lrt","selected","p0","p1")]
  colnames(orgSelectome) <- c("Ensembl.Gene.ID","Omega.0","LRT","Positive.Selection","P.0","P.1")

  ### Structure  
  #Gene.stable.ID, Transcript.stable.ID, Exon.region.start..bp., Exon.region.end..bp., 
  # Exon.rank.in.Transcript, CDS.start..within.cDNA, CDS.end..within.cDNA
  ## !!! CDS length -> CDS end - CDS start
  #Transcript ID needed ! If not transcript ID, I have many rows for 1 gene ID !!!
  #orgStructure <- read.table(paste("~/Gene structure information"), sep=",", header=TRUE)
  orgStructure <- read.table(paste("DrosoStructure.txt"), sep="\t", header=TRUE)
  orgStructure$CDS.Length <- orgStructure$CDS.end..within.cDNA. - orgStructure$CDS.start..within.cDNA.
  temp<-(orgStructure[orgStructure$CDS.Length<0,]); nrow(temp);temp<-na.omit(temp);nrow(temp)
 
  orgStructure <- orgStructure[,c("Gene.stable.ID","Transcript.stable.ID", "CDS.Length",
                                  "Exon.rank.in.transcript","Exon.region.start..bp.", "Exon.region.end..bp.")]
  colnames(orgStructure) <- c("Ensembl.Gene.ID", "Ensembl.Transcript.ID", "CDS.Length", 
                              "Exon.Rank.in.Transcript", "Exon.Chr.Start..bp.","Exon.Chr.End..bp.")
  
  if(!FPKM){
  ### Gene expression 
  #orgExpression <- read.table(paste("~/Expression information"), sep="\t", header=TRUE)
  #orgExpression loaded later (ComBat)
  #tab1 <- read.delim("drosoEMTAB379.tsv") #  (Affymetrix data from Bgee, exp1)
  #tab2 <- read.delim("drosoGSE3955.tsv") # (Affymetrix data from Bgee, exp2)
  #orgExpression <- rbind(tab1, tab2)
  } else if (FPKM){
    orgExpression <- read.table("droso_FPKMs.txt",header=T)
    colnames(orgExpression)[1] <- "Ensembl.Gene.ID"
    orgExpression[,-1] <- as.numeric(as.character(unlist((orgExpression[,-1]))))
    orgExpression$Mean.1d <- rowMeans(subset(orgExpression, select=c(Male.1d, Female.1d), na.rm=TRUE))
    orgExpression$Mean.5d <- rowMeans(subset(orgExpression, select=c(Male.5d, Female.5d), na.rm=TRUE))                                             
    orgExpression$Mean.30d <- rowMeans(subset(orgExpression, select=c(Male.30d, Female.30d), na.rm=TRUE))
    orgExpression <- orgExpression[,c(colnames(orgExpression)[regexpr("ale", colnames(orgExpression))<0])]
  }
} else if(organism=="Xenopus"){
  ### Gene information
  #Ensembl.Gene.ID, Associated.Gene.Name, X..GC.content
  #orgGenes <- read.delim(paste("~/Gene strucuture data"))
  orgGenes <- read.delim(paste("XenoGenesEnsV80.txt"),sep="\t") 
  #check import : colnames(orgGenes)
  
  ### Paralog information
  #Ensembl.Gene.ID, "Organism".Paralog.Ensembl.Gene.ID, Homology.Type
  #orgParalogs <- read.table(paste("~/Paralog information"), sep="\t", header=TRUE)
  orgParalogs <- read.table(paste("XenoParalogsEnsV80.txt"), sep="\t", header=TRUE)
  #check import : colnames(orgParalogs)
  
  ### PPI information (protein-protein interaction) 
  #sciName, kingdom, taxID, protein, locus, scorecutoff, connectivity, toppercentile
  #The number of direct neighbors of genes in protein-protein network
  #orgConnectivity <- read.table(paste("~/PPI information"), sep="\t", header=TRUE)
  orgConnectivity <- read.table(paste("connectivity.txt"), sep="\t", header=TRUE)
  orgConnectivity <- subset(orgConnectivity, regexpr("Xenopus tropicalis", sciName)>0) 
  #check if this is the good dataset : colnames(orgConnectivity)
  
  ### Protein gene connection
  #Ensembl.Gene.ID, Ensembl.Protein.ID
  # => need transcript ID for merging with connectivity
  #orgProtein <- read.table(paste("~/Protein information"), sep="\t", header=TRUE)
  orgProtein <- read.table(paste("XenoProteinEnsV80.txt"), sep="\t", header=TRUE)
  #check if this is the good dataset : colnames(orgProtein)
  
  if(phylAge){
    ### Phyletic ages of genes from OGEE
    #sciName, kingdom, locus, phyleticage, phyleticageTaxID, taxID, sciName.1, scorecutoff
    #orgPhyleticage <- read.table(paste("~/Phyletic age information"), sep="\t", header=TRUE)
    #orgPhyleticage <- read.table(paste("phyleticage.txt"), sep="\t", header=TRUE)
    #orgPhyleticage <- subset(orgPhyleticage, regexpr("Xenopus tropicalis", sciName)>0) #0 
    
    # Phyletic age from Ensembl
    orgPhyleticage <- read.table("ageXeno.txt", sep="\t")    
  }
  
  if(essentiality){
    ### Essentiality of genes
    #omim_pheno_accession, omim_gene_descriptio, omim_morbid_description, mgi_accession, ens_mouse_id, 
    #ens_human_id, orthology_type, human_omim_desc_essentiality, mouse_ontology_essentiality
    #orgEssentiality <- read.table(paste("~/Essentiality information"), sep="\t", header=TRUE)
    orgEssentiality <- read.table(paste("Ogee_essentiality_tmXenopus.txt"), sep="\t", header=TRUE)
    #only text mining available (and only 6 genes)
  }            
  
  ### Omega from Selectome
  #Ensembl.Gene.ID, Ensembl.Transcript.ID, Omega.0, Omega.2, LRT, P.value, Q.value, 
  #Positive.Selection, P.0, P.1, P2 
  #orgSelectome <- read.table(paste("~/Evolutionary rate information"), sep="\t", header=TRUE)
  orgSelectome <- read.table(paste("Xenopus_selectome_Euteleostomi.txt"), sep="\t", header=TRUE)  
  orgSelectome <- orgSelectome[,c("gene","omega0","lrt","selected","p0","p1")]
  colnames(orgSelectome) <- c("Ensembl.Gene.ID","Omega.0","LRT","Positive.Selection","P.0","P.1")
  
  
  ### Structure  
  #Ensembl.Gene.ID, Ensembl.Transcript.ID, CDS.Length, Exon.Rank.in.Transcript, Exon.Chr.Start..bp., 
  #Exon.Chr.End..bp.
  #Transcript ID needed ! If not transcript ID, I have many rows for 1 gene ID !
  #orgStructure <- read.table(paste("~/Gene structure information"), sep=",", header=TRUE)
  orgStructure <- read.table(paste("XenoStructureEnsV80.txt"), sep="\t", header=TRUE)
  
  ### Gene expression 
  #orgExpression <- read.table(paste("~/Expression information"), sep="\t", header=TRUE)
  orgExpression <- read.delim("xenoGSE37452.tsv") #  (RNA-seq from Bgee)
  orgExpression <- orgExpression[,c("Gene.ID","Stage.ID","Stage.name","RPKM")] 
  
} else if(organism=="Caeno"){
  ### Gene information
  #Ensembl.Gene.ID, Associated.Gene.Name, X..GC.content
  #orgGenes <- read.delim(paste("~/Gene strucuture data"))
  orgGenes <- read.delim(paste("CaenoGenes.txt"),sep="\t") 
  orgGenes <- orgGenes[,c("Gene.stable.ID", "Gene.name", "X..GC.content")]
  colnames(orgGenes) <- c("Ensembl.Gene.ID", "Associated.Gene.Name", "X..GC.content")
  
  ### Paralog information
  #Ensembl.Gene.ID, "Organism".Paralog.Ensembl.Gene.ID, Homology.Type
  #orgParalogs <- read.table(paste("~/Paralog information"), sep="\t", header=TRUE)
  orgParalogs <- read.table(paste("CaenoParalogs.txt"), sep="\t", header=TRUE)
  orgParalogs <- orgParalogs[,c("Gene.stable.ID", "Paralog.gene.stable.ID", "Homology.type")]
  colnames(orgParalogs) <- c("Ensembl.Gene.ID", "Paralog.Ensembl.Gene.ID","Homology.Type")
  
  ### PPI information (protein-protein interaction) 
  #sciName, kingdom, taxID, protein, locus, scorecutoff, connectivity, toppercentile
  #The number of direct neighbors of genes in protein-protein network
  #orgConnectivity <- read.table(paste("~/PPI information"), sep="\t", header=TRUE)
  orgConnectivity <- read.table(paste("connectivity.txt"), sep="\t", header=TRUE)
  orgConnectivity <- subset(orgConnectivity, regexpr("Caenorhabditis elegans", sciName)>0)
  
  ### Protein gene connection
  #Ensembl.Gene.ID, Ensembl.Protein.ID
  # => need transcript ID for merging with connectivity
  #orgProtein <- read.table(paste("~/Protein information"), sep="\t", header=TRUE)
  orgProtein <- read.table(paste("CaenoProtein.txt"), sep="\t", header=TRUE)
  orgProtein <- orgProtein[,c("Gene.stable.ID","Protein.stable.ID")]
  colnames(orgProtein) <- c("Ensembl.Gene.ID", "Ensembl.Protein.ID")
  
  if(phylAge){
    ### Phyletic ages of genes from OGEE
    #sciName, kingdom, locus, phyleticage, phyleticageTaxID, taxID, sciName.1, scorecutoff
    #orgPhyleticage <- read.table(paste("~/Phyletic age information"), sep="\t", header=TRUE)
    #orgPhyleticage <- read.table(paste("phyleticage.txt"), sep="\t", header=TRUE)
    #orgPhyleticage <- subset(orgPhyleticage, regexpr("Caenorhabditis elegans", sciName)>0) 
    
    # with the phyletic age retrieved from Ensembl
    orgPhyleticage <- read.table("ageCaeno.txt", sep="\t")
  }
  
  if(essentiality){
    ### Essentiality of genes
    #omim_pheno_accession, omim_gene_descriptio, omim_morbid_description, mgi_accession, ens_mouse_id, 
    #ens_human_id, orthology_type, human_omim_desc_essentiality, mouse_ontology_essentiality
    #only experiment taken (not text mining)
    #orgEssentiality <- read.table(paste("~/Essentiality information"), sep="\t", header=TRUE)
    orgEssentiality <- read.table(paste("Ogee_essentiality_expCaeno.txt"), sep="\t", header=TRUE) 
  }  
  
  #No Selectome data available for C. elegans ! -> calculate Omega with dN/dS
  #Omega.0 : purfiying selection, only if dN/d <1 else NA
  orgOrthologs <- read.table(paste("CaenoOrthologs.txt"), header=TRUE, sep="\t") 
  orgOrthologs <- orgOrthologs[,c("Gene.stable.ID", "dN", "dS", "Homology.type")]
  colnames(orgOrthologs) <- c("Ensembl.Gene.ID", "dN", "dS", "Homology.Type")
  #Omega is calculated as dN/dS if dS is > 0, otherwise Omega = 0
  orgOrthologs$Omega <- ifelse(orgOrthologs$dS>0,orgOrthologs$dN/orgOrthologs$dS,0)  
  orgOrthologs$Omega.0 <- ifelse(orgOrthologs$Omega<1, orgOrthologs$Omega, NA)
  
  ### Structure  
  #Ensembl.Gene.ID, Ensembl.Transcript.ID, CDS.Length, Exon.Rank.in.Transcript, Exon.Chr.Start..bp., 
  #Exon.Chr.End..bp.
  #Transcript ID needed ! If not transcript ID, I have many rows for 1 gene ID !!!
  #orgStructure <- read.table(paste("~/Gene structure information"), sep=",", header=TRUE)
  orgStructure <- read.table(paste("CaenoStructure.txt"), sep="\t", header=TRUE)
  orgStructure$CDS.Length <- orgStructure$CDS.end..within.cDNA. - orgStructure$CDS.start..within.cDNA.
  temp<-(orgStructure[orgStructure$CDS.Length<0,]); nrow(temp);temp<-na.omit(temp);nrow(temp)
  orgStructure <- orgStructure[,c("Gene.stable.ID","Transcript.stable.ID", "CDS.Length",
                                  "Exon.rank.in.transcript","Exon.region.start..bp.", "Exon.region.end..bp.")]
  colnames(orgStructure) <- c("Ensembl.Gene.ID", "Ensembl.Transcript.ID", "CDS.Length", 
                              "Exon.Rank.in.Transcript", "Exon.Chr.Start..bp.","Exon.Chr.End..bp.")
  ### Gene expression 
  #orgExpression <- read.table(paste("~/Expression information"), sep="\t", header=TRUE)
  if(!FPKM){ 
    #we had RNA-seq data from Bgee, but very few stages -> so in fact we only used 
    #the "FPKM" data from Li et al. (2014)
    orgExpression <- read.delim("caenoSRP000401.tsv") #  (RNA-seq from Bgee)    
    orgExpression <- orgExpression[,c("Gene.ID","Stage.ID","Stage.name","RPKM")] 
  } else if(FPKM){
    orgExpression <- read.delim("caeno_FPKMs.txt")
    colnames(orgExpression)[1] <- "Ensembl.Transcript.ID"
    devStageNames <- colnames(orgExpression)[-1]
    orgExpression[,-1] <- as.numeric(as.character(unlist((orgExpression[,-1]))))
    geneTransc <- read.delim("CaenoGenesTransc.txt") #Ens81, no filters, unique results only
    orgExpression <- merge(orgExpression, geneTransc, by="Ensembl.Transcript.ID")
    orgExpression <- orgExpression[,c("Ensembl.Gene.ID", devStageNames)]
  }
} 

#####
# Reshape the datasets--------------------------------------------------------------------------
####

### Data for gene expression
if(!twoDataSets & !FPKM){
    orgExpression <- orgExpression[,c("Gene.ID","Stage.ID","Stage.name",values)] 
    colnames(orgExpression) <- c("Ensembl.Gene.ID","Stage.ID","Stage.name",values)
}

if(organism!="Caeno"){
  cat("\nOverall ",nrow(orgSelectome)," transcripts with available evolutionary rate.","\nSummary :",sep="")
  summary(orgSelectome)
}

cat("\nOverall ",nrow(orgParalogs)," transcripts in paralogs data.","\nSummary :",sep="")
summary(orgParalogs)

### Take only within-species paralogs
orgParalogs <- orgParalogs[regexpr("within_species",orgParalogs$Homology.Type)>0,] 

### Calculate paralog numbers for each Ensembl.Gene.ID
sumParalogs <- count(orgParalogs,"Ensembl.Gene.ID")          # new df with 2 col ("Ensembl.Gene.ID","freq")
names(sumParalogs) <- c("Ensembl.Gene.ID","Paralogs.Number") #replace colname "freq" 
orgParalogs <- orgParalogs[,c("Ensembl.Gene.ID"),drop=FALSE]  
orgParalogs <- merge(orgParalogs, sumParalogs, by=c("Ensembl.Gene.ID"), all.x=TRUE, sort=FALSE)
orgParalogs <- orgParalogs[,c("Ensembl.Gene.ID","Paralogs.Number")]
orgParalogs <- unique(orgParalogs)
cat("\nOverall ",nrow(orgParalogs)," transcripts with within species paralogs.","\nSummary :",sep="")
summary(orgParalogs)

### Connectivity
orgConnectivity <- orgConnectivity[,c("locus", "protein","connectivity")]
colnames(orgConnectivity) <- c("Ensembl.Gene.ID", "Ensembl.Protein.ID","Connectivity")
orgConnectivity <- merge(orgConnectivity, orgProtein, by=c("Ensembl.Gene.ID","Ensembl.Protein.ID"), all.x=TRUE, sort=FALSE)
orgConnectivity <- na.omit(orgConnectivity)
orgConnectivity <- orgConnectivity[,c("Ensembl.Gene.ID","Connectivity")]
cat("\nSummary for connectivity data :")
summary(orgConnectivity)

### Phyletic age
if(phylAge){  
  colnames(orgPhyleticage) <- c("Ensembl.Gene.ID", "Taxon")
  orgPhyleticage <- merge(orgPhyleticage, ageTable, by="Taxon")
  orgPhyleticage$Taxon <- NULL
      #orgPhyleticage <- orgPhyleticage[, c("locus", "phyleticage")]  # (was for OGEE data)
      colnames(orgPhyleticage) <- c("Ensembl.Gene.ID", "Phyletic.Age")
      cat("\nSummary for phyletic age data :")
      summary(orgPhyleticage)
}

### Essentiality
if(essentiality){
  summary(orgEssentiality)
  if(organism == "Mus"){
    orgEssentiality <- orgEssentiality[,c("ens_mouse_id", "mouse_ontology_essentiality")]
    colnames(orgEssentiality) <- c("Ensembl.Gene.ID", "Essentiality")
    orgEssentiality$Essentiality <- ifelse(orgEssentiality$Essentiality == "yes", 1, 0)
  } else {
    orgEssentiality <- orgEssentiality[,c("locus","essential")]
    colnames(orgEssentiality) <- c("Ensembl.Gene.ID","Essentiality")
    orgEssentiality$Essentiality <- ifelse(orgEssentiality$Essentiality == "Y", 1, 0)
  }
  cat("\nSummary for essentiality data :")
  summary(orgEssentiality)
}

cat("\nOverall ",nrow(orgGenes)," transcripts", " in ",length(unique(orgGenes$Ensembl.Gene.ID))," genes for ", organismName," (",dataSource,").","\nSummary :",sep="")
summary(orgGenes)

#####*****Collect all files for gene properties in one*****#####
# Affymetrix Bgee data : only with Ensembl.Gene.Id matching (not Ensembl.Transcript.ID)
#total <- merge(orgGenes,orgOrthologs,by=c("Ensembl.Gene.ID"), all.x=TRUE, sort=FALSE)
total <- merge(orgGenes,orgParalogs,by=c("Ensembl.Gene.ID"), all.x=TRUE, sort=FALSE)
total$Paralogs.Number <- ifelse(is.na(total$Paralogs.Number), 0, total$Paralogs.Number)
if(phylAge){
  total <- merge(total,orgPhyleticage,by=c("Ensembl.Gene.ID"), all.x=TRUE, sort=FALSE)
}
if(essentiality){
  total <- merge(total,orgEssentiality,by=c("Ensembl.Gene.ID"), all.x=TRUE, sort=FALSE)
}
total <- merge(total,orgConnectivity, by=c("Ensembl.Gene.ID"), all.x=TRUE, sort=FALSE)
if(organism!="Caeno"){
  total <- merge(total,orgSelectome,by=c("Ensembl.Gene.ID"), all.x=TRUE, sort=FALSE)
}
if(organism == "Caeno"){
  total <- merge(total,orgOrthologs,by=c("Ensembl.Gene.ID"), all.x=TRUE, sort=FALSE)
}

# totalSave <- total

#####*****Genes structure data*****#####

cat("\nOverall ",nrow(orgStructure)," exons for ", organismName," (",dataSource,").","\nSummary :",sep="")
summary(orgStructure)
structure <- orgStructure

#####Compute values for gene structure
### Calculate the length of each exon
structure$Exon.Length <- ifelse(structure$Exon.Chr.End..bp.>0,structure$Exon.Chr.End..bp. - structure$Exon.Chr.Start..bp.+1,structure$Exon.Chr.End..bp.)

### Calculations for Introns
#For each transcript -> sum all exons' length 
exonLength <- aggregate(Exon.Length ~ Ensembl.Transcript.ID, FUN="sum", data=structure) #1 T.ID=sum Exons
names(exonLength) <- c("Ensembl.Transcript.ID","Exon.Total.Length")

#Calculate the number of all exons in the transcript (=find the highest Exon.rank.)
exonNumber <- aggregate(Exon.Rank.in.Transcript ~ Ensembl.Transcript.ID, FUN="max", data=structure)
names(exonNumber) <- c("Ensembl.Transcript.ID","Exon.Number")

#Find the start of the first exon in transcript (=find the smallest Exon.Chr.Start..bp.)
# = find the start of the transcript
transcriptStart <- aggregate(Exon.Chr.Start..bp. ~ Ensembl.Transcript.ID, FUN="min", data=structure)
names(transcriptStart) <- c("Ensembl.Transcript.ID","Transcript.Start")

#Find the end of the last exon in transcript (=find the highest Exon.Chr.End..bp.)
# = find the end of the transcript
transcriptEnd <- aggregate(Exon.Chr.End..bp. ~ Ensembl.Transcript.ID, FUN="max", data=structure)
names(transcriptEnd) <- c("Ensembl.Transcript.ID","Transcript.End")

#Calculate the longest transcript for each gene
maxCDS <- aggregate(CDS.Length ~ Ensembl.Gene.ID, FUN="max", data=structure)
names(maxCDS) <- c("Ensembl.Gene.ID","Max.CDS.Length")

### Put all the calculated data to one table
structure <- merge(structure,exonLength,by=c("Ensembl.Transcript.ID"), all.x=TRUE, sort=FALSE) 
structure <- merge(structure,exonNumber,by=c("Ensembl.Transcript.ID"), all.x=TRUE, sort=FALSE) 
structure <- merge(structure,transcriptStart,by=c("Ensembl.Transcript.ID"), all.x=TRUE, sort=FALSE) 
structure <- merge(structure,transcriptEnd,by=c("Ensembl.Transcript.ID"), all.x=TRUE, sort=FALSE) 
structure <- merge(structure,maxCDS,by=c("Ensembl.Gene.ID"), all.x=TRUE, sort=FALSE)

#Calculate the length of the transcript
structure$Transcript.Length <- ifelse(structure$Transcript.End>0,structure$Transcript.End - structure$Transcript.Start+1,structure$Transcript.End)

#Calculate the number of introns
structure$Intron.Number <- ifelse(structure$Exon.Number>0,structure$Exon.Number-1,structure$Exon.Number)

#Calculate the length of the Introns
structure$Intron.Length <- ifelse(structure$Transcript.Length>0,(structure$Transcript.Length - structure$Exon.Total.Length),structure$Transcript.Length)

#Intron Length is the mean length of introns in transcript
cat("\nIntron Length is the mean length of introns in the transcript.")
structure$Intron.Length <- ifelse(structure$Intron.Number>0,structure$Intron.Length/structure$Intron.Number,structure$Intron.Number)

#Label if transcript is longest for this gene
structure$Max.CDS.Length <- ifelse(structure$CDS.Length==structure$Max.CDS.Length,TRUE,FALSE)

#Choose only columns and rows that are needed
structure <- structure[,c('Ensembl.Gene.ID','Ensembl.Transcript.ID','CDS.Length','Intron.Length','Intron.Number',"Max.CDS.Length")]
#Many rows are the same, because before there was a row for each Exon, so just delete duplicates
structure <- unique(structure) 

# The structure data is the only table for which we have transcript ID. 
# So before merging with the other files, we need to select one transcript for each gene.

### Choice of the transcript : select one unique Ensemb.Gene.ID for the structure table
# 1st criterion : Max.CDS.length
# 2nd criterion : longest introns
# (3d criterion : the more introns)
# 3d criterion : random

## 1. Choose the transcripts with the maximal CDS length
structure1 <-subset(structure, Max.CDS.Length==TRUE) 
freq<-count(structure1,"Ensembl.Gene.ID")

oneTrans <- freq$Ensembl.Gene.ID[which(freq$freq==1),drop=FALSE]  
names(oneTrans)<-"Ensembl.Gene.ID"
structure1a <- subset(structure1, Ensembl.Gene.ID %in% oneTrans)
nrow(structure1a)==length(unique(structure1a$Ensembl.Gene.ID)) #=> TRUE
#df with the the genes that have only 1 transcript that has Max.CDS.length
#1a ok for final merging

manyTrans <- freq$Ensembl.Gene.ID[which(freq$freq>1),drop=FALSE] 
names(manyTrans)<-"Ensembl.Gene.ID"
structure1b <- subset(structure1, Ensembl.Gene.ID %in% manyTrans)
nrow(structure1b)==length(unique(structure1b$Ensembl.Gene.ID)) 

## 2. Maximal intron length 
#Maximal intron length for the rest 
maxIntron <- aggregate(Intron.Length ~ Ensembl.Gene.ID, FUN="max", data=structure1b) 
names(maxIntron) <- c("Ensembl.Gene.ID","Max.Intron") 
structure2 <- merge(structure1b, maxIntron, by="Ensembl.Gene.ID", all.x=TRUE)
structure2 <- structure2[structure2$Intron.Length==structure2$Max.Intron,] 

freq<-count(structure2,"Ensembl.Gene.ID")

oneTrans <- freq$Ensembl.Gene.ID[which(freq$freq==1),drop=FALSE]  
names(oneTrans)<-"Ensembl.Gene.ID"
structure2a <- subset(structure2, Ensembl.Gene.ID %in% oneTrans) 
structure2a <- structure2a[,-length(colnames(structure2a))] #remove Max.Intron col
nrow(structure2a)==length(unique(structure2a$Ensembl.Gene.ID)) #=> TRUE
#df with the the genes that have only 1 transcript that has Max.Intron
#2a ok for final merging

manyTrans <- freq$Ensembl.Gene.ID[which(freq$freq>1),drop=FALSE] 
names(manyTrans)<-"Ensembl.Gene.ID"
structure2b <- subset(structure2, Ensembl.Gene.ID %in% manyTrans)  
nrow(structure2b)==length(unique(structure2b$Ensembl.Gene.ID)) #=> FALSE

## 3. The one with the more introns ?
temp<-structure2b
length(!duplicated(temp)) 
temp<-temp[,-2] 
nrow(unique(temp))
# by removing the transcript id and take only the unique rows -> as many as rows as unique Gene.id
# => they differ only by the transcript id

## 4. Random one for the rest. 
structure3 <- structure2b
temp <- split(1:nrow(structure3),structure3$Ensembl.Gene.ID) 
temp2 <- sapply(temp,function(x){x <- x[1]}) 
structure3 <- structure3[temp2,] 
structure3 <- structure3[!is.na(structure3$Ensembl.Gene.ID),] 
structure3 <- structure3[,-(length(colnames(structure3)))]
nrow(structure3)==length(unique(structure3$Ensembl.Gene.ID)) #=> TRUE
# 3 ok for final merging

#Merge the df
endStruct <- rbind(structure1a, structure2a, structure3) 
nrow(endStruct)==length(unique(endStruct$Ensembl.Gene.ID)) #=> TRUE => all unique !


#####*****Merge tables for gene properties and gene structure*****#####
#Merge the two tables (with gene properties and gene structure)
# tempT <- total
total <- merge(total,endStruct,by=c("Ensembl.Gene.ID"), all.x=TRUE)


#####*****Choice of used genes*****#####
#If many expression values for one gene at one developmental stage -> mean values
cat("\n If many replicates for one gene at the same developmental stage, then the mean of expression is chosen.")

if(!twoDataSets & !FPKM){
  if(expDataSource == "Bgee"){
    orgExpression<-aggregate(Log.of.normalized.signal.intensity~Ensembl.Gene.ID+Stage.name,data=orgExpression,na.action=na.omit,FUN=mean)  
  }else if(expDataSource == "BgeeRNA"){
    orgExpression<-aggregate(RPKM~Ensembl.Gene.ID+Stage.name,data=orgExpression,na.action=na.omit,FUN=mean)  
  }
}

#For the RNA-seq data: RPKM values that are <1 => considered as not expressed 
if(!twoDataSets){
  if(!FPKM){
    if(expDataSource == "BgeeRNA"){
      orgExpression$RPKM[orgExpression$RPKM<1] <- 0
      min(orgExpression$RPKM[orgExpression$RPKM>0]) #check : must be > 1.000113 #Caeno : 1.000028
    }
  } else if(FPKM){
    orgExpression[,-1] <- apply(orgExpression[,-1],c(1,2),function(x){ifelse(x<1,0,x)})
  }
}

#####*****Change the shape of the data table*****#####
#From:
# Ensembl.Gene.ID  Stage.name  Log.of.normalized.signal.intensity
# *gene 1*          *stage1*      *valueA*
# *gene 1*          *stage2*      *valueB*
# *gene 1*          *stage3*      *valueC*
# *gene 2*          *stage1*      *valueD*
# *gene 2*          *stage2*      *valueE*
# *gene 2*          *stage3*      *valueF*

# ... to:
# Ensemble.Gene.ID    Stage1      Stage2          Stage3
# *gene1*           *valueA*      *valueB*      *valueC*      
# *gene2*           *valueD*      *valueE*      *valueF*

# orgExpT <- orgExpression
#use a pre-existing function (reshape2::dcast) (!efficient)
if(!twoDataSets & !FPKM){
  orgExpression <- dcast(orgExpression, Ensembl.Gene.ID~Stage.name, value.var=values)
}

# colnames(orgExpression)

#For mouse and drosophila, import dataset after Combat Normalization
if(twoDataSets){
  orgExpression <- importAfterComBat(organism)
}
nDevStages<-(ncol(orgExpression)-1)

### Order of the columns should be arranged according to the timing of developmental stages
#(i.e. from zygote to adult)
if(organism=="Mus"){
  orgExpression<-orgExpression[,c("Ensembl.Gene.ID","zygote stage", "Theiler stage 02 (mouse)", "Theiler stage 03 (mouse)",
                                  "Theiler stage 04 (mouse)", "Theiler stage 05 (mouse)", "Theiler stage 06 (mouse)",
                                  "neurula stage", "Theiler stage 13 (mouse)", "Theiler stage 15 (mouse)", 
                                  "Theiler stage 17 (mouse)", "Theiler stage 20 (mouse)", "Theiler stage 22 (mouse)",
                                  "Theiler stage 24 (mouse)","Theiler stage 26 (mouse)","post-juvenile adult stage")]
  
} else if(organism=="Danio"){
  orgExpression <- orgExpression[,c("Ensembl.Gene.ID","zygote stage","Gastrula:Shield (Danio)", "Gastrula:75%-epiboly (Danio)", 
                                    "Gastrula:90%-epiboly (Danio)", "Gastrula:Bud (Danio)", "Segmentation:5-9 somites (Danio)",
                                    "Segmentation:14-19 somites (Danio)", "Pharyngula:Prim-5 (Danio)", "Pharyngula:Prim-15 (Danio)",
                                    "Hatching:Long-pec (Danio)", "Larval:Day 4 (Danio)", "Larval:Day 5 (Danio)",
                                    "Larval:Days 14-20 (Danio)","Juvenile:Days 30-44 (Danio)")] 
} else if(organism=="Drosophila" && !FPKM){
  orgExpression <- orgExpression[,c("Ensembl.Gene.ID","embryonic stage 2 (Drosophila)", "embryonic stage 4 (Drosophila)",
                                    "embryonic stage 5 (Drosophila)", "early extended germ band stage (Drosophila)",
                                    "embryonic stage 10 (Drosophila)", "late extended germ band stage (Drosophila)", 
                                    "embryonic stage 17 (Drosophila)", "first instar larval stage (Drosophila)", 
                                    "second instar larval stage (Drosophila)")]                                 
} else if(organism=="Xenopus"){
  orgExpression<-orgExpression[,c("Ensembl.Gene.ID","2 cell stage", "4 cell stage", "8 cell stage",
                                  "NF stage 5 (16-cell) (Xenopus)", "NF stage 6 (32-cell) (Xenopus)",
                                  "NF stage 8 (Xenopus)", "NF stage 9 (Xenopus)", "NF stage 10 (Xenopus)",
                                  "gastrula stage",
                                  "NF stage 15 (Xenopus)", "NF stage 19 (Xenopus)", "NF stage 28 (Xenopus)",
                                  "NF stage 33 and 34 (Xenopus)", "NF stage 40 (Xenopus)")]
} else if (organism=="Caeno" && !FPKM){
  orgExpression<-orgExpression[,c("Ensembl.Gene.ID",
                                  "blastula stage", "gastrula stage", "embryo stage",
                                  "nematode larval stage", "fully formed stage")]
} 

#orgExpressionSave <- orgExpression

# ---------------------------------------------------------------------------------------
# EXPRESSION DISTRIBUTION
# ---------------------------------------------------------------------------------------

#####*****Change the names of the stages for the output*****#####
devStageNames<-colnames(orgExpression)[-1] #remove gene id
devStagePrintNames<-devStageNames
  ##PrintNames : "Theiler stage 02 (mouse)"  -> "Theiler stage 02" :
if(organism=="Mus"){
  for (i in 1:length(devStagePrintNames)){
    if(regexpr("mouse",devStagePrintNames[i])>0){
      devStagePrintNames[i]<-substring(devStagePrintNames[i],1,nchar(devStagePrintNames[i])-8)
    }
  }
} else if(organism=="Danio"){
  for (i in 1:length(devStagePrintNames)){
    if(regexpr("Danio",devStagePrintNames[i])>0){
      devStagePrintNames[i]<-substring(devStagePrintNames[i],1,nchar(devStagePrintNames[i])-8)
    }
  }
} else if(organism=="Drosophila" && !FPKM){
  for (i in 1:length(devStagePrintNames)){
    if(regexpr("Drosophila",devStagePrintNames[i])>0){
      devStagePrintNames[i]<-substring(devStagePrintNames[i],1,nchar(devStagePrintNames[i])-13)
    }
  }
} else if(organism=="Xenopus"){
  for (i in 1:length(devStagePrintNames)){
    if(regexpr("Xenopus",devStagePrintNames[i])>0){
      devStagePrintNames[i]<-substring(devStagePrintNames[i],1,nchar(devStagePrintNames[i])-10)
    }
  }
} 

#####*****Draw expression distribution*****#####
#Use the drawExpressionDistribution() function defined in"scriptFunctions.R"

devStagePrintNamesSave <- devStagePrintNames

### With "raw" data (Bgee are already log-transf)
# NB: named"RPKM" but for the microarray data = log.of.normalized.signal.intensity from Bgee
orgRPKM<-na.omit(orgExpression)
drawExpressionDistribution(values)
dev.copy2pdf(file=paste(folderAnalysis, organism,"DevStagesOriginalExpressionLogIntensities", expDataSource,".pdf", sep=""),onefile=TRUE)
dev.off()

if(expDataSource == "Bgee" && !FPKM){
  #No log2-transformation: Bgee expression values are already log-transf. and > 0
  # -> juste perform quantile normalization (without log-transf.)
  #Quantile normalization : 
  orgRPKM <- na.omit(orgExpression)
  tempRPKM <- orgRPKM[,devStageNames]
  tempRPKM_m <- as.matrix(tempRPKM)
  tempRPKM <- normalize.quantiles(tempRPKM_m)
  orgRPKM[,devStageNames] <- data.frame(tempRPKM)
  drawExpressionDistribution(paste("QN", values))
  dev.copy2pdf(file=paste(folderAnalysis, organism,"OriginalExpressionRPKM", expDataSource,"ExtendedQN.pdf", sep=""),onefile=TRUE)
  dev.off()

} else if(expDataSource == "BgeeRNA" || FPKM){
  #With log2-transformation :  
  orgRPKM <- na.omit(orgExpression)
  x2 <- orgRPKM[,devStageNames]  #do not take the Ensembl.Gene.ID column 
  # x2 <- x2 * expNorm #except 0, min is 4.1 -> >1 => not useful because <1 -> 0 before
  x2 <- log2(x2)
  x2[x2==-Inf] <- log2(1.0001)
  orgRPKM[,devStageNames] <- data.frame(x2)   
  #replace the modified columns in the original dataframe -> Ensembl.Gene.ID column "back" in the dataframe  
  drawExpressionDistribution("log2-RPKM")
  dev.copy2pdf(file=paste(folderAnalysis, organism,"DevStagesOriginalExpressionlogRPKM", expDataSource,".pdf", sep=""),onefile=TRUE)
  #dev.off()
  
  #With quantile normalization :
  orgRPKM <- na.omit(orgExpression)
  tempRPKM <- orgRPKM[,devStageNames]
  tempRPKM <- log2(tempRPKM)
  tempRPKM_m <- as.matrix(tempRPKM)
  tempRPKM <- normalize.quantiles(tempRPKM_m)
  tempRPKM[tempRPKM==-Inf] <- log2(1.0001)
  orgRPKM[,devStageNames] <- data.frame(tempRPKM)
  drawExpressionDistribution("Quantile normalized signal intensity")
  dev.copy2pdf(file=paste(folderAnalysis, organism,"OriginalExpressionRPKM", expDataSource,"ExtendedQN.pdf", sep=""),onefile=TRUE)
  #dev.off()  
}

orgExpression<-orgRPKM  #continue with the normalized data

#####*****Mean, median and maximal expression*****#####
### Mean
cat("\n Mean is calculated taking in account tissues with 0 expression. 2+0+4=2",sep="")

# Define own functions: if juste built-in mean function -> NaN but we want NA
fmean <- function(x){
  #x <- subset(x,x>0)
  if(!all(is.na(x))){
    res <- mean(x, na.rm=TRUE)
  } else {
    res <- NA
  }
  return(res)
}
orgExpression$Mean.Expression <- apply(orgExpression[,devStageNames],1,fmean)

### Median
cat("\n Median is calculated taking in account tissues with 0 expression. 2+0+4=2",sep="")

fmedian <- function(x){
  #x <- subset(x,x>0)
  if(!all(is.na(x))){
    res <- median(x, na.rm=TRUE)
    #*a logical value indicating whether NA values should be stripped before the computation proceeds
  }else {
    res <- NA
  }
  return(res)
}

orgExpression$Median.Expression <- apply(orgExpression[,devStageNames],1,fmedian)

### Maximal
# Function to calculate maximal expression over all stages
fmax <- function(x){
  if(!all(is.na(x))){
  res <- max(x, na.rm=TRUE)
  } else {
    res <- NA
  }
  return(res)
}
orgExpression$Max.Expression <- apply(orgExpression[,devStageNames],1,fmax)

#####*****Stage specificity*****#####
### Calculated in the same way as tissue specificity ?????
#Function to calculate Tau: x^=xi/max(xi), tau=sum(1-x^)/(n-1) = sum(1-xi/max(xi))/n-1
ftau <- function(x){
  if(!all(is.na(x))){
    x <- (1-(x/x[length(x)]))  #1-x^  #x[length(x)] => Max.Expression column == max(xi)
    res <- sum(x, na.rm=TRUE)  #sum(1-x^)
    res <- res/(length(x)-1)   #sum(1-x^)/(n-1)
  } else {
    res <- NA
  }
  return(res)
}
orgExpression$Tau <- apply(orgExpression[,c(devStageNames,paste("Max.Expression",sep=""))],1,ftau)

cat("\n Expression data are avalable for: ", nrow(orgExpression)," genes.",sep="")
cat("\n Summary of expression data after normalisation and calculating tau: ")
summary(orgExpression)
# orgExpressionSave2 <- orgExpression

### Change the names of the stages for the output
#e.g.: replace blank & - by .or append ".stage" (needed for regexpr)
# devStagePrintNamesSave2 <- devStagePrintNames
if(organism=="Mus"){
  devStageOutNames<-c()
  for(i in devStagePrintNames){
    x<-gsub(" ",".",i)
    x<-gsub("-",".",x)
    
    if(regexpr("stage",x)<0){
      x<-paste(x,".stage",sep="")
    }
    # if there are blank spaces and -, we may have added many "." -> remove repeated "." :
    x<-gsub('(\\.)\\1+', '\\1', x)
    devStageOutNames<-append(devStageOutNames,x)
  }
  for(i in 1:length(devStageOutNames)){
    if(regexpr("stage",devStageOutNames[i])<0){
      devStageOutNames[i]<-paste(devStageOutNames[i],".stage",sep="")
    }
  }
} else if(organism=="Danio"){
  devStageOutNames<-devStagePrintNames
  for(i in 1:length(devStageOutNames)){
    if(regexpr("stage",devStageOutNames[i])<0){
      devStageOutNames[i] <- paste0("stage",i,".",devStageOutNames[i])      
    } 
     devStageOutNames[i] <-gsub(":",".",devStageOutNames[i])
     devStageOutNames[i] <-gsub("-",".",devStageOutNames[i])
    devStageOutNames[i]<-gsub('%', '.', devStageOutNames[i])
    devStageOutNames[i]<-gsub(' ', '.', devStageOutNames[i])
  }
} else if(organism=="Drosophila" && !FPKM){
  devStageOutNames<-devStagePrintNames
  for(i in 1:length(devStageOutNames)){
    if(regexpr("stage",devStageOutNames[i])<0){
      devStageOutNames[i]<-paste(devStageOutNames[i],".stage",sep="") #normally all have already "stage"
    } 
    devStageOutNames[i]<-gsub(' ', '.', devStageOutNames[i])
  }
} else if(organism=="Xenopus"){
  devStageOutNames <- vector()
  for(i in devStagePrintNames[regexpr("^[0-9]",devStagePrintNames)>0]){ #sentences beginning with numbers
    i <- paste(word(i,3), word(i,2), word(i,1))
    devStageOutNames <- append(devStageOutNames, i)
  }
  temp <- devStagePrintNames[regexpr("^[0-9]",devStagePrintNames)<0]  
  temp<-gsub("\\(", ".", temp)
  temp<-gsub("\\)", ".", temp)
  temp<-gsub("-", ".", temp)
  
  devStageOutNames <- append(devStageOutNames, temp)
  devStageOutNames<-gsub(" ", ".", devStageOutNames)

} else if(organism=="Caeno" && !FPKM){
  devStageOutNames <- devStagePrintNames
  devStageOutNames <- gsub(" ", ".", devStageOutNames)

} else if(FPKM){  #later we need to have "stage" in the name of the developmental stage
  devStageOutNames <- devStagePrintNames
  devStageOutNames <- paste0(devStageOutNames, ".stage")
}
colnames(orgExpression)[(1+1):(1+length(devStageNames))]<-devStageOutNames #the first col is gene ID -> not change it

write.table(orgExpression, file=paste(folderAnalysis, organism, "Expression", expDataSource,".txt",sep=""),row.names = FALSE, col.names=TRUE, quote = FALSE)
summary(orgExpression)

write.table(total, file=paste(folderAnalysis, organism, "TotalSum", expDataSource,"Parameters", ".txt",sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)

#####*****Add calculated values to the main table and save files*****#####
#Add calculated expression values to the main table
# tempT <- total
# orgExpressionSave3 <- total
total <- merge(total,orgExpression,by=c("Ensembl.Gene.ID"), all.x=TRUE, incomparables = NA, sort=FALSE)

### Choose only useful columns
if(organism=="Caeno"){  #special case: no Selectome data (-> save without P.1, LRT and Positive.Selection)
  if(essentiality){
    if(phylAge){
      totalDevStages <- total[,c("Ensembl.Gene.ID","CDS.Length", "Intron.Length", "Intron.Number", 
                                 "Omega.0","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", 
                                 "X..GC.content", "Paralogs.Number","Connectivity", "Phyletic.Age", "Essentiality",
                                 devStageOutNames)]
    } else if(!phylAge){
      totalDevStages <- total[,c("Ensembl.Gene.ID", "CDS.Length", "Intron.Length", "Intron.Number", 
                                 "Omega.0","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", 
                                 "X..GC.content", "Paralogs.Number","Connectivity", "Essentiality",
                                 devStageOutNames)]
    }
  }else if(!essentiality){
    if(phylAge){
      totalDevStages <- total[,c("Ensembl.Gene.ID","CDS.Length", "Intron.Length", "Intron.Number", 
                                 "Omega.0","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", 
                                 "X..GC.content", "Paralogs.Number","Connectivity", "Phyletic.Age",
                                 devStageOutNames)]
    } else if(!phylAge){
      totalDevStages <- total[,c("Ensembl.Gene.ID", "CDS.Length", "Intron.Length", "Intron.Number", 
                                 "Omega.0","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", 
                                 "X..GC.content", "Paralogs.Number","Connectivity",
                                 devStageOutNames)]
    }
  }
} else{
  if(essentiality){
    if(phylAge){
        totalDevStages <- total[,c("Ensembl.Gene.ID","CDS.Length", "Intron.Length", "Intron.Number", 
                                    "Omega.0", "LRT", "P.1", "Positive.Selection","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", 
                                   "X..GC.content", "Paralogs.Number","Connectivity", "Phyletic.Age", "Essentiality",
                                 devStageOutNames)]
    } else if(!phylAge){
      totalDevStages <- total[,c("Ensembl.Gene.ID", "CDS.Length", "Intron.Length", "Intron.Number", 
                                 "Omega.0", "LRT", "P.1", "Positive.Selection","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", 
                                 "X..GC.content", "Paralogs.Number","Connectivity", "Essentiality",
                                 devStageOutNames)]
    }
  } else if(!essentiality){
    if(phylAge){
      totalDevStages <- total[,c("Ensembl.Gene.ID","CDS.Length", "Intron.Length", "Intron.Number", 
                                 "Omega.0", "LRT", "P.1", "Positive.Selection","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", 
                                 "X..GC.content", "Paralogs.Number","Connectivity", "Phyletic.Age",
                                 devStageOutNames)]
    } else if(!phylAge){
      totalDevStages <- total[,c("Ensembl.Gene.ID", "CDS.Length", "Intron.Length", "Intron.Number", 
                                 "Omega.0", "LRT", "P.1", "Positive.Selection","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", 
                                 "X..GC.content", "Paralogs.Number","Connectivity",
                                 devStageOutNames)]
    }
  }
}

### Choose only useful colums
# totalSave <- total
if(organism=="Caeno"){  #we don't have the values from Selectome
  if(essentiality){
    if(phylAge){
      total <- total[,c("Ensembl.Gene.ID", "X..GC.content", "Omega.0", "CDS.Length", 
                        "Intron.Length", "Intron.Number","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", "Paralogs.Number", 
                        "Connectivity", "Phyletic.Age", "Essentiality")]
    } else if(!phylAge){
      total <- total[,c("Ensembl.Gene.ID", "X..GC.content", "Omega.0", 
                        "CDS.Length", "Intron.Length", "Intron.Number","Max.Expression", "Mean.Expression", 
                        "Median.Expression", "Tau", "Paralogs.Number","Connectivity", "Essentiality")]
    }
  } else if(!essentiality){
    if(phylAge){
      total <- total[,c("Ensembl.Gene.ID", "X..GC.content", "Omega.0", "CDS.Length", 
                        "Intron.Length", "Intron.Number","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", "Paralogs.Number", 
                        "Connectivity", "Phyletic.Age")]
    } else if(!phylAge){
      total <- total[,c("Ensembl.Gene.ID", "X..GC.content", "Omega.0", 
                        "CDS.Length", "Intron.Length", "Intron.Number","Max.Expression", "Mean.Expression", 
                        "Median.Expression", "Tau", "Paralogs.Number","Connectivity")]
    }
  }  
  
} else{
  if(essentiality){
    if(phylAge){
      total <- total[,c("Ensembl.Gene.ID", "X..GC.content", "Omega.0", "LRT", "P.1", "Positive.Selection", "CDS.Length", 
                        "Intron.Length", "Intron.Number","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", "Paralogs.Number", 
                        "Connectivity", "Phyletic.Age", "Essentiality")]
    } else if(!phylAge){
      total <- total[,c("Ensembl.Gene.ID", "X..GC.content", "Omega.0", "LRT", "P.1", "Positive.Selection", 
                      "CDS.Length", "Intron.Length", "Intron.Number","Max.Expression", "Mean.Expression", 
                      "Median.Expression", "Tau", "Paralogs.Number","Connectivity", "Essentiality")]
    }
  } else if(!essentiality){
    if(phylAge){
      total <- total[,c("Ensembl.Gene.ID", "X..GC.content", "Omega.0", "LRT", "P.1", "Positive.Selection", "CDS.Length", 
                        "Intron.Length", "Intron.Number","Max.Expression", "Mean.Expression", "Median.Expression", "Tau", "Paralogs.Number", 
                        "Connectivity", "Phyletic.Age")]
    } else if(!phylAge){
      total <- total[,c("Ensembl.Gene.ID", "X..GC.content", "Omega.0", "LRT", "P.1", "Positive.Selection", 
                        "CDS.Length", "Intron.Length", "Intron.Number","Max.Expression", "Mean.Expression", 
                        "Median.Expression", "Tau", "Paralogs.Number","Connectivity")]
    }
  }
}

cat("\nOverall ",nrow(total)," genes.","\nSummary :",sep="")
summary(total)
cat("\n Overall ",nrow(totalDevStages)," genes for each stage."," Summary:",sep="")
summary(totalDevStages)

### Save the results
write.table(total,file=paste(folderAnalysis, organism, "Table", expDataSource, ".txt",sep=""),row.names = FALSE, col.names=TRUE, quote = FALSE)
write.table(totalDevStages,file=paste(folderAnalysis, organism, "TableStages", expDataSource,".txt",sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)

# ---------------------------------------------------------------------------------------
# COMPUTE CORRELATIONS : Cytoscape
# ---------------------------------------------------------------------------------------
#####*****Make partial correlation with glm model for Cytoscape****#####
### Cytoscape : correlation only with the parameters, not with the expresssion at each developmental stages
#Load the data
data <- read.table(paste(folderAnalysis,organism, "Table", expDataSource, ".txt",sep=""), header=TRUE)
#partial <- FALSE #still defined if run from the beginning
#corMethod <- "spearman"   #"pearson" - still defined if run from the beginning
cat("\nSummary of the data in the first step : ", sep="")
summary(data)

### Select only needed columns
if(organism=="Caeno"){      #we don't have the values from Selectome
  if(phylAge){
    data <- data[,c("Omega.0","CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                    "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number", "Phyletic.Age")]
    colnames(data) <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                        "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number", "Phyletic.Age")
    parametersNames <- c("Omega", "CDS length", "Intron \n length", "Intron \n number", "Median \n expression", 
                         "Maximal \n expression", "Tau", "%GC content", "Paralogs \n number", "Phyletic age")
  } else if(!phylAge){
    data <- data[,c("Omega.0","CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                    "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number")]
    colnames(data) <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                        "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number")
    parametersNames <- c("Omega", "CDS length", "Intron \n length", "Intron \n number", "Median \n expression", 
                         "Maximal \n expression", "Tau", "%GC content", "Paralogs \n number")
  }
} else {
  if(phylAge){
    data <- data[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                    "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number", "Phyletic.Age")]
    colnames(data) <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                        "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number", "Phyletic.Age")
    parametersNames <- c("Omega", "LRT", "P 1", "CDS length", "Intron \n length", "Intron \n number", "Median \n expression", 
                         "Maximal \n expression", "Tau", "%GC content", "Paralogs \n number", "Phyletic age")
  } else if(!phylAge){
    data <- data[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                    "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number")]
    colnames(data) <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                        "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number")
    parametersNames <- c("Omega", "LRT", "P 1", "CDS length", "Intron \n length", "Intron \n number", "Median \n expression", 
                         "Maximal \n expression", "Tau", "%GC content", "Paralogs \n number")
  }
}

cat("\nTo calculate correlation ", corMethod, " correlation was used.",sep="")

if (partial){
  cat("\nPartial correlation was performed.",sep="")
  part <- "Partial"
} else {
  cat("\nNormal correlation was performed.",sep="")
  part <- "Normal"
}

#Delete all genes with NA and not known parameters
data <- na.omit(data)
data <- data[data$Max.Expression>0.00015,]  # genes that are not expressed in any stages
cat("\nAll the genes with unknown parameters are removed from the analysis. ", nrow(data), " are left for the analysis.", "\nSummary: ", sep="")
summary(data)
geneData <- data   

#Normalization of the data (function in "scriptFunctions.R")
geneData <-normalizeData(geneData)

cat("\nThe data are normalized. Log2 of all paremeters is taken, exept GC content.", sep="")
cat("Summary of the data after normalization :",sep="")
summary(geneData)
cat("\nGraphical representation of the data is saved in the file \"Parameters\".", sep="")
# geneDataSave <- geneData

### Draw histograms
#x11(height=12, width=18)
#x11()
pdf(paste(folderAnalysis, organism,"Parameters", expDataSource,".pdf", sep=""))
par(mfrow=c(6,4),cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
histogram(data$Omega[which(data$Omega < quantile(data$Omega,0.95))], main=paste("Omega"),xlab="Omega",ylab="Percentage of total",breaks=30) 
histogram(geneData$Omega, main=paste("log2(Omega)"),xlab="log2(Omega)", col=my.col[3],ylab="Percentage of total",breaks=30) 
histogram(data$CDS.Length[which(data$CDS.Length < quantile(data$CDS.Length,0.95))], main=paste("CDS length"), xlab="CDS length",ylab="Percentage of total",breaks=30) 
histogram(geneData$CDS.Length, main=paste("log2(CDS length)"), xlab="log2(CDS length)", col=my.col[3],ylab="Percentage of total",breaks=30) 
histogram(data$Intron.Number[which(data$Intron.Number<quantile(data$Intron.Number,0.95))], main=paste("Intron number"), xlab="Intron number",ylab="Percentage of total",breaks=30) 
histogram(geneData$Intron.Number, main=paste("log2(Intron number)"), xlab="log2(Intron number)", col=my.col[3],ylab="Percentage of total",breaks=30) 
histogram(data$Intron.Length[which(data$Intron.Length<quantile(data$Intron.Length,0.95))], main=paste("Intron length"), xlab="Intron length",ylab="Percentage of total",breaks=30) 
histogram(geneData$Intron.Length, main=paste("log2(Intron length)"), xlab="log2(Intron length)", col=my.col[3],ylab="Percentage of total",breaks=30) 
histogram(data$Tau, main=paste("Tau"), xlab="Tau",ylab="Percentage of total",breaks=30) 
histogram(geneData$Tau, main=paste("log2(Tau)"), xlab="log2(Tau)", col=my.col[3],ylab="Percentage of total",breaks=30) 
histogram(data$Paralogs.Number[which(data$Paralogs.Number<quantile(data$Paralogs.Number,0.95))], main=paste("Paralogs Number"), xlab="Paralogs Number",ylab="Percentage of total",breaks=30) 
histogram(geneData$Paralogs.Number, main=paste("log2(Paralogs Number)"), xlab="log2(Paralogs Number)", col=my.col[3],ylab="Percentage of total",breaks=30) 
if(organism!="Caeno"){ #no Selectome data available for Caeno
  histogram(data$P.1[which(data$P.1<quantile(data$P.1,0.95))], main=paste("P.1"), xlab="P.1",ylab="Percentage of total",breaks=30) 
  histogram(geneData$P.1, main=paste("(P.1)^1/4"), xlab="(P.1)^1/4", col=my.col[3],ylab="Percentage of total",breaks=30) 
  histogram(data$LRT, main=paste("LRT"), xlab="LRT",ylab="Percentage of total",breaks=30) 
  histogram(geneData$LRT[which(data$LRT>quantile(data$LRT,0.95))], main=paste("(LRT)^1/4"), xlab="(LRT)^1/4", col=my.col[3],ylab="Percentage of total",breaks=30) 
}
histogram(data$Median.Expression, main=paste("Median Expression"), xlab="Median Expression",ylab="Percentage of total",breaks=30) 
histogram(data$Max.Expression, main=paste("Max Expression"), xlab="Max Expression",ylab="Percentage of total",breaks=30) 
histogram(data$X..GC.content, main=paste("%GC Content"), xlab="%GC Content",ylab="Percentage of total",breaks=30) 
histogram(geneData$X..GC.content, main=paste("log2(X..GC.content)"), xlab="log2(X..GC.content)", col=my.col[3],ylab="Percentage of total",breaks=50)
if(phylAge){  
  histogram(data$Phyletic.Age, main=paste("Phyletic Age"), xlab="Phyletic Age",ylab="Percentage of total",breaks=30)
}
#dev.copy2pdf(file=paste(folderAnalysis, organism,"Parameters", expDataSource,".pdf", sep=""),onefile=TRUE)#,paper="A4r"
dev.off()

### Save names of the used variables 
variableNames <- colnames(geneData)

### Calculating correlations
#*A simple way to compute the sample partial correlation for some data is 
#*to solve the two associated linear regression problems, get the residuals, 
#*and calculate the correlation between the residuals.

if (partial==TRUE){
  x <- data.frame(x1=NULL,x2=NULL,corValue=NULL,pValue=NULL,significant=NULL)
  for(j in variableNames){
    #j is the name of variable for which the correlation is calculated
    variablesToUse <- variableNames[variableNames != j] #all other variables
    t <- length(variableNames)-1
    
    for(n in c(1:t)){
      #*left part of the formula (~)
      j2 <- variablesToUse[n]
      variablesToUse2 <- variablesToUse[variablesToUse != j2]
      fmodel <-"geneData$"
      fmodel <- paste(fmodel,j,"~",sep="")
      fmodel2 <-"geneData$"
      fmodel2 <- paste(fmodel2,j2,"~",sep="")
      
      for(i in variablesToUse2){
        #*right part of the formula (~)
        fmodel <- paste(fmodel,"geneData$",i,"+",sep="")
        fmodel2 <- paste(fmodel2,"geneData$",i,"+",sep="")
      }
      
      fmodel <- substr(fmodel, 1, nchar(fmodel)-1)    # delete last character "+"
      fmodel2 <- substr(fmodel2, 1, nchar(fmodel2)-1) # delete last character "+"
      
      fmx <- glm(fmodel, na.action = na.exclude)
      fmy <- glm(fmodel2, na.action = na.exclude)
      
      xres <- resid(fmx)
      yres <- resid(fmy)
      
      ct <- cor.test(xres, yres, method=corMethod)          
      s <- ct$estimate
      coeff <- ct$p.value
      x <- rbind(x, data.frame(x1=j,x2=j2,corValue=s,pValue=coeff))
    }
  }
} else { #partial == FALSE
  x <- data.frame(x1=NULL,x2=NULL,corValue=NULL,pValue=NULL,significant=NULL)
  
  for(j in variableNames){
    #j is the name of variable for which the correlation is calculated
    variablesToUse <- variableNames[variableNames != j] #all other variables
    t <- length(variableNames)-1
    
    for(n in c(1:t)){
      j2 <- variablesToUse[n]
      variablesToUse2 <- variablesToUse[variablesToUse != j2]
      ct <- cor.test(geneData[,j], geneData[,j2], method=corMethod)
      s <- ct$estimate
      coeff <- ct$p.value
      x <- rbind(x, data.frame(x1=j,x2=j2,corValue=s,pValue=coeff,))
    }
  }
}

#For the test of significance -> adjust the p-value for multiple comparisons
# chosen method : Benjamini & Hochberg ; threshold : 0.05 for the adjusted p-val
x$pValAdjust <- p.adjust(x$pValue,method="BH")
x$significant <- as.integer(x$pValAdjust < 0.05)

row.names(x) <- c(1:nrow(x))
sameCorRowNumbers <- vector("numeric")

for(i in rownames(x)){
  sameCor <- x[with(x,x$x2 == x[i,]$x1 & x$x1 == x[i,]$x2),]
  
  if(as.integer(rownames(sameCor))>as.integer(i)){
    sameCorRowNumbers <- append(sameCorRowNumbers, as.integer(rownames(sameCor)))
  }
}

x <- x[-sameCorRowNumbers,]
saveX <- x

cat("\n Result of the correlation is saved in \"Cor_Original\".",sep="")
write.table(x,file=paste(folderAnalysis, organism, part, corMethod, "Cor", expDataSource, "Original.txt", sep=""),row.names = FALSE,quote = FALSE)

xp <- x[x$x1=="Omega",]
saveXp <- xp
xp <- xp[xp$significant >0,]
saveXp <- xp
xp$var <- xp$corValue*xp$corValue
v <- sum(xp$var)*100

print(paste("The variance of Omega is explained to ",v, "% through used parameters",sep=""))  

# ---------------------------------------------------------------------------------------
# PLOT THE CORRELATIONS: Cytoscape
# ---------------------------------------------------------------------------------------
### Making the file for Cytoscape
x$abs <- abs(x$corValue)*20
x$sign <- sign(x$corValue)
saveX2 <-x
x <- x[,c("x1","x2","significant","abs","sign")]
cat("\n Result of the correlation for Cytoscape representation is saved in \"Cor_List\".",sep="")
write.table(x, file=paste(folderAnalysis, organism, part, corMethod, "Cor", expDataSource,"List.txt",sep=""),row.names = FALSE,quote = FALSE)

### Draw graph in Cytoscape
x<-read.table(paste(folderAnalysis, organism, part, corMethod, "Cor", expDataSource,"List.txt",sep=""),header=T)
graphC <- x
cy <- CytoscapeConnection()
pluginVersion(cy)
#Initialize
g <- new ("graphNEL", edgemode = "undirected")
g <- initNodeAttribute (g, "nodeType", "char", "undefined")
g <- initNodeAttribute (g, "label", "char", "undefined")
g <- initEdgeAttribute (g, "edgeType", "char", "undefined")
g <- initEdgeAttribute (g, "significant", "char", "undefined")
g <- initEdgeAttribute (g, "sign", "char", "undefined")

#Add nodes and edges
#g <- addNode("title.node", g)
parameters <- unique(union(levels(graphC$x1),levels(graphC$x2)))
for (p in parameters){
  g <- addNode(p, g)
}
for (n in 1:length(rownames(graphC))){
  g <- addEdge(as.character(graphC[n,1]), as.character(graphC[n,2]), g)
}

#Add node and edge attribues
for (p in parameters){
  nodeData(g, p, "nodeType") <- p
  nodeData(g, p, "label") <- p
}
nodeData(g, parameters, "label") <- parametersNames
for (n in 1:length(rownames(graphC))){
  edgeData(g, as.character(graphC[n,1]), as.character(graphC[n,2]), "edgeType") = as.character(graphC[n,4])
  edgeData(g, as.character(graphC[n,1]), as.character(graphC[n,2]), "sign") = graphC[n,5]
  edgeData(g, as.character(graphC[n,1]), as.character(graphC[n,2]), "significant") = graphC[n,3]#edgeData(g, as.character(graphC[n,1]), as.character(graphC[n,2]), "label") = round(graphC[n,4]/20*graphC[n,5], digits=2)
}

#Create a CytoscapeWindow, after first making sure that no prior window of the same name
cy <- CytoscapeConnection()
setDefaultBackgroundColor(cy, my.col[1])
window.title = 'Correlation'
if (window.title %in% as.character (getWindowList(cy)))
  deleteWindow (cy, window.title)
cw <- new.CytoscapeWindow (window.title, g)

#Set window and network sizes
setWindowSize (cw, 1200, 1200)
fitContent (cw)
setZoom (cw, 0.9 * getZoom (cw))

#Send graph to Cytoscape
displayGraph (cw)

#Set default settings for the graph
setDefaultEdgeColor(cw, my.col[2])
lockNodeDimensions(cw, FALSE)
setNodeShapeDirect(cw, parameters, "ellipse")
setNodeFontSizeDirect(cw, parameters, 10)
setNodeColorDirect(cw, parameters, my.col[7])
setNodeWidthDirect(cw, parameters, 85)
setNodeHeightDirect(cw, parameters, 40)

#Ask Cytoscape to layout the graph
layoutNetwork (cw, 'attribute-circle')

#Instruct Cytoscape to use each node's 'label' attribute as the value for the visible label it draws on the node
setNodeLabelRule (cw, 'label')
setEdgeLineWidthRule(cw, "edgeType", as.character(graphC$abs), as.numeric(graphC$abs))
setEdgeColorRule(cw, "sign", c("-1", "1"), c("blue", "red"), mode="lookup")
setEdgeOpacityRule(cw, "significant", c("1", "0"), c("175", "0"), mode="lookup")
#setEdgeLabelRule(cw, "label")

#Ask Cytoscape to redraw the graph using these rules
redraw (cw)
#saveLayout(cw,'CorrelationLayout13')     #manually change the order of parameters
#restoreLayout(cw, 'CorrelationLayout13')
  #if you change manually the plot in the cytoscape window
  # => then "saveLyout(cw, 'name_you_want')
  #the next time you run, if you want your layout previously defined then
  #comment the saveLayout line and uncomment the restoreLayout
  #-> it will give the same layout as you chose manually

#Display the current graph.
fitContent (cw)

#Save the graph
filename=paste(folderAnalysis,organism, part, corMethod, "Cor", expDataSource, "Cyt.png", sep="")
####  folder ~Bureau/stage15... doesn't work, use /home/user/Bureau/stage15...
filename <- sub("~","/home/user",filename)   #=> maybe to change according to your path
saveImage (cw, filename, 'png')

# ---------------------------------------------------------------------------------------
# COMPUTE CORRELATIONS : Circos plot (circlize package)
# ---------------------------------------------------------------------------------------
#####*****Organism partial correlation with GLM model for Circlize*****#####
cat("\nCalculating partial correlation with glm model for each tissue separately.",sep="")

##### Load the data
dataOrg <- read.table(paste(folderAnalysis, organism, "TableStages", expDataSource,".txt",sep=""), header=TRUE)
cat("\nSummary of the data (", nrow(dataOrg), " genes) in the first step : ", sep="")
summary(dataOrg)
dataOrg <- dataOrg[dataOrg$Max.Expression>0.00015,]

cat("\nTo calculate correlation ", corMethod, " correlation was used.",sep="")
if (partial){
  cat("Partial correlation was performed.",sep="")
  part <- "Partial"
} else {
  cat("Normal correlation was performed.",sep="")
  part <- "Normal"
}

cat("\nThe analysis is done for ", length(devStageNames), " developmental stages.", sep="")

##### Select needed columns 
if(organism == "Caeno"){ #no Selectome data for Caeno (LRT and P.1 not available)
  if(phylAge){
    dataOrg <- dataOrg[,c("Omega.0", "CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", "Phyletic.Age", devStageOutNames)]
    parameterNames <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number", "Phyletic.Age")
  } else if (!phylAge){
    dataOrg <- dataOrg[,c("Omega.0","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", devStageOutNames)]
    parameterNames <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number")  
  }  
} else {
  if(phylAge){
    dataOrg <- dataOrg[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", "Phyletic.Age", devStageOutNames)]
    parameterNames <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number", "Phyletic.Age")
  } else if (!phylAge){
    ### we don't have phyletic age for the moment !!! and no Omega
    dataOrg <- dataOrg[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", devStageOutNames)]
    parameterNames <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number")  
  }
}

colnames(dataOrg) <- c(parameterNames, devStageOutNames)
geneDataOrg <- na.omit(dataOrg)
cat("\nAll the genes with unknown parameters are removed from the analysis. ", nrow(geneDataOrg)," are left for the analysis.", " Summary: ", sep="")
summary(geneDataOrg)

##### Normalization of the data (function defined in "scriptFunctions.R")
geneDataOrg <- normalizeData2(geneDataOrg)

##### Draw representation of the expression data
cat("\nGraphical representation of the expression data is saved in the file \"DevStagesExpression\".", sep="")
x11()
dev.new(height=9, width=12)
par(cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
palette(rev(rich.colors(length(devStageNames)+2)))

plot(density(geneDataOrg[,devStageOutNames[1]],n=1000), main = "Expression values among different stages",xlab="Log of normalized signal intensity",col=(1), lwd=3)

for(i in c(2:length(devStageNames))){
  lines(density(geneDataOrg[,devStageOutNames[i]],n = 1000), col=(i), lwd=3)
}
legend("topright",devStagePrintNames,col=(1:length(devStageOutNames)),lty="solid", lwd=3)

dev.copy2pdf(device=quartz, file=paste(folderAnalysis, organism,"StageExpression", expDataSource,".pdf", sep=""),onefile=TRUE)#,paper="A4r"
#dev.off()

cat("\nOverall ",nrow(geneDataOrg)," genes were used for analysis.","\nSummary:",sep="")
summary(geneDataOrg)

##### Calculation correlation
variableNames <- devStageOutNames
x <- data.frame(x1=NULL,x2=NULL,corValue=NULL,pValue=NULL,significant=NULL)
for(j in variableNames){
  #j is the name of variable for which the correlation is calculated
  variablesToUse <- parameterNames #Names of other variables
  t = length(variablesToUse)
  for(n in c(1:t)){
    j2 <- variablesToUse[n]
    variablesToUse2 <- variablesToUse[variablesToUse != j2]
    
    if(partial==TRUE){
      fmodel <-"geneDataOrg$"
      fmodel <- paste(fmodel,j,"~",sep="")
      fmodel2 <-"geneDataOrg$"
      fmodel2 <- paste(fmodel2,j2,"~",sep="")
      for(i in variablesToUse2){
        fmodel <- paste(fmodel,"geneDataOrg$",i,"+",sep="")
        fmodel2 <- paste(fmodel2,"geneDataOrg$",i,"+",sep="")
      }
      fmodel <- substr(fmodel, 1, nchar(fmodel)-1) #delet last character "+"
      fmodel2 <- substr(fmodel2, 1, nchar(fmodel2)-1) #delet last character "+"
      fmx <- glm(fmodel, na.action = na.exclude)
      fmy <- glm(fmodel2, na.action = na.exclude)
      xres <- resid(fmx)
      yres <- resid(fmy)
      ct <- cor.test(xres, yres, method=corMethod)      
    } else {
      ct <- cor.test(geneDataOrg[,j], geneDataOrg[,j2], method=corMethod)
    }
    s <- ct$estimate
    coeff <- ct$p.value
    x <- rbind(x, data.frame(x1=j,x2=j2,corValue=s,pValue=coeff))
  }
}

saveX3a <- x
variableNames <- parameterNames
for(j in variableNames) {
  #j is the name of variable for which the correlation is calculated
  variablesToUse <- variableNames[variableNames != j] #Names of other variables
  t = length(variablesToUse)
  
  for(n in c(1:t)){
    j2 <- variablesToUse[n]
    variablesToUse2 <- variablesToUse[variablesToUse != j2]
    
    if(partial==TRUE){
      fmodel <-"geneDataOrg$"
      fmodel <- paste(fmodel,j,"~",sep="")
      fmodel2 <-"geneDataOrg$"
      fmodel2 <- paste(fmodel2,j2,"~",sep="")
      
      for(i in variablesToUse2){
        fmodel <- paste(fmodel,"geneDataOrg$",i,"+",sep="")
        fmodel2 <- paste(fmodel2,"geneDataOrg$",i,"+",sep="")
      }
      
      fmodel <- substr(fmodel, 1, nchar(fmodel)-1) #delet last character "+"
      fmodel2 <- substr(fmodel2, 1, nchar(fmodel2)-1) #delet last character "+"
      fmx <- glm(fmodel, na.action = na.exclude)
      fmy <- glm(fmodel2, na.action = na.exclude)
      xres <- resid(fmx)
      yres <- resid(fmy)
      ct <- cor.test(xres, yres, method=corMethod)
    } else {
      ct <- cor.test(geneDataOrg[,j], geneDataOrg[,j2], method=corMethod)
    }
    s <- ct$estimate
    coeff <- ct$p.value
    x <- rbind(x, data.frame(x1=j,x2=j2,corValue=s,pValue=coeff))
  }
}

#For the test of significance -> adjust the p-value for multiple comparisons
#Chosen method : Benjamini & Hochberg ; threshold : 0.05 for the adjusted p-val
x$pValAdjust <- p.adjust(x$pValue,method="BH")
x$significant <- as.integer(x$pValAdjust < 0.05)

saveX3b <- x

row.names(x) <- c(1:nrow(x))
x$x2 <- factor(x$x2,levels=levels(x$x1))
sameCorRowNumbers <- vector("numeric")

for(i in rownames(x)){
  sameCor <- x[with(x,x$x2 == x[i,]$x1 & x$x1 == x[i,]$x2),]   #remove correlation with itself
  if(nrow(sameCor)>0){
    if(as.integer(rownames(sameCor))>as.integer(i)){
      sameCorRowNumbers <- append(sameCorRowNumbers, as.integer(rownames(sameCor)))  #remove the 2n of a1xa2 [...] a2xa1
    }
  }
}

x <- x[-sameCorRowNumbers,]
saveX4 <- x    # 99rows

#Save the file
cat("\n Correlation table",sep="")
cat("\n Result of the correlation is saved in \"CorStages_Original\".",sep="")
write.table(x,file=paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""),row.names = FALSE,quote = FALSE)

# ---------------------------------------------------------------------------------------
# PLOT THE CORRELATIONS: circlize 
# ---------------------------------------------------------------------------------------
### Reshape data: from dataframe to matrix
x <- read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""), header=T)
corrTab <- x
corrTab <- corrTab[,c("x1","x2","corValue")]
corrTab$abs <- abs(corrTab$corValue)
corrTab <- corrTab[,c("x1","x2","abs")]

#Function defined in "scriptFunctions.R")
#-> takes the data frame in argument and returns a symmetric matrix
mat <- dfToMatrix(corrTab)  
diag(mat) <- 0      # set the diagonal values to 0 because we don't want to draw it

#Transparency : 0.8 transparency if not significant, else 0.3 transparency
xb <- read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""), header=T)
corrTabTransp <- xb
corrTabTransp <- corrTabTransp[,c("x1","x2","significant")]
matTransp <- dfToMatrix(corrTabTransp)
#0 not significant, 1 significant : if x==0 -> 1/2 transparent, if x==1 => 0 transparency
matTransp <- apply(matTransp,c(1,2), function(x) {ifelse(x==0, 0.8,0.3)})  
diag(matTransp) <- 1  #transparent on the diagonal  (correlation with itself -> transparent)

#Color: blue if negative correlation, red if positive correlation
xc <- read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""), header=T)
corrTabSign <- xc
corrTabSign$sign <- sign(corrTabSign$corValue)  #-1 if negative, 0 if zero, +1 if positive
corrTabSign <- corrTabSign[,c("x1","x2","sign")]
matSign <- dfToMatrix(corrTabSign)
matSign <- apply(matSign,c(1,2), function(x) { ifelse(x<0,"#1874CD","#FF3030")})
matSign <- apply(matSign,c(1,2), function(x) {toString(x)})

### Circlize the matrix and draw it
# matSave <- mat

#Format the names for the graph
if(organism=="Mus"){
  paramNames <- colnames(mat)[regexpr("stage",colnames(mat))<0]
  paramNames <- gsub("\\.", " ",paramNames)
  paramNames <- gsub(" {2}", "%", paramNames) # the only 2 blank spaces -> % GC content
  colnames(mat)[regexpr("stage",colnames(mat))<0] <- paramNames
  rownames(mat)[regexpr("stage",colnames(mat))<0] <- paramNames
  
  stageNames <- colnames(mat)[regexpr("stage",colnames(mat))>0]
  stageNames <- gsub("\\.adult\\.stage", "\nadult stage", stageNames)
  stageNames <- sub("\\.stage","\nstage",stageNames,perl=FALSE)
  stageNames <- sub("stage\\.","stage ",stageNames,perl=FALSE)
  stageNames <- gsub("\\.", "-",stageNames)
  colnames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames
  rownames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames

} else if (organism=="Danio"){
  paramNames <- colnames(mat)[regexpr("stage",colnames(mat))<0]
  paramNames <- gsub("\\.{2}", "%", paramNames)
  paramNames <- gsub("\\.", " ",paramNames)
  colnames(mat)[regexpr("stage",colnames(mat))<0] <- paramNames
  rownames(mat)[regexpr("stage",colnames(mat))<0] <- paramNames
  
  stageNames <- colnames(mat)[regexpr("stage",colnames(mat))>0]
  stageNames <- sub("\\.", "\n",stageNames)
  stageNames <- gsub("\\.{2}", "% ", stageNames)
  stageNames <- gsub("([0-9]).([0-9])","\\1-\\2",stageNames)
  stageNames <- gsub("\\.", " ",stageNames)
  stageNames <- gsub("Segmentation ", "Segmentation\n", stageNames)
  stageNames <- gsub("Gastrula ", "Gastrula\n", stageNames)
  stageNames <- gsub("Hatching ", "Hatching\n", stageNames)
  stageNames <- gsub("Pharyngula ", "Pharyngula\n", stageNames)
  colnames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames
  rownames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames
} else if (organism=="Drosophila" & !FPKM){
  paramNames <- colnames(mat)[regexpr("stage",colnames(mat))<0]
  paramNames <- gsub("\\.{2}", "%", paramNames)
  paramNames <- gsub("\\.", " ",paramNames)
  colnames(mat)[regexpr("stage",colnames(mat))<0] <- paramNames
  rownames(mat)[regexpr("stage",colnames(mat))<0] <- paramNames
  
  stageNames <- colnames(mat)[regexpr("stage",colnames(mat))>0]
  stageNames <- gsub("\\.", " ",stageNames)
  stageNames <- sub(" germ band","\ngerm band", stageNames)
  stageNames <- sub("embryonic ","embryonic\n", stageNames)
  stageNames <- sub("instar larval","instar\nlarval", stageNames)
  colnames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames
  rownames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames
} else if(organism=="Xenopus"){
  paramNames <- colnames(mat)[regexpr("stage",colnames(mat))<0]
  paramNames <- gsub("\\.", " ",paramNames)
  paramNames <- gsub(" {2}", "%", paramNames) # the only 2 blank spaces -> % GC content
  colnames(mat)[regexpr("stage",colnames(mat))<0] <- paramNames
  rownames(mat)[regexpr("stage",colnames(mat))<0] <- paramNames
  
  stageNames <- colnames(mat)[regexpr("stage",colnames(mat))>0]
  stageNames <- gsub("\\.{2}", "-",stageNames)  
  stageNames <- gsub("\\.", " ",stageNames)  
  stageNames <- sub("(\\d)-(\\d{2}) cell ","\n\\1-\\2 cell",stageNames)
  stageNames <- sub("stage cell (\\d)","\\1 cell stage",stageNames)
  colnames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames
  rownames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames
} else if(organism=="Caeno"  & !FPKM){
  paramNames <- colnames(mat)[regexpr("stage",colnames(mat))<0]
  paramNames <- gsub("\\.", " ",paramNames)
  paramNames <- gsub(" {2}", "%", paramNames) # the only 2 blank spaces -> % GC content
  colnames(mat)[regexpr("stage",colnames(mat))<0] <- paramNames
  rownames(mat)[regexpr("stage",colnames(mat))<0] <- paramNames
  
  stageNames <- colnames(mat)[regexpr("stage",colnames(mat))>0]
  stageNames <- gsub("\\.", " ",stageNames)  
  stageNames <- sub(" stage", "\nstage",stageNames)
  colnames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames
  rownames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames  

} else if(organism=="Drosophila" && FPKM){
  stageNames <- colnames(mat)[regexpr("stage",colnames(mat))>0]
  stageNames <- gsub(".stage","", stageNames)
  stageNames <- gsub("([0-9])\\.([0-9])", "\\1\\-\\2", stageNames)
  stageNames <- gsub("\\."," ", stageNames)
  stageNames <- gsub("Embryo", "Embryo ", stageNames)
  stageNames <- gsub("L3PS", "L3PS ", stageNames)
  colnames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames
  rownames(mat)[regexpr("stage",rownames(mat))>0] <- stageNames  
  paramNames <- colnames(mat)[!(colnames(mat)%in% stageNames)]  

} else if(organism=="Caeno" && FPKM){
  stageNames <- colnames(mat)[regexpr("stage",colnames(mat))>0]
  stageNames <- gsub(".stage","", stageNames)
  stageNames <- gsub("([0-9])\\.([0-9])", "\\1\\-\\2", stageNames)
  stageNames <- gsub("\\_"," ", stageNames)
  
  colnames(mat)[regexpr("stage",colnames(mat))>0] <- stageNames
  rownames(mat)[regexpr("stage",rownames(mat))>0] <- stageNames
  paramNames <- colnames(mat)[!(colnames(mat)%in% stageNames)]  
}

#Set color
grid.col=NULL
grid.col[paramNames] <- colorRampPalette(c("darkolivegreen","darkolivegreen4"), alpha = TRUE)(length(paramNames))
grid.col[stageNames] <- colorRampPalette(c("yellow","yellow3"), alpha = TRUE)(length(stageNames))

#Graphical settings for circlize plot
x11()
png(filename=paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Circlize.png",sep=""),
    width=1000, height=1000)
par(bg = "black")
# circos.par(start.degree = 90, gap.degree=5)

circos.par(start.degree = 90, gap.degree=5, "canvas.xlim" = c(-1.2,1.2), "canvas.ylim" = c(-1.2,1.2))
  #for some organisms with longer stage name, need to adjust canvas.xlim and canvas.ylim

#Draw "chord diagram"
chordDiagram(mat, symmetric=TRUE, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1),
             order=c(stageNames,paramNames), row.col=matSign,
             grid.col=grid.col, transparency=matTransp)

first <- paramNames[1]
first <- append(first, stageNames[1])
last <- paramNames[length(paramNames)]
last <- append(last, stageNames[length(stageNames)])
for(i in seq_along(first)) {
  start.degree = get.cell.meta.data("cell.start.degree", sector.index = first[i], track.index = 1)
  end.degree = get.cell.meta.data("cell.end.degree", sector.index = last[i], track.index = 1)
  rou1 = get.cell.meta.data("cell.bottom.radius", sector.index = first[i], track.index = 1)
  rou2 = get.cell.meta.data("cell.top.radius", sector.index = last[i], track.index = 1)
  draw.sector(start.degree, end.degree, rou1, rou2, border = NA, col = "#444444")
}
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% stageNames) {          #or (clearer) : if(regexpr("stage",sector.name)>0 not ok for FPKM
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5),col="white", cex=1.2)
  } else {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",
                niceFacing = TRUE, adj = c(0.5, 0),col="white", cex=1.2)
  }
}, bg.border = NA)
dev.off()
circos.clear()

# ---------------------------------------------------------------------------------------
# Other graphs and tests
# ---------------------------------------------------------------------------------------
#### Evolution of 1 parameter through development 
#Parameters to set if this part of the code is run separately
organism <- organism
folder <- folder
folderAnalysis <- folderAnalysis
part <- part
corMethod <- corMethod
expDataSource <- expDataSource
FPKM <- FPKM
x <- "Omega"  #Parameter to plot
data<-read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""),header=T)
data <- data[which(regexpr("stage", data$x1)>0 & regexpr(x, data$x2)>0),] 
data$significant <- as.factor(data$significant) #need as factor for the colour below

#Color as follows :
# green -> zygote, cleavage
# orange -> blastula, gastrula
# blue -> neurula
# red -> organogenesis, segmentation
# black -> post-embryonic


colLabel <- graphParamConti(organism)$colLabel
speciesName <- graphParamConti(organism)$speciesName

        # if(organism=="Mus"){
        #   colLabel <- c(rep("green",3), rep("orange",3), rep("blue", 1), rep("red",5), rep("black", 3))
        # }else if(organism=="Danio"){
        #   colLabel <- c(rep("green",1), rep("orange", 2), rep("blue", 2), rep("red",4), rep("black", 5))
        # }else if(organism=="Drosophila" & !FPKM){
        #   colLabel <- c(rep("green",1), rep("orange", 2), rep("blue",2), rep("red", 2),rep("black", 2))
        # }else if(organism=="Caeno" & !FPKM){
        #   colLabel <- c(rep("orange", 1), rep("blue",1), rep("red", 1),rep("black", 2))
        # }else if(organism=="Xenopus"){
        #   colLabel <- c(rep("green",5), rep("orange", 4), rep("blue",2), rep("red", 3))
        # }else if(organism=="Drosophila" & FPKM){
        #   colLabel <- c(rep("green",1), rep("orange",1), rep("blue", 1), rep("red",9), rep("black", 15))
        # } else if(organism=="Caeno" & FPKM){
        #   colLabel <- c(rep("green",2), rep("orange",8), rep("blue", 4), rep("red",10), rep("black", 11))
        # }
data$x1 <- factor(as.character(data$x1),levels=data$x1)
tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA))
x11()
ggplot(data, aes(x=x1, y=corValue,group=1,colour=significant)) +
  #   geom_point(aes(colour=significant),size=4) + geom_line(colour="black")+
  geom_point(size=4) + geom_line(colour="black")+
  #   scale_colour_gradient(low="white", high="chartreuse4",guide="none")+
  #   scale_colour_manual(breaks=c("0","1"), values=c("chocolate","green"),guide="none")+
  scale_colour_manual(values=c("0"="chocolate","1"="green"),guide="none")+
  scale_y_continuous("Spearmans' (partial) correlation coefficient")+
#   scale_x_discrete(paste0(x," - ",organism))+
  scale_x_discrete(x)+
  ggtitle(tit)+
  theme(title = element_text(lineheight=0.8, face="bold"),
        axis.title.y = element_text(face="bold", colour="#990000", size=10),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
  theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())

ggsave(filename=paste0(folderAnalysis,x,"Plot",organism, ".png"))



#### Plot 1 parameter with continuous x-axis  --------------------------------------------
plotOneParaConti <- function(x){
  data<-read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""),header=T)
  data <- data[which(regexpr("stage", data$x1)>0 & regexpr(x, data$x2)>0),] 
  data$significant <- as.factor(data$significant) #need as factor for the colour below
  
  #class defined elsewhere
  colLabel <- graphParamConti(organism)$colLabel
  speciesName <- graphParamConti(organism)$speciesName
  stageNum <- graphParamConti(organism)$stageNum
  tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA))
  data <- data[c(1:length(colLabel)),]
  stagesName <- data$x1[1:length(colLabel)]
  data$x1 <- as.numeric(as.character(stageNum))
  
  x11()
  ggplot(data, aes(x=x1, y=corValue,group=1,colour=significant)) +
    #   geom_point(aes(colour=significant),size=4) + geom_line(colour="black")+
    geom_point(size=4) + geom_line(colour="black")+
    #   scale_colour_gradient(low="white", high="chartreuse4",guide="none")+
    #   scale_colour_manual(breaks=c("0","1"), values=c("chocolate","green"),guide="none")+
    scale_colour_manual(values=c("0"="chocolate","1"="green"),guide="none")+
    scale_y_continuous("Spearmans' (partial) correlation coefficient")+
    #   scale_x_discrete(paste0(x," - ",organism))+
    scale_x_continuous(x, breaks=data$x1, labels=stagesName)+
    ggtitle(tit)+
    theme(title = element_text(lineheight=0.8, face="bold"),
          axis.title.y = element_text(face="bold", colour="#990000", size=10),
          axis.text.y = element_text(colour="black"),
          axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
    theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())
  
  ggsave(filename=paste0(folderAnalysis,x,"ContinuousPlot",organism, ".png"))
}
# x <- "Omega"  #Parameter to plot
plotOneParaConti("Omega")
plotOneParaConti("Phyletic.Age")
plotOneParaConti("Paralogs.Number")
plotOneParaConti("CDS.Length")
if(organism!="Caeno"){
  plotOneParaConti("P.1")
  plotOneParaConti("LRT")
}

##### Essentiality and Omega in Violin Plot & T-test --------------------------------------
if(essentiality){
    data <- read.table(paste(folderAnalysis,organism, "Table", expDataSource, ".txt",sep=""), header=TRUE)
    data <- data[data$Max.Expression>0.00015,]
    if(organism=="Caeno"){ #we don't have the Selectome values
      if(phylAge){
        data <- data[,c("Omega.0", "CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                        "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number", "Phyletic.Age", "Essentiality")]
        
        colnames(data) <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                            "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number","Phyletic.Age", "Essentiality")
        
      } else if(!phylAge){
        data <- data[,c("Omega.0","CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                        "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number", "Essentiality")]
        
        colnames(data) <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                            "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number", "Essentiality")  
      }
    } else {
      if(phylAge){
        data <- data[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                        "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number", "Phyletic.Age", "Essentiality")]
        
        colnames(data) <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                            "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number","Phyletic.Age", "Essentiality")
      
      } else if(!phylAge){
        data <- data[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                        "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number", "Essentiality")]
        
        colnames(data) <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", "Median.Expression", 
                            "Max.Expression", "Tau", "X..GC.content", "Paralogs.Number", "Essentiality")  
      }
    }
    
    data <- na.omit(data)
    geneData<-data
    
    ##### Normalization of the data                           #=> cf. function ?
    
    geneData <- normalizeData(geneData)

    summary(geneData)
    
    ##### Names of the variables used
    variableNames <- colnames(geneData)
    variablesToUse <- variableNames[variableNames != c("Omega") & variableNames != c("Essentiality")]
    #following does not work :
    #variablesToUse <- variableNames[variableNames != c("Omega", "Essentiality")] #all other variables 
    
    fmodel <-"geneData$"
    fmodel <- paste(fmodel,"Omega","~",sep="")
    
    for(i in variablesToUse){
    fmodel <- paste(fmodel,"geneData$",i,"+",sep="")
    }
    
    fmodel <- substr(fmodel, 1, nchar(fmodel)-1) #delet last character "+"
    fmx <- glm(fmodel, na.action = na.exclude)
    xres <- resid(fmx)
    
    dataPlot <- data.frame(Omega=xres, Essentiality=geneData$Essentiality)
    dataP <- as.matrix(dataPlot[,c("Omega", "Essentiality")])
    
    dev.off() #-> added : palette() did not work, dev.off()  and after it runs
    dev.new(height=8, width=8) 
    X11()
    png(paste0(folderAnalysis, organism, "OmegaEssentiality.png"))
    palette(rainbow(9))
    trellis.par.set(list(background=list(col=my.col[1]), add.text=list(col=my.col[2], cex=1.5),axis.line=list(col=my.col[2]), axis.text=list(col=my.col[2], cex=1.5),
                       par.main.text=list(col=my.col[2], cex=1.3), par.xlab.text=list(col=my.col[2], cex=1.5), par.ylab.text=list(col=my.col[2], cex=1.7), plot.line=list(col=my.col[2]), dot.line=list(lwd=1,
                                                                                                                                                                                                        lty=2, col="#4B4B4B"))) #trellis.par.get()
    
    bwplot(as.numeric(dataP[,1])~dataP[,2], xlab="", ylab="residuals of log2(Omega)", main=paste("Distribution of Omega according to essentiality",sep=""), horizontal=FALSE, col = c("#00BFFF"), fill=c("blue"),panel = function(x,y,..., box.ratio, col, pch){
    panel.violin(x=x, y=y,..., cut = 0, varwidth = TRUE, box.ratio = 4*box.ratio, col=col)
    panel.bwplot(x=x, y=y, ..., varwidth = TRUE ,box.ratio = .5, pch='|', notch=TRUE)
    },
    par.settings = list(box.rectangle=list(col=my.col[2], lwd=2), plot.symbol = list(pch='.', cex = 0.1, col=my.col[2]), box.umbrella=list(col=my.col[2])), scales=list(x=list(rot=10,
                                                                                                                                                                             labels=c("Not Essential", "Essential"))))
    dev.off()
    dev.copy2pdf(device=quartz, file=paste(folderAnalysis, organism,"OmegaEssentialityNewAres.pdf", sep=""),onefile=TRUE)#,paper="A4r"
    
    
    dataPF <- data.frame(dataPlot)
    dataPF$Omega <- as.numeric(dataPF$Omega)
    
    t.test(Omega~Essentiality, data=dataPF)
    t.test(Omega~Essentiality, data=dataPF, alternative = "greater")
}
# ################################################
# ################################################
# # boxplot for the experiment for mus & drosophila (merge 2 datasets)

# if(organism=="Mus"){
# orgExpression1 <- read.delim("musEMEXP51.tsv") #  (Affymetrix data1 from Bgee)
# orgExpression2 <- read.delim("musEMTAB368.tsv") #  (Affymetrix data2 from Bgee)
# 
# orgExpression1 <- orgExpression1[,c("Gene.ID","Stage.ID","Stage.name","Log.of.normalized.signal.intensity")] 
# colnames(orgExpression1) <- c("Ensembl.Gene.ID","Stage.ID","Stage.name","Log.of.normalized.signal.intensity")
# 
# orgExpression2 <- orgExpression2[,c("Gene.ID","Stage.ID","Stage.name","Log.of.normalized.signal.intensity")] 
# colnames(orgExpression2) <- c("Ensembl.Gene.ID","Stage.ID","Stage.name","Log.of.normalized.signal.intensity")
# 
# dat2 <- rbind(orgExpression1,orgExpression2)
# 
# orgExpression <- dat2
# orgExpression <- dcast(orgExpression, Ensembl.Gene.ID~Stage.name, value.var="Log.of.normalized.signal.intensity")
# orgExpression<-orgExpression[,c("Ensembl.Gene.ID","zygote stage", "Theiler stage 02 (mouse)", "Theiler stage 03 (mouse)",
#                                 "Theiler stage 04 (mouse)", "Theiler stage 05 (mouse)", "Theiler stage 06 (mouse)",
#                                 "neurula stage", "Theiler stage 13 (mouse)", "Theiler stage 15 (mouse)", 
#                                 "Theiler stage 17 (mouse)", "Theiler stage 20 (mouse)", "Theiler stage 22 (mouse)",
#                                 "Theiler stage 24 (mouse)","Theiler stage 26 (mouse)","post-juvenile adult stage")]
# stageNames <-colnames(orgExpression[-1])
# 
# x11()
# png(paste0("/home/user/Bureau/stage15/mouse/",organism,"BOXbeforeQNbwplot.png"))
# bwplot(Log.of.normalized.signal.intensity ~ Stage.name, data=dat2, 
#        scales = list(x = list(rot=90)))
# dev.off()
# # 
# # boxplot(dat[-c(1,16:20)])
# 
# orgRPKM <- na.omit(orgExpression)
# tempRPKM <- orgRPKM[,-1]
# tempRPKM_m <- as.matrix(tempRPKM)
# tempRPKM <- normalize.quantiles(tempRPKM_m)
# orgRPKM[,-1] <- data.frame(tempRPKM)
# x11()
# png(paste0("/home/user/Bureau/stage15/mouse/",organism,"BOXafterQNplot.png"))
# boxplot(orgRPKM[,-1], las=2, cex.lab=0.5)
# dev.off()
# }else if (organism=="Drosophila"){
#   orgExpression1 <- read.delim("drosoEMTAB379.tsv") #  (Affymetrix data from Bgee, exp1)
# orgExpression2 <- read.delim("drosoGSE3955.tsv") # (Affymetrix data from Bgee, exp2)
#   
#   orgExpression1 <- orgExpression1[,c("Gene.ID","Stage.ID","Stage.name","Log.of.normalized.signal.intensity")] 
#   colnames(orgExpression1) <- c("Ensembl.Gene.ID","Stage.ID","Stage.name","Log.of.normalized.signal.intensity")
#   
#   orgExpression2 <- orgExpression2[,c("Gene.ID","Stage.ID","Stage.name","Log.of.normalized.signal.intensity")] 
#   colnames(orgExpression2) <- c("Ensembl.Gene.ID","Stage.ID","Stage.name","Log.of.normalized.signal.intensity")
#   
#   dat2 <- rbind(orgExpression1,orgExpression2)
#   
#   orgExpression <- dat2
#   orgExpression <- dcast(orgExpression, Ensembl.Gene.ID~Stage.name, value.var="Log.of.normalized.signal.intensity")
#   
#   orgExpression <- orgExpression[,c("Ensembl.Gene.ID","embryonic stage 2 (Drosophila)", "embryonic stage 4 (Drosophila)",
#                                     "embryonic stage 5 (Drosophila)", "early extended germ band stage (Drosophila)",
#                                     "late extended germ band stage (Drosophila)", "embryonic stage 10 (Drosophila)", 
#                                     "embryonic stage 17 (Drosophila)", "first instar larval stage (Drosophila)", 
#                                     "second instar larval stage (Drosophila)")]   
#   stageNames <-colnames(orgExpression[-1])
#   
#   x11()
#   png(paste0("/home/user/Bureau/stage15/drosophila/BOXBeforeQNbwplot.png"))
#   bwplot(Log.of.normalized.signal.intensity ~ Stage.name, data=dat2, 
#          scales = list(x = list(rot=90)))
#   dev.off()
#   
#   orgRPKM <- na.omit(orgExpression)
#   tempRPKM <- orgRPKM[,-1]
#   tempRPKM_m <- as.matrix(tempRPKM)
#   tempRPKM <- normalize.quantiles(tempRPKM_m)
#   orgRPKM[,-1] <- data.frame(tempRPKM)
# x11()  
#   png(paste0("/home/user/Bureau/stage15/drosophila/BOXafterQNplot.png"))
# boxplot(orgRPKM[,-1], las=2, cex.lab=0.5)
#   dev.off()
#  
#}

# ---------------------------------------------------------------------------------------
# Which regression curve has the best fit ? 
# ---------------------------------------------------------------------------------------
#First- or second-degree polynomial has the best fit ?
#Take only stages until birth (e.g. larval, post-juvenil stages are discarded)
#Criterion best fit : ANOVA test

#Parameters to set if this part of the code is run separately
folderAnalysis <- folderAnalysis
organism <- organism
part <- part
corMethod <- corMethod
expDataSource <- expDataSource
symbRNA <- symbRNA
FPKM <- FPKM

dataCorr <- read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""), header=T)

#Build model for the following parameters
param1 <- "Omega"
param2 <- "Paralogs.Number"
param3 <- "Phyletic.Age"

dataCorr <- dataCorr[,c("x1", "x2", "corValue", "significant")]
dataCorr1 <- dataCorr[regexpr("stage", dataCorr$x1)>0 & regexpr(param1, dataCorr$x2)>0,]
dataCorr2 <- dataCorr[regexpr("stage", dataCorr$x1)>0 & regexpr(param2, dataCorr$x2)>0,]
dataCorr3 <- dataCorr[regexpr("stage", dataCorr$x1)>0 & regexpr(param3, dataCorr$x2)>0,]

stageName <- as.vector(dataCorr1$x1)
stageName <- as.vector(dataCorr2$x1)
stageName <- as.vector(dataCorr3$x1)

#Set parameters using a class defined in "scriptFunctions.R"
stageNum <- graphParamModel(organism)$stageNum
sc <- graphParamModel(organism)$sc
speciesName <- graphParamModel(organism)$speciesName

                  #If not use the class:
                  # if(organism == "Mus"){  #scale : dpc
                  #   #stages : 
                  #   stageNum <- c(0, 1, 2, 3, 4, 4.5, 7.5, 8.5, 9.5, 10.5, 12, 14, 16, 18)
                  #   sc <- "Days post conception"
                  #   speciesName <- "M. musculus"
                  # } else if(organism == "Danio"){  #scale : hours of dvpmt at 28.5°C
                  # #     stageNum <- c(0, 6, 8, 9, 10, 11.66, 16, 24, 30, 48, 96, 120, 336, 720)
                  #     stageNum <- c(0, 6, 8, 9, 10, 11.66, 16, 24, 30, 48)  
                  #   sc <- "Hours of development"
                  #   speciesName <- "D. rerio"
                  # } else if(organism == "Drosophila" & !FPKM){ #scale : minutes after fertilization
                  #     stageNum <- c(45, 105, 150, 270, 290, 440, 1200)
                  # #     stageNum <- (stageNum%/%60) + ((stageNum%%60)/60)
                  #   sc <- "Minutes after fertilization"
                  #   speciesName <- "D. melanogaster"
                  # } else if(organism == "Xenopus"){  #hours post-fertilization
                  #   stageNum <- c(1.5, 2, 2.25, 2.75, 3, 5, 7, 8, 9, 17.5, 20.75, 32.5, 44.5, 66)
                  #   sc <- "Hours post-fertilization"
                  #   speciesName <- "X. tropicalis"
                  # } else if(organism == "Drosophila" & FPKM){  #hours post-fertilization
                  # #     stageNum <- c(seq(1, 23,2), 36.5, 60.5, 96)
                  # #     stageNum <- c(seq(1, 23,2), 36.5, 60.5) #try without last stage
                  #   stageNum <- c(seq(1, 23,2)) #with only 1 post-embryonic stage
                  #   stageNum <- stageNum*60  #have it in minutes like other drosophila data
                  #   sc <- "Minutes after fertilization"
                  #   speciesName <- "D. melanogaster"
                  # } else if(organism == "Caeno" & FPKM){  #hours post-fertilization
                  #   stageNum <- c(seq(0,240,30),seq(300, 720,30)) 
                  #   sc <- "Minutes after fertilization"
                  #   speciesName <- "C. elegans"
                  # }

dataCorr1 <- dataCorr1[c(1:length(stageNum)),]
dataCorr2 <- dataCorr2[c(1:length(stageNum)),]
dataCorr3 <- dataCorr3[c(1:length(stageNum)),]

dataCorr1$stageNum <- as.numeric(as.character(stageNum))
dataCorr2$stageNum <- as.numeric(as.character(stageNum))
dataCorr3$stageNum <- as.numeric(as.character(stageNum))

dataCorr1 <- dataCorr1[dataCorr1$significant==1,] #fit the model only on significant values
dataCorr2 <- dataCorr2[dataCorr2$significant==1,] #fit the model only on significant values
dataCorr3 <- dataCorr3[dataCorr3$significant==1,] #fit the model only on significant values

sink(paste(folderAnalysis, organism, part, corMethod, "Models.txt",sep=""))

#1st degree models
lm1d1 <- lm(corValue ~ stageNum, data=dataCorr1)
lm2d1 <- lm(corValue ~ stageNum, data=dataCorr2)
lm3d1 <- lm(corValue ~ stageNum, data=dataCorr3)
 
#2nd degree models
lm1d2 <- lm(corValue ~ stageNum + I(stageNum^2), data=dataCorr1)
lm2d2 <- lm(corValue ~ stageNum + I(stageNum^2), data=dataCorr2)
lm3d2 <- lm(corValue ~ stageNum + I(stageNum^2), data=dataCorr3)

#Title for the graphs
tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA))

#Compare both models for param1
cat(paste("***** Models for", param1, "*****"))
paramName <- substitute(bold(param1), list(param1=param1))
cat("\n")
anova(lm1d1, lm1d2)
cat("\n")
if(anova(lm1d1, lm1d2)$"Pr(>F)"[[2]]<0.05){
  cat("Second-degree polynomial significantly better.\n")
  a<-summary(lm1d2)$coef[,1][[1]]
  b<-summary(lm1d2)$coef[,1][[2]]
  c<-summary(lm1d2)$coef[,1][[3]]
  #find the minimum 
  quadratic <- function(x){ eval(a) + eval(b)*x + eval(c)*x^2}
  if(c>0){
  min<-optimize(quadratic, interval = range(dataCorr1$stageNum))$minimum
  cat(paste0("On the range of the values (", param1,"), the minimum is at : ", min," ", sc, ".\n" ))    
  }else if(c<0){
    max<-optimize(quadratic, interval = range(dataCorr1$stageNum), maximum=T)$maximum
    cat(paste0("On the range of the values (", param1,"), the maximum is at : ", max," ", sc, ".\n" ))    
  }
  png(paste(folderAnalysis, organism, param1, "ModelCurve.png",sep=""))
  plot(corValue~stageNum, data=dataCorr1, ylab="Spearman's part. corr. coeff.", xlab=sc,
       xlim=range(stageNum), ylim=range(corValue))
  title (main=tit, sub=paramName)
  curve(a+b*x+c*x^2, min(dataCorr1$stageNum), max(dataCorr1$stageNum), col="red", add=T)
  dev.off()
  cat("Summary of 2nd degree polynomial regression model :\n")
  summary(lm1d2)
} else if (anova(lm1d1, lm1d2)$"Pr(>F)"[[2]]>=0.05){
  png(paste(folderAnalysis, organism, param1, "Model.png",sep=""))
  plot(corValue~stageNum, data=dataCorr1, xlab=sc, ylab=paste0("Spearman's part. corr. coeff."),
       xlim=range(stageNum), ylim=range(corValue))
  points(predict(lm1d1)~dataCorr1$stageNum, type="l", col="red", lwd=1)
  title (main=tit, sub=paramName)
  dev.off()
  cat("Second-degree polynomial not significantly better.\n")
  cat("Summary of 1st degree polynomial regression model :\n")
  summary(lm1d1)
}

#Compare both models for param2
cat(paste("\n\n***** Models for", param2, "*****"))
paramName <- substitute(bold(param2), list(param2=param2))
cat("\n")
anova(lm2d1, lm2d2)
cat("\n")

if(anova(lm2d1, lm2d2)$"Pr(>F)"[[2]]<0.05){
  cat("Second-degree polynomial significantly better.\n")
  a<-summary(lm2d2)$coef[,1][[1]]
  b<-summary(lm2d2)$coef[,1][[2]]
  c<-summary(lm2d2)$coef[,1][[3]]
  #find the minimum 
  quadratic <- function(x){ eval(a) + eval(b)*x + eval(c)*x^2}
  if(c>0){
  min<-optimize(quadratic, interval = range(dataCorr2$stageNum))$minimum
  cat(paste0("On the range of the values (", param2,"), the minimum is at : ", min," ", sc, ".\n" ))    
  }else if(c<0){
    max<-optimize(quadratic, interval = range(dataCorr2$stageNum), maximum=T)$maximum
    cat(paste0("On the range of the values (", param2,"), the maximum is at : ", max," ", sc, ".\n" ))    
  }
  png(paste(folderAnalysis, organism, param2, "ModelCurve.png",sep=""))
  plot(corValue~stageNum, data=dataCorr2, ylab=paste0("Spearman's part. corr. coeff."), xlab=sc,
       xlim=range(stageNum), ylim=range(corValue))
  title (main=tit, sub=paramName)
  curve(a+b*x+c*x^2, min(dataCorr2$stageNum), max(dataCorr2$stageNum), col="red", add=T)
  dev.off()
  
  cat("Summary of 2nd degree polynomial regression model :\n")
  summary(lm2d2)
} else if (anova(lm2d1, lm2d2)$"Pr(>F)"[[2]]>=0.05){
  png(paste(folderAnalysis, organism, param2, "Model.png",sep=""))
  plot(corValue~stageNum, data=dataCorr2, xlab=sc, ylab=paste0("Spearman's part. corr. coeff."),
       xlim=range(stageNum), ylim=range(corValue))
  points(predict(lm2d1)~dataCorr2$stageNum, type="l", col="red", lwd=1)
  title (main=tit, sub=paramName)
  dev.off()
  cat("Second-degree polynomial not significantly better.\n")
  cat("Summary of 1st degree polynomial regression model :\n")
  summary(lm2d1)
}

#Compare both models for param3
paramName <- substitute(bold(param3), list(param3=param3))
cat(paste("\n\n***** Models for", param3, "*****"))
cat("\n")
anova(lm3d1, lm3d2)
cat("\n")

if(anova(lm3d1, lm3d2)$"Pr(>F)"[[2]]<0.05){
  cat("Second-degree polynomial significantly better.\n")
  a<-summary(lm3d2)$coef[,1][[1]]
  b<-summary(lm3d2)$coef[,1][[2]]
  c<-summary(lm3d2)$coef[,1][[3]]
  #find the minimum 
  quadratic <- function(x){ eval(a) + eval(b)*x + eval(c)*x^2}
  if(c>0){
    min<-optimize(quadratic, interval = range(dataCorr3$stageNum))$minimum
    cat(paste0("On the range of the values (", param3,"), the minimum is at : ", min," ", sc, ".\n" ))    
  }else if(c<0){
    max<-optimize(quadratic, interval = range(dataCorr3$stageNum), maximum=TRUE)$maximum
    cat(paste0("On the range of the values (", param3,"), the maximum is at : ", max," ", sc, ".\n" ))    
  }
  png(paste(folderAnalysis, organism, param3, "ModelCurve.png",sep=""))
  plot(corValue~stageNum, data=dataCorr3, ylab=paste0("Spearman's part. corr. coeff."), xlab=sc,
       xlim=range(stageNum), ylim=range(corValue))
  title (main=tit, sub=paramName)
  curve(a+b*x+c*x^2, min(dataCorr3$stageNum), max(dataCorr3$stageNum), col="red", add=T)
  dev.off()
  
  cat("Summary of 2nd degree polynomial regression model :\n")
  summary(lm3d2)
} else if (anova(lm3d1, lm3d2)$"Pr(>F)"[[2]]>=0.05){
  png(paste(folderAnalysis, organism, param3, "Model.png",sep=""))
  plot(corValue~stageNum, data=dataCorr3, xlab=sc, ylab=paste0("Spearman's part. corr. coeff."),
       xlim=range(stageNum), ylim=range(corValue))
  points(predict(lm3d1)~dataCorr3$stageNum, type="l", col="red", lwd=1)
#   plot(predict(lm3d1)~dataCorr3$stageNum, type="l", col="red", lwd=1)
  title (main=tit, sub=paramName)
  dev.off()
  cat("Second-degree polynomial not significantly better.\n")
  cat("Summary of 1st degree polynomial regression model :\n")
  summary(lm3d1)
}

sink()

# -------------------------------
# One plot with all models
# -------------------------------
dataCorr <- read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""), header=T)

param1 <- "Omega"
param2 <- "Paralogs.Number"
param3 <- "Phyletic.Age"

dataCorr <- dataCorr[,c("x1", "x2", "corValue", "significant")]
dataCorr1 <- dataCorr[regexpr("stage", dataCorr$x1)>0 & regexpr(param1, dataCorr$x2)>0,]
dataCorr2 <- dataCorr[regexpr("stage", dataCorr$x1)>0 & regexpr(param2, dataCorr$x2)>0,]
dataCorr3 <- dataCorr[regexpr("stage", dataCorr$x1)>0 & regexpr(param3, dataCorr$x2)>0,]

dataCorr1 <- dataCorr1[c(1:length(stageNum)),]
dataCorr2 <- dataCorr2[c(1:length(stageNum)),]
dataCorr3 <- dataCorr3[c(1:length(stageNum)),]

stageName <- as.vector(dataCorr1$x1)
stageNum <- graphParamModel(organism)$stageNum
sc <- graphParamModel(organism)$sc
speciesName <- graphParamModel(organism)$speciesName

#Negative or positive slope ? 
ifelse(all(dataCorr1$corValue<0), sign1<-"neg", sign1<-"pos")
ifelse(all(dataCorr2$corValue<0), sign2<-"neg", sign2<-"pos")
ifelse(all(dataCorr3$corValue<0), sign3<-"neg", sign3<-"pos")
#works here so because all ...$corValue have same signe

dataCorr1$stageNum <- as.numeric(as.character(stageNum))
dataCorr2$stageNum <- as.numeric(as.character(stageNum))
dataCorr3$stageNum <- as.numeric(as.character(stageNum))

dataCorr1 <- dataCorr1[dataCorr1$significant==1,] #fit the model only on significant values
dataCorr2 <- dataCorr2[dataCorr2$significant==1,] #fit the model only on significant values
dataCorr3 <- dataCorr3[dataCorr3$significant==1,] #fit the model only on significant values

#To plot all curves on the same graph : divide by the max. in absolute terms values 
dataCorr1$relativeCor <- abs(dataCorr1$corValue/max(abs(dataCorr1$corValue)))
dataCorr2$relativeCor <- abs(dataCorr2$corValue/max(abs(dataCorr2$corValue)))
dataCorr3$relativeCor <- abs(dataCorr3$corValue/max(abs(dataCorr3$corValue)))

#1st degree models
lm1d1 <- lm(relativeCor ~ stageNum, data=dataCorr1)
lm2d1 <- lm(relativeCor ~ stageNum, data=dataCorr2)
lm3d1 <- lm(relativeCor ~ stageNum, data=dataCorr3)

#2nd degree models
lm1d2 <- lm(relativeCor ~ stageNum + I(stageNum^2), data=dataCorr1)
lm2d2 <- lm(relativeCor ~ stageNum + I(stageNum^2), data=dataCorr2)
lm3d2 <- lm(relativeCor ~ stageNum + I(stageNum^2), data=dataCorr3)

tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA))

png(paste(folderAnalysis, organism, "AllModels.png",sep=""))
plot(1, type="n", axes=T, xlab=sc, ylab="Relative values of the coefficients", 
     xlim=range(stageNum), ylim=c(0,1))
title(main = tit)

#Compare both models for param1
if(anova(lm1d1, lm1d2)$"Pr(>F)"[[2]]<0.05){
  a<-summary(lm1d2)$coef[,1][[1]]
  b<-summary(lm1d2)$coef[,1][[2]]
  c<-summary(lm1d2)$coef[,1][[3]]
  #find the minimum 
  quadratic <- function(x){ eval(a) + eval(b)*x + eval(c)*x^2}
  curve(a+b*x+c*x^2, min(dataCorr1$stageNum), max(dataCorr1$stageNum), col="red", add=T,
        lty=ifelse(sign1=="neg", 2, 1))
} else if (anova(lm1d1, lm1d2)$"Pr(>F)"[[2]]>=0.05){
  points(predict(lm1d1)~dataCorr1$stageNum, type="l", col="red", lwd=1, 
         lty=ifelse(sign1=="neg", 2, 1))
}

#Compare both models for param2
if(anova(lm2d1, lm2d2)$"Pr(>F)"[[2]]<0.05){
  a<-summary(lm2d2)$coef[,1][[1]]
  b<-summary(lm2d2)$coef[,1][[2]]
  c<-summary(lm2d2)$coef[,1][[3]]
  #find the minimum 
  quadratic <- function(x){ eval(a) + eval(b)*x + eval(c)*x^2}
  curve(a+b*x+c*x^2, min(dataCorr2$stageNum), max(dataCorr2$stageNum), col="blue", add=T,
        lty=ifelse(sign2=="neg", 2, 1))
} else if (anova(lm2d1, lm2d2)$"Pr(>F)"[[2]]>=0.05){
 points(predict(lm2d1)~dataCorr2$stageNum, type="l", col="blue", lwd=1,
        lty=ifelse(sign2=="neg", 2, 1))
}

#Compare both models for param3
if(anova(lm3d1, lm3d2)$"Pr(>F)"[[2]]<0.05){
  a<-summary(lm3d2)$coef[,1][[1]]
  b<-summary(lm3d2)$coef[,1][[2]]
  c<-summary(lm3d2)$coef[,1][[3]]
  #find the minimum 
  quadratic <- function(x){ eval(a) + eval(b)*x + eval(c)*x^2}
  curve(a+b*x+c*x^2, min(dataCorr3$stageNum), max(dataCorr3$stageNum), col="green", add=T,
        lty=ifelse(sign3=="neg", 2, 1))
} else if (anova(lm3d1, lm3d2)$"Pr(>F)"[[2]]>=0.05){
  points(predict(lm3d1)~dataCorr3$stageNum, type="l", col="green", lwd=1,
         lty=ifelse(sign3=="neg", 2, 1))
}
legend("bottomright", c(param1, param2, param3), lty=c(rep(1,3)), lwd=c(rep(2.5,3)), col=c("red", "blue", "green"))

dev.off()


# Try: models with Log axis---------------------------------------------------------------------------
            # folderAnalysis <- folderAnalysis
            # organism <- organism
            # part <- part
            # corMethod <- corMethod
            # expDataSource <- expDataSource
            # param <- "Omega"
            # dataCorr <- read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "Original.txt",sep=""), header=T)
            # dataCorr <- dataCorr[,c("x1", "x2", "corValue")]
            # dataCorr <- dataCorr[regexpr("stage", dataCorr$x1)>0 & regexpr(param, dataCorr$x2)>0,]  
            # stageName <- as.vector(dataCorr$x1)
            # #"class" defined elsewhere
            # stageNum <- graphParamContiCut(organism)$stageNum
            # sc <- graphParamContiCut(organism)$sc
            # speciesName <- graphParamContiCut(organism)$speciesName
            # 
            #             # if(organism == "Mus"){  #scale : dpc
            #             #   stageNum <- c(1, 1, 2, 3, 4, 4.5, 7.5, 8.5, 9.5, 10.5, 12, 14, 16, 18, 42) #for the log
            #             #   sc <- "days post conception"
            #             #   speciesName <- "M. musculus"
            #             # } else if(organism == "Danio"){  #scale : hours of dvpmt at 28.5°C
            #             #   stageNum <- c(1, 6, 8, 9, 10, 11.66, 16, 24, 30, 48, 96, 120, 336, 720) #try with log
            #             #   sc <- "hours of development"
            #             #   speciesName <- "D. rerio"
            #             # } else if(organism == "Drosophila" & !FPKM){ #scale : minutes after fertilization
            #             #   stageNum <- c(42.5, 110, 155, 275, 290, 450, 1200, 2190, 3630)
            #             #   sc <- "minutes after fertilization"
            #             #   speciesName <- "D. melanogaster"
            #             # } else if(organism == "Xenopus"){  #hours post-fertilization
            #             #   stageNum <- c(1.5, 2, 2.25, 2.75, 3, 5, 7, 8, 9, 17.5, 20.75, 32.5, 44.5, 66)
            #             #   sc <- "hours post-fertilization"
            #             #   speciesName <- "X. tropicalis"
            #             # } else if(organism == "Drosophila" & FPKM){  #hours post-fertilization
            #             #   stageNum <- c(seq(1, 23,2), 36.5, 60.5, 96)
            #             #   stageNum <- stageNum*60  #have it in minutes like other drosophila data
            #             #   sc <- "minutes after fertilization"
            #             #   speciesName <- "D. melanogaster"
            #             # } else if(organism == "Caeno" & FPKM){  #hours post-fertilization
            #             #   stageNum <- c(seq(1,240,30),seq(300, 720,30), 850) #try with log
            #             #   sc <- "minutes after fertilization"
            #             #   speciesName <- "C. elegans"
            #             # }
            # 
            # dataCorr <- dataCorr[c(1:length(stageNum)),]
            # dataCorr$stageNum <- as.numeric(as.character(stageNum))
            # dataCorr$stageNum <- log10(dataCorr$stageNum)
            # 
            # sink(paste(folderAnalysis, organism, part, corMethod, "ModelsLog.txt",sep=""))
            # #try the 1st degree model
            # lm1 <- lm(corValue ~ stageNum, data=dataCorr)
            # cat("For the model of the form : y= ax+b, r² = ", summary(lm1)$r.squared)
            # plot(corValue~stageNum, data=dataCorr);points(predict(lm1)~dataCorr$stageNum, type="l", col="red", lwd=1)
            # 
            # if(eval(summary(lm1)$coef[,4][2])<0.05){
            #   cat("\nThe effect of the stages is statistically significant (p-val = ", eval(summary(lm1)$coef[,4][2]),").",sep="")
            # }else if(eval(summary(lm1)$coef[,4][2])>=0.05){
            #   cat("\nThe effect of the stages is statistically not significant (p-val = ", eval(summary(lm1)$coef[,4][2]),").", sep="")
            # }
            # #try the 2nd degree model
            # lm2 <- lm(corValue ~ stageNum + I(stageNum^2), data=dataCorr)
            # cat("\n\nFor the model of the form : y= ax²+b, r² = ", summary(lm2)$r.squared)
            # 
            # if(eval(summary(lm2)$coef[,4][2])<0.05 & eval(summary(lm2)$coef[,4][3])<0.05){
            #   cat("\nThe effect of 1st degree term is statistically significant  (p-val =", eval(summary(lm2)$coef[,4][2]),")
            #       and the effect of 2nd degree is statistically significant (p-val =", eval(summary(lm2)$coef[,4][3]),").\n\n",sep="")
            #   
            # }else if(eval(summary(lm2)$coef[,4][2])>=0.05 & eval(summary(lm2)$coef[,4][3])>=0.05){
            #   cat("\nThe effects of 1st (p-val = ",eval(summary(lm2)$coef[,4][2]),") and 2nd degree terms (p-val = ", 
            #       eval(summary(lm2)$coef[,4][3]), ") \nare not statistically significant.\n\n",sep="")
            #   
            # }else if(eval(summary(lm2)$coef[,4][2])<0.05 & eval(summary(lm2)$coef[,4][3])>=0.05){
            #   cat("\nThe effect of 1st degree term is statistically significant (p-val = ", eval(summary(lm2)$coef[,4][2]),")
            #       but not the effect of the 2nd degree term (p-val = ", eval(summary(lm2)$coef[,4][3]),").\n\n", sep="")
            #   
            # }else if(eval(summary(lm2)$coef[,4][2])>=0.05 & eval(summary(lm2)$coef[,4][3])<0.05){
            #   cat("\nThe effect of 2nd degree term is statistically significant (p-val = ", eval(summary(lm2)$coef[,4][3]),")
            #       but not the effect of the 1st degree term (p-val = ", eval(summary(lm2)$coef[,4][2]),").\n\n", sep="")
            # }       
            # #compare both models :
            # anova(lm1, lm2)
            # cat("\n\n")
            # summary(lm1)
            # cat("\n\n")
            # summary(lm2)
            # 
            # if(anova(lm1, lm2)$"Pr(>F)"[[2]]<0.05){
            #   a<-summary(lm2)$coef[,1][[1]]
            #   b<-summary(lm2)$coef[,1][[2]]
            #   c<-summary(lm2)$coef[,1][[3]]
            #   #find the minimum 
            #   quadratic <- function(x){ eval(a) + eval(b)*x + eval(c)*x^2}
            #   min<-optimize(quadratic, interval = range(dataCorr$stageNum))$minimum
            #   
            #   cat("On the range of the values, the minimum is at :", min, sc )
            #   cat("\nThe nearest stage : ", stageName[which.min(abs(dataCorr$stageNum - min))])
            # }
            #   sink()
            #   
            #   png(paste(folderAnalysis, organism, part, corMethod, "ModelLog.png",sep=""))
            #   plot(corValue~stageNum, data=dataCorr, ylab="Spearman's part. corr. coeff. (Omega)", 
            #        xlab=paste("Log of",sc))
            #   points(predict(lm1)~dataCorr$stageNum, type="l", col="darkblue", lwd=1)
            #   points(predict(lm2)~dataCorr$stageNum, type="l", col="red", lwd=1)
            #   dev.off()
            #   
            #   if(anova(lm1, lm2)$"Pr(>F)"[[2]]<0.05){
            #     png(paste(folderAnalysis, organism, part, corMethod, "ModelLogCurve.png",sep=""))
            #     plot(corValue~stageNum, data=dataCorr, ylab=paste0("Spearman's part. corr. coeff. (",param,")"), 
            #          xlab=paste("Log of",sc))
            #     curve(a+b*x+c*x^2, min(dataCorr$stageNum), max(dataCorr$stageNum), col="red", add=T)
            #     dev.off()
            # }

# ---------------------------------------------------------------------------------------
# Calculate the correlations between Tau of orthologs (reference : Mus) (speciation events)
# (mouse only)
# ---------------------------------------------------------------------------------------
if(organism=="Mus"){
  #data with orthologs information
  musZebra <- read.delim("mus_zebra_ortho.txt")
  musDroso <- read.delim("mus_droso_ortho.txt")
  musXeno <- read.delim("mus_xeno_ortho.txt")
  
  #select onyl one2one orthology
  musZebra <- musZebra[regexpr("one2one",musZebra$Homology.Type)>0,] 
  musXeno <- musXeno[regexpr("one2one",musXeno$Homology.Type)>0,]  
  musDroso <- musDroso[regexpr("one2one",musDroso$Homology.Type)>0,] 
  
  o <- "Mus"
  tableMus <- read.table(paste(folderAnalysis, o, "TableStages", expDataSource,".txt",sep=""), header=TRUE)
  o <- "Danio" ; f <- "/home/user/Bureau/stage15/zebrafish/"
  tableZebra <- read.table(paste(f, o, "TableStages", expDataSource,".txt",sep=""), header=TRUE)
  o <- "Xenopus"; f <- "/home/user/Bureau/stage15/xenopus/"
  tableXeno <- read.table(paste(f, o, "TableStages", exp,".txt",sep=""), header=TRUE)
  o <- "Drosophila" ;  f <- "/home/user/Bureau/stage15/drosophila/"
  tableDroso <- read.table(paste(f, o, "TableStages", expDataSource,".txt",sep=""), header=TRUE)
  
  tableMus <- tableMus[,c("Ensembl.Gene.ID", "Tau")]
  colnames(tableMus) <- c("Ensembl.Gene.ID", "Tau.Mus")
  
  tableDroso <- tableDroso[,c("Ensembl.Gene.ID", "Tau")]
  colnames(tableDroso) <- c("Drosophila.Ensembl.Gene.ID", "Tau.Droso")
  
  tableZebra <- tableZebra[,c("Ensembl.Gene.ID", "Tau")]
  colnames(tableZebra) <- c("Zebrafish.Ensembl.Gene.ID", "Tau.Zebra")
  
  tableXeno<- tableXeno[,c("Ensembl.Gene.ID", "Tau")]
  colnames(tableXeno) <- c("Xenopus.Ensembl.Gene.ID", "Tau.Xeno")
  
  # Correlations Mus - Zebrafish
  musZebra <- merge(tableMus, musZebra, by="Ensembl.Gene.ID")
  musZebra <- merge(tableZebra, musZebra, by="Zebrafish.Ensembl.Gene.ID") 
  (ct <- cor.test(musZebra$Tau.Mus, musZebra$Tau.Zebra, method="pearson"))
  coeffMZ <- ct$estimate
  divergenceMZ <- 429.6
  numberMZ <- nrow(na.omit(musZebra))
  
  # Correlations Mus - Drosophila
  musDroso <- merge(tableMus, musDroso, by="Ensembl.Gene.ID")  
  musDroso <- merge(tableDroso, musDroso, by="Drosophila.Ensembl.Gene.ID")
  ct <- (cor.test(musDroso$Tau.Mus, musDroso$Tau.Droso, method="pearson"))
  coeffMD <- ct$estimate
  divergenceMD <- 847
  numberMD <- nrow(na.omit(musDroso))
    
  # Correlations Mus - Xenopus
  musXeno <- merge(tableMus, musXeno, by="Ensembl.Gene.ID")
  musXeno <- merge(tableXeno, musXeno, by="Xenopus.Ensembl.Gene.ID") 
  ct <- (cor.test(musXeno$Tau.Mus, musXeno$Tau.Xeno, method="pearson"))
  coeffMX <- ct$estimate
  divergenceMX <- 355.7
  numberMX <- nrow(na.omit(musXeno))
  
  coeffMM <- 1
  divergenceMM <- 1
  numberMM <- nrow(na.omit(tableMus))
  
  orthoTable <- data.frame(corValue=c(coeffMM, coeffMX, coeffMZ, coeffMD), 
                           age=c(divergenceMM, divergenceMX, divergenceMZ, divergenceMD),
                           number=c(numberMM, numberMX, numberMZ, numberMD))
  
  orthoOrg <- data.frame(Ancestor=c("X. tropicalis", "D. rerio",  "D. melanogaster"),
                         Age=c(divergenceMX, divergenceMZ, divergenceMD))
  
  # ---------------------------------------------------------------------------------------
  # Correlation between Tau of Paralogs (within species paralogs, duplication events)
  # ---------------------------------------------------------------------------------------
  
  musPara <- read.delim("MusParalogsYoungestCoupleEnsV75.txt")  #"Ensembl.Gene.ID", "Paralog.E.G.ID", "Ancestor"  
  musAnces <- read.delim("MusAncestorDifference.txt")           #"Organisms", "Age"
  colnames(musAnces) <-c("Ancestor", "Age")
  para <- merge(musPara, musAnces, by="Ancestor") 
  
  #For the label of the x-axis (needed after): 
  time <- rbind(musAnces, orthoOrg)
  time <- time[order(time$Age),]
  time$color <- ifelse(time$Ancestor %in% musAnces$Ancestor, "grey", "red") 
  
  tableMus <- read.table(paste(folderAnalysis, organism, "TableStages", expDataSource,".txt",sep=""), header=TRUE) 
  tableRef <- merge(tableMus, para, by="Ensembl.Gene.ID") 
  tableRef <- tableRef[,c("Tau", "Ensembl.Gene.ID", "Paralog.Ensembl.Gene.ID", "Age")]  #Tau for the Ensembl.Gene.ID
  
  colnames(tableRef) <- c("Tau.Para1", "Ensembl.Gene.ID.Para1", "Ensembl.Gene.ID", "Age") #Tau for the Paralog.E.G.ID
  tablePara <- merge(tableRef, tableMus, by="Ensembl.Gene.ID") 
  tablePara <- tablePara[,c("Tau.Para1", "Ensembl.Gene.ID.Para1", "Ensembl.Gene.ID", "Tau", "Age")]
  colnames(tablePara) <- c("Tau.Para1", "Ensembl.Gene.ID.Para1", "Ensembl.Gene.ID.Para2", "Tau.Para2", "Age")
  
  tableParSave <- tablePara 
  tablePara <- na.omit(tablePara) 
  
  paraSum <- plyr::count(tablePara, "Age")  #can be masked by dplyr::count !
  
  coeffPara <- list()
  agePara <- list()
  tablePara$Age <- as.factor(tablePara$Age)
  for(i in levels(tablePara$Age)){
    x <- subset(tablePara, tablePara$Age==i)$Tau.Para1
    y <- subset(tablePara, tablePara$Age==i)$Tau.Para2
    if(length(x)>2 & length(y)>2){
      corValue <- cor.test(x, y, method="pearson")$estimate
      coeffPara <- append(coeffPara, corValue)
      agePara <- append(agePara, i)
    }
  }
  paraTable <- data.frame(corValue=unlist(coeffPara), Age=unlist(agePara))
  
  paraTable$Age <- as.numeric(as.character(paraTable$Age))        #corValue, Age
  paraTable[1,2] <- 1  #replace age=0 for mus (logarithmic scale)
  
  paraTable <- merge(paraTable, paraSum, by="Age")
  paraTable <- paraTable[,c("corValue", "Age", "freq")]
  colnames(paraTable) <- c("corValue", "age", "number")
  
  # ---------------------------------------------------------------------------------------
  # Plot orthologs and paralogs of mouse
  # ---------------------------------------------------------------------------------------
  time$Age[1] <- 1
  ggplot() +
  theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())+
  geom_point(data=orthoTable, colour="red", aes(x=age, y=corValue,group=1, size=number, guide="none"), 
             show_guide=F) +
  scale_size_continuous(range=c(2,10))+
  geom_smooth(data=orthoTable, colour="red", aes(x=age, y=corValue), method="lm",se=F)+
   geom_line(data=orthoTable, colour="black", aes(x=age, y=corValue))+
   geom_point(data=paraTable, aes(x=age, y=corValue, group=1, size=(number)), colour="blue", show_guide=F)+
  geom_line(data=paraTable, colour="black", aes(x=age, y=corValue))+
  geom_smooth(data=paraTable, colour="blue", aes(x=age, y=corValue), method="lm", se=F)+
  scale_y_continuous("Pearson correlation coefficient")+
  # scale_x_continuous("Age (million years)",breaks=time$Age, labels=time$Ancestor)+
  scale_x_log10("Age (million years)",breaks=time$Age, labels=time$Ancestor)+
  theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
          axis.text.y = element_text(colour="black"),
          axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=time$color))
  ggsave(filename=paste0(folderAnalysis,"OrthoParaCorrrelation.png"))
}

# ---------------------------------------------------------------------------------------
# Take only more or less specific genes to calculate correlations
# ---------------------------------------------------------------------------------------
### Load the data
dataOrg <- read.table(paste(folderAnalysis, organism, "TableStages", expDataSource,".txt",sep=""), header=TRUE)
dataOrg <- dataOrg[dataOrg$Tau>0.3,]
if (partial){
  cat("Partial correlation was performed.",sep="")
  part <- "Partial"
} else {
  cat("Normal correlation was performed.",sep="")
  part <- "Normal"
}

###Select only needed columns
if(organism=="Caeno"){
  if(phylAge){
    dataOrg <- dataOrg[,c("Omega.0", "CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", "Phyletic.Age", devStageOutNames)]
    parameterNames <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number", "Phyletic.Age")
  } else if (!phylAge){
    ### we don't have phyletic age for the moment !!! and no Omega
    dataOrg <- dataOrg[,c("Omega.0","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", devStageOutNames)]
    parameterNames <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number")  
  }  
} else {
  if(phylAge){
    dataOrg <- dataOrg[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", "Phyletic.Age", devStageOutNames)]
    parameterNames <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number", "Phyletic.Age")
  } else if (!phylAge){
    ### we don't have phyletic age for the moment !!! and no Omega
    dataOrg <- dataOrg[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", devStageOutNames)]
    parameterNames <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number")  
  }
}

colnames(dataOrg) <- c(parameterNames, devStageOutNames)
geneDataOrg <- na.omit(dataOrg)


### Normalization of the data (function defined in "scriptFunctions.R")
geneDataOrg<-normalizeData2(geneDataOrg)
  
### Calculation correlation
variableNames <- devStageOutNames
x <- data.frame(x1=NULL,x2=NULL,corValue=NULL,pValue=NULL,significant=NULL)
for(j in variableNames){
  #j is the name of variable for which the correlation is calculated
  variablesToUse <- parameterNames #Names of other variables
  t = length(variablesToUse)
  for(n in c(1:t)){
    j2 <- variablesToUse[n]
    variablesToUse2 <- variablesToUse[variablesToUse != j2]
    
    if(partial==TRUE){
      fmodel <-"geneDataOrg$"
      fmodel <- paste(fmodel,j,"~",sep="")
      fmodel2 <-"geneDataOrg$"
      fmodel2 <- paste(fmodel2,j2,"~",sep="")
      for(i in variablesToUse2){
        fmodel <- paste(fmodel,"geneDataOrg$",i,"+",sep="")
        fmodel2 <- paste(fmodel2,"geneDataOrg$",i,"+",sep="")
      }
      fmodel <- substr(fmodel, 1, nchar(fmodel)-1) #delet last character "+"
      fmodel2 <- substr(fmodel2, 1, nchar(fmodel2)-1) #delet last character "+"
      fmx <- glm(fmodel, na.action = na.exclude)
      fmy <- glm(fmodel2, na.action = na.exclude)
      xres <- resid(fmx)
      yres <- resid(fmy)
      ct <- cor.test(xres, yres, method=corMethod)      
    } else {
      ct <- cor.test(geneDataOrg[,j], geneDataOrg[,j2], method=corMethod)
    }
    s <- ct$estimate
    coeff <- ct$p.value
    x <- rbind(x, data.frame(x1=j,x2=j2,corValue=s,pValue=coeff))
  }
}

saveX3a <- x
variableNames <- parameterNames
for(j in variableNames) {
  #j is the name of variable for which the correlation is calculated
  variablesToUse <- variableNames[variableNames != j] #Names of other variables
  t = length(variablesToUse)
  
  for(n in c(1:t)){
    j2 <- variablesToUse[n]
    variablesToUse2 <- variablesToUse[variablesToUse != j2]
    
    if(partial==TRUE){
      fmodel <-"geneDataOrg$"
      fmodel <- paste(fmodel,j,"~",sep="")
      fmodel2 <-"geneDataOrg$"
      fmodel2 <- paste(fmodel2,j2,"~",sep="")
      
      for(i in variablesToUse2){
        fmodel <- paste(fmodel,"geneDataOrg$",i,"+",sep="")
        fmodel2 <- paste(fmodel2,"geneDataOrg$",i,"+",sep="")
      }
      
      fmodel <- substr(fmodel, 1, nchar(fmodel)-1) #delet last character "+"
      fmodel2 <- substr(fmodel2, 1, nchar(fmodel2)-1) #delet last character "+"
      fmx <- glm(fmodel, na.action = na.exclude)
      fmy <- glm(fmodel2, na.action = na.exclude)
      xres <- resid(fmx)
      yres <- resid(fmy)
      ct <- cor.test(xres, yres, method=corMethod)
    } else {
      ct <- cor.test(geneDataOrg[,j], geneDataOrg[,j2], method=corMethod)
    }
    s <- ct$estimate
    coeff <- ct$p.value
    x <- rbind(x, data.frame(x1=j,x2=j2,corValue=s,pValue=coeff))
  }
}

#For the test of significance -> adjust the p-value for multiple comparisons
# chosen method : Benjamini & Hochberg ; threshold : 0.05 for the adjusted p-val
x$pValAdjust <- p.adjust(x$pValue,method="BH")
x$significant <- as.integer(x$pValAdjust < 0.05)

saveX3b <- x

row.names(x) <- c(1:nrow(x))
x$x2 <- factor(x$x2,levels=levels(x$x1))
sameCorRowNumbers <- vector("numeric")

for(i in rownames(x)){
  sameCor <- x[with(x,x$x2 == x[i,]$x1 & x$x1 == x[i,]$x2),]   #remove correlation with itself
  if(nrow(sameCor)>0){
    if(as.integer(rownames(sameCor))>as.integer(i)){
      sameCorRowNumbers <- append(sameCorRowNumbers, as.integer(rownames(sameCor)))  #remove the 2n of a1xa2 [...] a2xa1
    }
  }
}
  
x <- x[-sameCorRowNumbers,]
saveX4 <- x    # 99rows

cat("\n Correlation table",sep="")

cat("\n Result of the correlation is saved in \"CorStages_Original\".",sep="")

write.table(x,file=paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "specific.txt",sep=""),row.names = FALSE,quote = FALSE)

#####*****Take only broadly (accross developmental stages) expressed genes*****#####
dataOrg <- read.table(paste(folderAnalysis, organism, "TableStages", expDataSource,".txt",sep=""), header=TRUE)
dataOrg <- dataOrg[dataOrg$Tau <= 0.3,]
if (partial){
  cat("Partial correlation was performed.",sep="")
  part <- "Partial"
} else {
  cat("Normal correlation was performed.",sep="")
  part <- "Normal"
}

###Select only needed columns
if(organism=="Caeno"){
  if(phylAge){
    dataOrg <- dataOrg[,c("Omega.0", "CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", "Phyletic.Age", devStageOutNames)]
    parameterNames <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number", "Phyletic.Age")
  } else if (!phylAge){
    ### we don't have phyletic age for the moment !!! and no Omega
    dataOrg <- dataOrg[,c("Omega.0","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", devStageOutNames)]
    parameterNames <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number")  
  }  
} else {
  if(phylAge){
    dataOrg <- dataOrg[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", "Phyletic.Age", devStageOutNames)]
    parameterNames <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number", "Phyletic.Age")
  } else if (!phylAge){
    ### we don't have phyletic age for the moment !!! and no Omega
    dataOrg <- dataOrg[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", devStageOutNames)]
    parameterNames <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number")  
  }
}

colnames(dataOrg) <- c(parameterNames, devStageOutNames)
geneDataOrg <- na.omit(dataOrg)

##### Normalization of the data (function defined in "scriptFunctions.R")
geneDataOrg<-normalizeData2(geneDataOrg)

##### Calculation correlation
variableNames <- devStageOutNames
x <- data.frame(x1=NULL,x2=NULL,corValue=NULL,pValue=NULL,significant=NULL)
for(j in variableNames){
  #j is the name of variable for which the correlation is calculated
  variablesToUse <- parameterNames #Names of other variables
  t = length(variablesToUse)
  for(n in c(1:t)){
    j2 <- variablesToUse[n]
    variablesToUse2 <- variablesToUse[variablesToUse != j2]
    
    if(partial==TRUE){
      fmodel <-"geneDataOrg$"
      fmodel <- paste(fmodel,j,"~",sep="")
      fmodel2 <-"geneDataOrg$"
      fmodel2 <- paste(fmodel2,j2,"~",sep="")
      for(i in variablesToUse2){
        fmodel <- paste(fmodel,"geneDataOrg$",i,"+",sep="")
        fmodel2 <- paste(fmodel2,"geneDataOrg$",i,"+",sep="")
      }
      fmodel <- substr(fmodel, 1, nchar(fmodel)-1) #delet last character "+"
      fmodel2 <- substr(fmodel2, 1, nchar(fmodel2)-1) #delet last character "+"
      fmx <- glm(fmodel, na.action = na.exclude)
      fmy <- glm(fmodel2, na.action = na.exclude)
      xres <- resid(fmx)
      yres <- resid(fmy)
      ct <- cor.test(xres, yres, method=corMethod)      
    } else {
      ct <- cor.test(geneDataOrg[,j], geneDataOrg[,j2], method=corMethod)
    }
    s <- ct$estimate
    coeff <- ct$p.value
    x <- rbind(x, data.frame(x1=j,x2=j2,corValue=s,pValue=coeff))
  }
}

saveX3a <- x
variableNames <- parameterNames
for(j in variableNames) {
  #j is the name of variable for which the correlation is calculated
  variablesToUse <- variableNames[variableNames != j] #Names of other variables
  t = length(variablesToUse)
  
  for(n in c(1:t)){
    j2 <- variablesToUse[n]
    variablesToUse2 <- variablesToUse[variablesToUse != j2]
    
    if(partial==TRUE){
      fmodel <-"geneDataOrg$"
      fmodel <- paste(fmodel,j,"~",sep="")
      fmodel2 <-"geneDataOrg$"
      fmodel2 <- paste(fmodel2,j2,"~",sep="")
      
      for(i in variablesToUse2){
        fmodel <- paste(fmodel,"geneDataOrg$",i,"+",sep="")
        fmodel2 <- paste(fmodel2,"geneDataOrg$",i,"+",sep="")
      }
      
      fmodel <- substr(fmodel, 1, nchar(fmodel)-1) #delet last character "+"
      fmodel2 <- substr(fmodel2, 1, nchar(fmodel2)-1) #delet last character "+"
      fmx <- glm(fmodel, na.action = na.exclude)
      fmy <- glm(fmodel2, na.action = na.exclude)
      xres <- resid(fmx)
      yres <- resid(fmy)
      ct <- cor.test(xres, yres, method=corMethod)
    } else {
      ct <- cor.test(geneDataOrg[,j], geneDataOrg[,j2], method=corMethod)
    }
    s <- ct$estimate
    coeff <- ct$p.value
    x <- rbind(x, data.frame(x1=j,x2=j2,corValue=s,pValue=coeff))
  }
}

#For the test of significance -> adjust the p-value for multiple comparisons
# chosen method : Benjamini & Hochberg ; threshold : 0.05 for the adjusted p-val
x$pValAdjust <- p.adjust(x$pValue,method="BH")
x$significant <- as.integer(x$pValAdjust < 0.05)

saveX3b <- x

row.names(x) <- c(1:nrow(x))
x$x2 <- factor(x$x2,levels=levels(x$x1))
sameCorRowNumbers <- vector("numeric")

for(i in rownames(x)){
  sameCor <- x[with(x,x$x2 == x[i,]$x1 & x$x1 == x[i,]$x2),]   #remove correlation with itself
  if(nrow(sameCor)>0){
    if(as.integer(rownames(sameCor))>as.integer(i)){
      sameCorRowNumbers <- append(sameCorRowNumbers, as.integer(rownames(sameCor)))  #remove the 2n of a1xa2 [...] a2xa1
    }
  }
}

x <- x[-sameCorRowNumbers,]
saveX4 <- x    # 99rows

cat("\n Correlation table",sep="")

cat("\n Result of the correlation is saved in \"CorStages_Original\".",sep="")

write.table(x,file=paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "broad.txt",sep=""),row.names = FALSE,quote = FALSE)

#####*****Plot correlations with broad/specififc genes*****#####
specificCorr <-read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "specific.txt",sep=""),header=T)
broadCorr <-read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "broad.txt",sep=""),header=T)

#Function for plotting defined in "scriptFunctions.R"
plotByStages(broadCorr, "Omega", "broad")
plotByStages(specificCorr, "Omega", "specific")

# ---------------------------------------------------------------------------------------
# Correlation between tissue- and stage-specificity 
# (mouse only, tissue-specificity from Kryuchkova-Mostacci & Robinson-Rechavi 2015)
# ---------------------------------------------------------------------------------------
if(organism=="Mus"){
  tissuTau <- read.table("~/Bureau/stage15/article/DevelopmentMouseMusExpressionENCODE.txt", header=T, sep=" ")
  devStageNames <- colnames(tissuTau)[regexpr("Averaged", colnames(tissuTau))>0]
  tissuTau$Tau <- apply(tissuTau[,c(devStageNames,paste("Max.Expression",sep=""))],1,ftau)
  tissuTau <- tissuTau[,c("Ensembl.Gene.ID", "Tau")]
  colnames(tissuTau) <- c("Ensembl.Gene.ID", "Tau.Tissues") 
  
  stageTau <- read.table("~/Bureau/stage15/mouse/analysis_14.08_PA/MusExpressionBgee.txt", header=T, sep=" ")
  stageTau <- stageTau[,c("Ensembl.Gene.ID", "Tau")]
  colnames(stageTau) <- c("Ensembl.Gene.ID", "Tau.Stages")  
  
  compareTau <- merge(tissuTau, stageTau, by="Ensembl.Gene.ID")
  
  ct <- cor.test(compareTau$Tau.Stages, compareTau$Tau.Tissues, method="pearson")
  
  sink("correlationStagesTissuesSpecificity.txt")
  cat("\nPearson correlation coefficient between tissues and stages specificity :")
  cat("\n", ct$estimate, ", p-value :", ct$p.value)
  if(ct$p.value<0.05){
    cat("\nThe correlation is statistically significant.\n\n")
  }else if(ct$p.value>=0.05){
    cat("\nThe correlation is not statistically significant.\n\n")
  }
  ct
  sink()
}
# ---------------------------------------------------------------------------------------
# Compare % of identity with orthologous genes across development
# (pairs: mouse and other organisms)
# ---------------------------------------------------------------------------------------
#With Drosophila orthologous genes
orthoMD <- read.delim("~/Bureau/stage15/orthology/mus_droso.txt") 
expression <- read.table("~/Bureau/stage15/mouse/analysis_14.08_PA/MusExpressionBgee.txt", header=T, sep=" ")
orthoMD <- orthoMD[regexpr("one2one",orthoMD$Homology.Type)>0,] 
devStageNames <- colnames(expression)[regexpr("stage", colnames(expression))>0]
expression <- expression[,c("Ensembl.Gene.ID", devStageNames, "Tau")] 
expOrtho <- merge(orthoMD, expression, by="Ensembl.Gene.ID") 
ortholog <- "D. melanogaster"

sink(paste0(folderAnalysis, "orthology/corrTauX.OrthologyWith",ortholog,".txt"))
cat("Correlation test between Tau and % of identity with", ortholog, "orthologous genes :\n")
cor.test(expOrtho$Tau, expOrtho$X..Identity.with.respect.to.query.gene, method="pearson")
cat("\nNumber of orthologous genes used : ", nrow(expOrtho))
sink()

corrCoeff <- vector()
p.values <- vector()
plotOrthologous(ortholog) #function defined in "scriptFunctions.R"

#With Xenopus orthologous genes
orthoMD <- read.delim("~/Bureau/stage15/orthology/mus_xenopus.txt") 
expression <- read.table("~/Bureau/stage15/mouse/analysis_14.08_PA/MusExpressionBgee.txt", header=T, sep=" ")
orthoMD <- orthoMD[regexpr("one2one",orthoMD$Homology.Type)>0,] 
devStageNames <- colnames(expression)[regexpr("stage", colnames(expression))>0]
expression <- expression[,c("Ensembl.Gene.ID", devStageNames, "Tau")] 
expOrtho <- merge(orthoMD, expression, by="Ensembl.Gene.ID") 
ortholog <- "X. tropicalis"

sink(paste0(folderAnalysis, "orthology/corrTauX.OrthologyWith",ortholog,".txt"))
cat("Correlation test between Tau and % of identity with", ortholog, "orthologous genes :\n")
cor.test(expOrtho$Tau, expOrtho$X..Identity.with.respect.to.query.gene, method="pearson")
cat("\nNumber of orthologous genes used : ", nrow(expOrtho))
sink()
corrCoeff <- vector()
p.values <- vector()
plotOrthologous(ortholog)

#With Danio orthologous genes
orthoMD <- read.delim("~/Bureau/stage15/orthology/mus_zebra.txt")   
expression <- read.table("~/Bureau/stage15/mouse/analysis_14.08_PA/MusExpressionBgee.txt", header=T, sep=" ")
devStageNames <- colnames(expression)[regexpr("stage", colnames(expression))>0]
orthoMD <- orthoMD[regexpr("one2one",orthoMD$Homology.Type)>0,] 
expression <- expression[,c("Ensembl.Gene.ID", devStageNames, "Tau")] 
expOrtho <- merge(orthoMD, expression, by="Ensembl.Gene.ID") 
ortholog <- "D. rerio"

sink(paste0(folderAnalysis, "orthology/corrTauX.OrthologyWith",ortholog,".txt"))
cat("Correlation test between Tau and % of identity with", ortholog, "orthologous genes :\n")
cor.test(expOrtho$Tau, expOrtho$X..Identity.with.respect.to.query.gene, method="pearson")
cat("\nNumber of orthologous genes used : ", nrow(expOrtho))
sink()
corrCoeff <- vector()
p.values <- vector()
plotOrthologous(ortholog)

#With Caeno orthologous genes
orthoMD <- read.delim("~/Bureau/stage15/orthology/mus_caeno.txt")
expression <- read.table("~/Bureau/stage15/mouse/analysis_14.08_PA/MusExpressionBgee.txt", header=T, sep=" ")
devStageNames <- colnames(expression)[regexpr("stage", colnames(expression))>0]
orthoMD <- orthoMD[regexpr("one2one",orthoMD$Homology.Type)>0,] 
expression <- expression[,c("Ensembl.Gene.ID", devStageNames, "Tau")] 
expOrtho <- merge(orthoMD, expression, by="Ensembl.Gene.ID") 
ortholog <- "C. elegans"

sink(paste0(folderAnalysis, "orthology/corrTauX.OrthologyWith",ortholog,".txt"))
cat("Correlation test between Tau and % of identity with", ortholog, "orthologous genes :\n")
cor.test(expOrtho$Tau, expOrtho$X..Identity.with.respect.to.query.gene, method="pearson")
cat("\nNumber of orthologous genes used : ", nrow(expOrtho))
sink()
corrCoeff <- vector()
p.values <- vector()
plotOrthologous(ortholog)

# ---------------------------------------------------------------------------------------
# Correlation between the level of expression at each stage and the stage specificity
# ---------------------------------------------------------------------------------------
#Parameters to set if this part of the code is run separately
folderAnalysis <-folderAnalysis 
organism <-organism 
FPKM <- FPKM
expDataSource <-expDataSource 
expData <- read.table(paste(folderAnalysis, organism, "Expression", expDataSource,".txt",sep=""), header=T)
stageNames <- colnames(expData)[regexpr("stage", colnames(expData))>0]

c.Val <- vector()
p.Val <- vector()
for(i in stageNames){
  ct <- cor.test(expData[,i], expData$Tau, method="spearman")
  coeff <- ct$estimate
  p.val <- ct$p.value
  c.Val <- append(coeff, c.Val)
  p.Val <- append(p.val, p.Val)
}

tauEvo <- data.frame(stage=as.factor(stageNames), corrVal=c.Val, pVal=p.Val)
tauEvo$significant <- ifelse(tauEvo$pVal<0.05, 1, 0)


colLabel <- graphParamAll(organism)$colLabel
symbRNA <- graphParamAll(organism)$symbRNA
speciesName <- graphParamAll(organism)$speciesName

tauEvo <- tauEvo[c(1:length(colLabel)),]

tauEvo$stage <- factor(as.character(tauEvo$stage),levels=tauEvo$stage)
tauEvo$significant <- factor(as.character(tauEvo$significant),levels=tauEvo$significant)
tit <- substitute(italic(speciesName^{symbRNA}), list(speciesName=speciesName, symbRNA=symbRNA))

x11()
ggplot(tauEvo, aes(x=factor(stage), y=corrVal,group=1,colour=significant)) +
  geom_point(size=4) + geom_line(colour="black")+
  scale_colour_manual(values=c("0"="chocolate","1"="green"),guide="none")+
  scale_y_continuous("Spearmans' (partial) correlation coefficient")+
  scale_x_discrete("Stage specificity of expression")+
  ggtitle(tit)+
  theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
  theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())

ggsave(filename=paste0(folderAnalysis,organism,"Tau.png"))
dev.off()


# ---------------------------------------------------------------------------------------
# Plot connectivity
# ---------------------------------------------------------------------------------------
#For zebrafish : no connectivity for the genes for which we have values for the other parameters
#Parameters to set if this part of the code is run separately
expDataSource <- expDataSource
organism <- organism 
folderAnalysis <- folderAnalysis 
partial <- partial 
phylAge <-phylAge 
corMethod <- corMethod 
FPKM <- FPKM

dataOrg <- read.table(paste(folderAnalysis, organism, "TableStages", expDataSource,".txt",sep=""), header=T)
dataOrg <- dataOrg[dataOrg$Max.Expression>0.00015,]
devStageOutNames <- colnames(dataOrg)[regexpr("stage", colnames(dataOrg))>0]

##### All parameters
###Leave only needed columns

if(organism=="Caeno"){
  if(phylAge){
    dataOrg <- dataOrg[,c("Omega.0", "CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", "Phyletic.Age", "Connectivity", devStageOutNames)]
    parameterNames <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number", "Phyletic.Age", "Connectivity")
  } else if (!phylAge){
    ### we don't have phyletic age for the moment !!! and no Omega
    dataOrg <- dataOrg[,c("Omega.0","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", "Connectivity", devStageOutNames)]
    parameterNames <- c("Omega", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number", "Connectivity")  
  }  
} else {
  if(phylAge){
    dataOrg <- dataOrg[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", "Phyletic.Age", "Connectivity", devStageOutNames)]
    parameterNames <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number", "Phyletic.Age", "Connectivity")
  } else if (!phylAge){
    ### we don't have phyletic age for the moment !!! and no Omega
    dataOrg <- dataOrg[,c("Omega.0", "LRT", "P.1","CDS.Length", "Intron.Length", "Intron.Number", 
                          "X..GC.content", "Paralogs.Number", "Connectivity", devStageOutNames)]
    parameterNames <- c("Omega", "LRT", "P.1", "CDS.Length", "Intron.Length", "Intron.Number", 
                        "X..GC.content", "Paralogs.Number", "Connectivity")  
  }
}

colnames(dataOrg) <- c(parameterNames, devStageOutNames)
geneDataOrg <- na.omit(dataOrg)

##### Normalization of the data (function defined elsewhere)
geneDataOrg<-normalizeData2(geneDataOrg)

##### Calculation correlation
variableNames <- devStageOutNames
x <- data.frame(x1=NULL,x2=NULL,corValue=NULL,pValue=NULL,significant=NULL)
for(j in variableNames){
  #j is the name of variable for which the correlation is calculated
  variablesToUse <- parameterNames #Names of other variables
  t = length(variablesToUse)
  for(n in c(1:t)){
    j2 <- variablesToUse[n]
    variablesToUse2 <- variablesToUse[variablesToUse != j2]
    
    if(partial==TRUE){
      fmodel <-"geneDataOrg$"
      fmodel <- paste(fmodel,j,"~",sep="")
      fmodel2 <-"geneDataOrg$"
      fmodel2 <- paste(fmodel2,j2,"~",sep="")
      for(i in variablesToUse2){
        fmodel <- paste(fmodel,"geneDataOrg$",i,"+",sep="")
        fmodel2 <- paste(fmodel2,"geneDataOrg$",i,"+",sep="")
      }
      fmodel <- substr(fmodel, 1, nchar(fmodel)-1) #delet last character "+"
      fmodel2 <- substr(fmodel2, 1, nchar(fmodel2)-1) #delet last character "+"
      fmx <- glm(fmodel, na.action = na.exclude)
      fmy <- glm(fmodel2, na.action = na.exclude)
      xres <- resid(fmx)
      yres <- resid(fmy)
      ct <- cor.test(xres, yres, method=corMethod)      
    } else {
      ct <- cor.test(geneDataOrg[,j], geneDataOrg[,j2], method=corMethod)
    }
    s <- ct$estimate
    coeff <- ct$p.value
    x <- rbind(x, data.frame(x1=j,x2=j2,corValue=s,pValue=coeff))
  }
}

saveX3a <- x
variableNames <- parameterNames
for(j in variableNames) {
  #j is the name of variable for which the correlation is calculated
  variablesToUse <- variableNames[variableNames != j] #Names of other variables
  t = length(variablesToUse)
  
  for(n in c(1:t)){
    j2 <- variablesToUse[n]
    variablesToUse2 <- variablesToUse[variablesToUse != j2]
    
    if(partial==TRUE){
      fmodel <-"geneDataOrg$"
      fmodel <- paste(fmodel,j,"~",sep="")
      fmodel2 <-"geneDataOrg$"
      fmodel2 <- paste(fmodel2,j2,"~",sep="")
      
      for(i in variablesToUse2){
        fmodel <- paste(fmodel,"geneDataOrg$",i,"+",sep="")
        fmodel2 <- paste(fmodel2,"geneDataOrg$",i,"+",sep="")
      }
      
      fmodel <- substr(fmodel, 1, nchar(fmodel)-1) #delet last character "+"
      fmodel2 <- substr(fmodel2, 1, nchar(fmodel2)-1) #delet last character "+"
      fmx <- glm(fmodel, na.action = na.exclude)
      fmy <- glm(fmodel2, na.action = na.exclude)
      xres <- resid(fmx)
      yres <- resid(fmy)
      ct <- cor.test(xres, yres, method=corMethod)
    } else {
      ct <- cor.test(geneDataOrg[,j], geneDataOrg[,j2], method=corMethod)
    }
    s <- ct$estimate
    coeff <- ct$p.value
    x <- rbind(x, data.frame(x1=j,x2=j2,corValue=s,pValue=coeff))
  }
}

#For the test of significance -> adjust the p-value for multiple comparisons
# chosen method : Benjamini & Hochberg ; threshold : 0.05 for the adjusted p-val
x$pValAdjust <- p.adjust(x$pValue,method="BH")
x$significant <- as.integer(x$pValAdjust < 0.05)

saveX3b <- x

row.names(x) <- c(1:nrow(x))
x$x2 <- factor(x$x2,levels=levels(x$x1))
sameCorRowNumbers <- vector("numeric")

for(i in rownames(x)){
  sameCor <- x[with(x,x$x2 == x[i,]$x1 & x$x1 == x[i,]$x2),]   #remove correlation with itself
  if(nrow(sameCor)>0){
    if(as.integer(rownames(sameCor))>as.integer(i)){
      sameCorRowNumbers <- append(sameCorRowNumbers, as.integer(rownames(sameCor)))  #remove the 2n of a1xa2 [...] a2xa1
    }
  }
}

x <- x[-sameCorRowNumbers,]
saveX4 <- x    # 99rows
write.table(x,file=paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "OriginalWithConnect.txt",sep=""),row.names = FALSE,quote = FALSE)

### Plot variation correlation connectivity
data <- x
x <- "Connectivity"
#data<-read.table(paste(folderAnalysis, organism, part, corMethod, "CorStages", expDataSource, "OriginalWithConnect.txt",sep=""),header=T)

data <- data[which(regexpr("stage", data$x1)>0 & regexpr(x, data$x2)>0),] 
data$significant <- as.factor(data$significant) #need as factor for the colour below
speciesName <- graphParamAll(organism)$speciesName
symbRNA <- graphParamAll(organism)$symbRNA
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
  scale_x_discrete(paste0(x))+
  ggtitle(tit)+
  theme(axis.title.y = element_text(face="bold", colour="#990000", size=10),
        axis.text.y = element_text(colour="black"),
        axis.title.x = element_text(face="bold", colour="#CD5B45", size=20),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=10, colour=colLabel))+
  theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank())

ggsave(filename=paste0(folderAnalysis,x,organism,".png"))


# ---------------------------------------------------------------------------------------
# Draw expression according to phyletic age
# ---------------------------------------------------------------------------------------
#Parameters to set if this part of the code is run separately 
organism <- organism 
FPKM <- FPKM
expDataSource <- expDataSource 
values <- values 
twoDataSets <- twoDataSets

if (expDataSource == "Bgee" ) {
  values <- "Log.of.normalized.signal.intensity"
}else if (expDataSource == "BgeeRNA" ) {
  values <- "RPKM"
}

### Input data
if(organism=="Mus"){
  orgPhyleticage <- read.table("ageMus.txt", sep="\t")
  twoDataSets <- TRUE
} else if(organism=="Danio"){
  # with the phyletic age retrieved from Ensembl
  orgPhyleticage <- read.table("ageDanio.txt", sep="\t")
  ### Gene expression 
  #orgExpression <- read.table(paste("~/Expression information"), sep="\t", header=TRUE)
  orgExpression <- read.delim("danioETABM33.tsv") #  (Affymetrix data from Bgee)
} else if(organism=="Drosophila"){
  # with the phyletic age retrieved from Ensembl
  orgPhyleticage <- read.table("ageDroso.txt", sep="\t")
  if(!FPKM){
    twoDataSets <- TRUE
  }else if (FPKM){
    orgExpression <- read.table("droso_FPKMs.txt",header=T)
    colnames(orgExpression)[1] <- "Ensembl.Gene.ID"
    orgExpression[,-1] <- as.numeric(as.character(unlist((orgExpression[,-1]))))
    orgExpression$Mean.1d <- rowMeans(subset(orgExpression, select=c(Male.1d, Female.1d), na.rm=TRUE))
    orgExpression$Mean.5d <- rowMeans(subset(orgExpression, select=c(Male.5d, Female.5d), na.rm=TRUE))                                             
    orgExpression$Mean.30d <- rowMeans(subset(orgExpression, select=c(Male.30d, Female.30d), na.rm=TRUE))
    orgExpression <- orgExpression[,c(colnames(orgExpression)[regexpr("ale", colnames(orgExpression))<0])]
  }
} else if(organism=="Xenopus"){
 # with the phyletic age retrieved from Ensembl
  orgPhyleticage <- read.table("ageXeno.txt", sep="\t")
  ### Gene expression 
  #orgExpression <- read.table(paste("~/Expression information"), sep="\t", header=TRUE)
  orgExpression <- read.delim("xenoGSE37452.tsv") #  (RNA-seq from Bgee)
  orgExpression <- orgExpression[,c("Gene.ID","Stage.name","RPKM")]   
} else if(organism=="Caeno"){
  # with the phyletic age retrieved from Ensembl
  orgPhyleticage <- read.table("ageCaeno.txt", sep="\t")
  ### Gene expression 
  if(FPKM){
    orgExpression <- read.delim("caeno_FPKMs.txt")  
    colnames(orgExpression)[1] <- "Ensembl.Transcript.ID"
    devStageNames <- colnames(orgExpression)[-1]
    orgExpression[,-1] <- as.numeric(as.character(unlist((orgExpression[,-1]))))
    geneTransc <- read.delim("CaenoGenesTransc.txt") #Ens81, no filters, unique results only  
    orgExpression <- merge(orgExpression, geneTransc, by="Ensembl.Transcript.ID")
    orgExpression <- orgExpression[,c("Ensembl.Gene.ID", devStageNames)]
  }
} 

### Reshape the datasets
colnames(orgPhyleticage) <- c("Ensembl.Gene.ID", "Taxon")
orgPhyleticage <- merge(orgPhyleticage, ageTable, by="Taxon")

### Data for gene expression
if(!twoDataSets){
  if(!FPKM){
    orgExpression <- orgExpression[,c("Gene.ID","Stage.name",values)] 
    nDevStages<-length(levels(orgExpression$Stage.name))
    colnames(orgExpression) <- c("Ensembl.Gene.ID","Stage.name",values)
  } else if(FPKM){
    nDevStages<-length(colnames(orgExpression))-1
  }
}
#before and after aggregate : 7833 levels of Ensembl.Gene.ID, 54'831 rows
if(!FPKM & !twoDataSets){
  if(expDataSource == "Bgee"){
    orgExpression<-aggregate(Log.of.normalized.signal.intensity~Ensembl.Gene.ID+Stage.name,data=orgExpression,na.action=na.omit,FUN=mean)  
  }else if(expDataSource == "BgeeRNA"){
    orgExpression<-aggregate(RPKM~Ensembl.Gene.ID+Stage.name,data=orgExpression,na.action=na.omit,FUN=mean)  
  }
}

#For the RNA-seq : RPKM values that are <1 => considered as not expressed 
if(!twoDataSets){
  if(!FPKM){
    if(expDataSource == "BgeeRNA"){
      orgExpression$RPKM[orgExpression$RPKM<1] <- 0
      min(orgExpression$RPKM[orgExpression$RPKM>0]) #check : must be > 1.000113 #Caeno : 1.000028
    }
  } else if(FPKM){
    orgExpression[,-1] <- apply(orgExpression[,-1],c(1,2),function(x){ifelse(x<1,0,x)})
  }
}
require(reshape2)

### Change the shape of the dataset
if(!FPKM & !twoDataSets){
  orgExpression <- dcast(orgExpression, Ensembl.Gene.ID~Stage.name, value.var=values)
}

if(twoDataSets){
  orgExpression <- importAfterComBat(organism)
}

if(organism=="Mus"){
  orgExpression<-orgExpression[,c("Ensembl.Gene.ID","zygote stage", "Theiler stage 02 (mouse)", "Theiler stage 03 (mouse)",
                                  "Theiler stage 04 (mouse)", "Theiler stage 05 (mouse)", "Theiler stage 06 (mouse)",
                                  "neurula stage", "Theiler stage 13 (mouse)", "Theiler stage 15 (mouse)", 
                                  "Theiler stage 17 (mouse)", "Theiler stage 20 (mouse)", "Theiler stage 22 (mouse)",
                                  "Theiler stage 24 (mouse)","Theiler stage 26 (mouse)","post-juvenile adult stage")]  
} else if(organism=="Danio"){
  orgExpression <- orgExpression[,c("Ensembl.Gene.ID","zygote stage","Gastrula:Shield (Danio)", "Gastrula:75%-epiboly (Danio)", 
                                    "Gastrula:90%-epiboly (Danio)", "Gastrula:Bud (Danio)", "Segmentation:5-9 somites (Danio)",
                                    "Segmentation:14-19 somites (Danio)", "Pharyngula:Prim-5 (Danio)", "Pharyngula:Prim-15 (Danio)",
                                    "Hatching:Long-pec (Danio)", "Larval:Day 4 (Danio)", "Larval:Day 5 (Danio)",
                                    "Larval:Days 14-20 (Danio)","Juvenile:Days 30-44 (Danio)")] 
} else if(organism=="Drosophila" && !FPKM){
  orgExpression <- orgExpression[,c("Ensembl.Gene.ID","embryonic stage 2 (Drosophila)", "embryonic stage 4 (Drosophila)",
                                    "embryonic stage 5 (Drosophila)", "early extended germ band stage (Drosophila)",
                                    "embryonic stage 10 (Drosophila)", "late extended germ band stage (Drosophila)", 
                                    "embryonic stage 17 (Drosophila)", "first instar larval stage (Drosophila)", 
                                    "second instar larval stage (Drosophila)")]                                 
} else if(organism=="Xenopus"){
  orgExpression<-orgExpression[,c("Ensembl.Gene.ID","2 cell stage", "4 cell stage", "8 cell stage",
                                  "NF stage 5 (16-cell) (Xenopus)", "NF stage 6 (32-cell) (Xenopus)",
                                  "NF stage 8 (Xenopus)", "NF stage 9 (Xenopus)", "NF stage 10 (Xenopus)",
                                  "gastrula stage",
                                  "NF stage 15 (Xenopus)", "NF stage 19 (Xenopus)", "NF stage 28 (Xenopus)",
                                  "NF stage 33 and 34 (Xenopus)", "NF stage 40 (Xenopus)")]
  
} else if (organism=="Caeno" && !FPKM){
  orgExpression<-orgExpression[,c("Ensembl.Gene.ID",
                                  "blastula stage", "gastrula stage", "embryo stage",
                                  "nematode larval stage", "fully formed stage")]
} 
orgExpression <- orgExpression[,c(1:(length(graphParamModel(organism)$stageNum)+1))]
devStageNames <- colnames(orgExpression)[-1]
colnames(orgExpression) <- c("Ensembl.Gene.ID", as.character((graphParamModel(organism)$stageNum)))
orgExpression <- melt(orgExpression)
colnames(orgExpression) <- c("Ensembl.Gene.ID", "Stage.Num", "Expression")

### Collect files
totalAge <- merge (orgExpression, orgPhyleticage, by=c("Ensembl.Gene.ID"), all.x=T, sort=F)
totalAge <- na.omit(totalAge)

#### Plot 
ageName <- levels(as.factor(totalAge$Taxon))
ageY <- levels(as.factor(totalAge$Age))

totalAge$Stage.Num <- as.numeric(as.character(totalAge$Stage.Num))

#Draw all ages in one plot
x11()
png(paste0(folderAnalysis,"PhyleticExpression",organism, ".png"))
plot(1, type="n", axes=T, xlab=graphParamModel(organism)$sc, ylab="Relative values of mean expression", 
     xlim=range(totalAge$Stage.Num), ylim=c(-0.1,1))
for(i in ageName){
  data <- subset(totalAge[which(totalAge$Taxon==i),])
  data <- aggregate(Expression~Stage.Num, data=data, FUN=mean) 
  data$Expression <- data$Expression/max(data$Expression)
  points(data$Expression~data$Stage.Num, col=match(i, ageName))
  points(data$Expression~data$Stage.Num, type="l", col=match(i, ageName))
}

findAge <- function(i){
  x <- ageTable[which(ageTable$Taxon==i),]
  return (x[[2]])
}
ageY <- lapply(ageName, FUN=findAge)
ageNameY <- paste0(ageName, " (", ageY, ")")
legend("bottomright", c(ageNameY), lty=c(rep(1,length(ageName))),col=c(1:length(ageName)), cex=1.2)
# may need to adjust cex=... separately !

dev.off()

#Draw separate graph for each phyletic age
x11()
par(mfrow=c(6,3))
pdf(paste0(folderAnalysis,"PhyleticExpression",organism, ".pdf"))
for(i in ageName){
  data <- subset(totalAge[which(totalAge$Taxon==i),])
  data <- aggregate(Expression~Stage.Num, data=data, FUN=mean) 
  data$Expression <- data$Expression/max(data$Expression)
  plot(1, type="n", axes=T, xlab=graphParamModel(organism)$sc, ylab="Relative values of mean expression", 
       xlim=range(totalAge$Stage.Num), ylim=c(-0.1,1))
  points(data$Expression~data$Stage.Num, col=match(i, ageName))
  points(data$Expression~data$Stage.Num, type="l", col=match(i, ageName))

  ageY <- lapply(ageName, FUN=findAge)
  ageNameY <- paste0(ageName, "(", ageY, ")")
  
  legend("bottomright", c(paste0(i, " (",findAge(i),")")), lty=1, col=match(i, ageName), cex=0.8)
}
dev.off()

# ---------------------------------------------------------------------------------------
# PCA and heatmaps (1 dataset)
# ---------------------------------------------------------------------------------------
#Include comparison: without, QN normalization
#Function defined in "scriptFunctions.R"
#If you get error like "X1" not found, look at the function definition in "scriptFunctions.R"
# !!! check : this function (re)defines folder path !!! (maybe have to change it)
#For the size of legend text: may need to adjust cex=... separately !
if(!twoDataSets){
  plotHeatmapPCAoneSet(organism)
}

# ---------------------------------------------------------------------------------------
# PCA and heatmaps (2 datasets: mouse and drosophila microarray data)
# ---------------------------------------------------------------------------------------
#Include comparison: without, Combat, ComBat+QN, QN+ComBat normalization
#Function defined in "scriptFunctions.R"
#If you get error like "X1" not found, look at the function definition in "scriptFunctions.R"
# !!! check : this function (re)defines folder path !!! (maybe have to change it)
#For the size of legend text: may need to adjust cex=... separately !
if(twoDataSets){
  plotHeatmapPCAtwoSets(organism)
}
