library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma) #perform differential expression analysis using surrogate variables either directly or
               #with the limma package

# data should be a matrix with features (genes, transcripts, voxels) in the rows and samples
# in the columns.

# 2  types  of  variables  that  are  being  considered:  
# (1)  adjustment  variables; e.g. age/sex  of  the  patients, 
# a variable like the date the arrays were processed.
# (2) variables of interest; e.g. cancer  versus  control.

# Two Model matrices must be made:  
# (1)the full Model : the full Model includes terms for both the adjustment
# variables and the variables of interest.
# (2)and the null Model : Model matrix that includes terms for all of the 
# adjustment variables but not the variables of interest.  
#=>analyze the association between the variables of interest and gene expression,
# adjusting for the adjustment variables.  The Model matrices can be created using Model.matrix().

# For  the  bladder  cancer  study,  the  variable  of  interest  is  cancer  status.   
# To  begin  we  will  assume  no adjustment variables.
pheno <- pData(bladderEset) #the variable of interest (cancer status)
edata <- exprs(bladderEset) #expression data

#we  create  the  full  Model  matrix  -  including  both  the  adjustment  variables  
# and  the  variable  of interest  (cancer  status).   In  this  case  we  only  have  the  
# variable  of  interest.   Since  cancer  status  has multiple levels, we treat it as a factor variable
Mod <- Model.matrix(~as.factor(cancer), data=pheno)

#null Model contains only the adjustment variables.  Since we are not adjusting for any other variables
# in this analysis, only an intercept is included in the Model
Mod0 <- Model.matrix(~1,data=pheno)

# number of latent factors that need to  be  estimated
n.sv <- num.sv(edata,Mod,method="leek")
svobj <- sva(edata,Mod,Mod0) #svobj <- sva(edata,Mod,Mod0, n.sv=n.sv)

pValues <- f.pvalue(edata,Mod,Mod0)
qValues <- p.adjust(pValues,method="BH")

ModSv <- cbind(Mod,svobj$sv)
Mod0Sv <- cbind(Mod0,svobj$sv)
pValuesSv <- f.pvalue(edata,ModSv,Mod0Sv)
qValuesSv <- p.adjust(pValuesSv,method="BH")

fit <- lmFit(edata,ModSv)

contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)
# The next step is to calculate the test statistics using the eBayes function:
eb <- eBayes(fitContrasts)
topTableF(eb, adjust="BH")

batch <- pheno$batch

Modcombat <- Model.matrix(~1, data=pheno)

combat_edata <- ComBat(dat=edata, batch=batch, Mod=Modcombat, par.prior=TRUE, prior.plots=FALSE)

pValuesComBat <- f.pvalue(combat_edata,Mod,Mod0)
qValuesComBat <- p.adjust(pValuesComBat,method="BH")

ModBatch <- Model.matrix(~as.factor(cancer) + as.factor(batch),data=pheno)
Mod0Batch <- Model.matrix(~as.factor(batch),data=pheno)
pValuesBatch <- f.pvalue(edata,ModBatch,Mod0Batch)
qValuesBatch <- p.adjust(pValuesBatch,method="BH")

n.sv <- num.sv(edata,Mod,vfilter=2000,method="leek")
svobj <- sva(edata,Mod,Mod0,n.sv=n.sv,vfilter=2000)
