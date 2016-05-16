source("http://bioconductor.org/biocLite.R")
biocLite("ReadqPCR")
biocLite("NormqPCR")
biocLite("ddCt")
biocLite("HTqPCR")

require("ReadqPCR")
require("NormqPCR")
require("ddCt")
require("HTqPCR")

vignette("ReadqPCR")
vignette("NormqPCR")
vignette("ddCt")
vignette("HTqPCR")

path <- system.file("exData", package = "NormqPCR")
taqman.example <- file.path(path, "/example.txt")
qPCRBatch.taqman <- read.taqman(taqman.example)

sampleNames(qPCRBatch.taqman)
exprs(qPCRBatch.taqman)

dat = new("qPCRBatch")
exprs(dat) <- exprs(qPCRBatch.taqman)
sampleNames(dat) <- colnames(qPCRBatch.taqman)
dat_norm = deltaCt(dat, c("Actb.Rn00667869_m1"), combineHkgs=FALSE, calc="arith")
exprs(dat_norm)

contM <- cbind(c(0,0,1,1,0,0,1,1),c(1,1,0,0,1,1,0,0))
colnames(contM) <- c("interestingPhenotype","wildTypePhenotype")
rownames(contM) <- colnames(qPCRBatch.taqman)
res = deltaDeltaCt(dat, hkgs=c("Actb.Rn00667869_m1","B2m.Rn00560865_m1","Gapdh.Rn99999916_s1"), contrastM=contM, case="interestingPhenotype",control="wildTypePhenotype",paired=F,statCalc="geom",hkgCalc="arith",maxNACase=1,maxNAControl=1)
head(res)



## Read data:
dat = new("qPCRBatch")
exprs(dat) <- read.table("/Users/robk/Downloads/Ionita_qPCR_3T3L1exp2-rawdata.txt",header=T,stringsAsFactors=F,sep="\t",row.names=1)
rownames(exprs(dat)) = gsub(" ","_",rownames(exprs(dat)))
rownames(exprs(dat)) = gsub("/","",rownames(exprs(dat)))

## assess housekeepers:
housekeepers = c("Mouse_U6_snRNA","RNU43_snoRNA","HmMsRt_U1_snRNA")
selectHKs(dat[housekeepers,], method="geNorm", Symbols=housekeepers, minNrHK=2, log=TRUE)

## normalise:
hkgs.chosen <- housekeepers
qPCRBatch.dCq <- deltaCq(qPCRBatch=dat, hkgs=hkgs.chosen) # With more than one hkg
exprs(qPCRBatch.dCq)

# colnames(exprs(qPCRBatch.dCq)) # look at the sample names so you can make a and b (below) match up
# a <- 5:6 # case
# b <- 1:2 # control
# b <- 3:4 # control
# 
# pvals <- data.frame(row.names = featureNames(qPCRBatch.dCq)) # dataframe to contain p vals from t-tests
# for(i in row.names(pvals)) {
#   case <- exprs(qPCRBatch.dCq)[i,a]
#   control <- exprs(qPCRBatch.dCq)[i,b]
#   if(sum(is.na(case)) != 0 | sum(is.na(control)) != 0) { 
#     pvals[i,1] <- NA
#   }else{
#     pvals[i,1] <- wilcox.test(case,control)$p.value
#   }
# }


## ddCT DEX
contM = diag(3)[c(1,1,2,2,3,3),]
colnames(contM) <- c("cells","EVfree","EV")
rownames(contM) <- colnames(exprs(dat))
case="EV"; control="cells"
case="EVfree"; control="cells"
res = deltaDeltaCt(dat, hkgs=c("Mouse_U6_snRNA","RNU43_snoRNA","HmMsRt_U1_snRNA"), contrastM=contM, case=case,control=control,paired=F,statCalc="geom",hkgCalc="geom",maxNACase=0,maxNAControl=0)
head(res)

# Block G:  Plot of deltaDeltaCq values
plotVals <- log2(as.numeric(as.character(res[,6])))
plotVals[as.character(res[,6]) == "-                 "] <- min(plotVals, na.rm=TRUE)
plotVals[as.character(res[,6]) == "+                 "] <- max(plotVals, na.rm=TRUE)
minusOneSD <- log2(as.numeric(as.character(res[,7])))
plusOneSD <- log2(as.numeric(as.character(res[,8])))

colVec <- rep("blue", length(plotVals))
colVec[as.character(res[,6]) == "-                 " | as.character(res[,6]) == "+                 "] <- "red"

library(gplots)
#png("qPCRdDCq.png"); 
barplot2(plotVals, col = colVec, ci.u = plusOneSD, ci.l = minusOneSD, plot.ci = TRUE, ylab=paste("delta delta Cq (",control,"/",case,")",sep=""), xlab = "miRNAs")
#dev.off()
