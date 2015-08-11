##
## Tissue classification
##
#load(classifier.path)
#
#require(randomForest)
#
#
#tmp = exprs.miRNA.rpm
#rownames(tmp) = gsub("-","_",rownames(tmp))
#unobservedIDs = rownames(clf$importance)[!rownames(clf$importance) %in% rownames(tmp)]
#tmp = rbind(tmp, matrix(0.0,nrow=length(unobservedIDs),ncol=ncol(tmp),dimnames=list(unobservedIDs,colnames(tmp))))
#tmp = log10(tmp+0.0001)
#
#predict(clf, t(tmp))
#
#
#hist(clf$importance, breaks=seq(0,22,by=0.5))
#clf$importance[clf$importance > 1,1, drop=F]
#
#require(gplots)
#par(oma=c(2,0,0,5))
#heatmap.2(tmp[rownames(tmp) %in% rownames(clf$importance[clf$importance > 0.5,1, drop=F]), ], trace="none")
#
#
#
#tmp2 = trainingData
#colnames(tmp2) = miRNA_annotations.all$tissue
#tmp = exprs.miRNA.rpm
#keepIDs = intersect(rownames(tmp), rownames(tmp2))
#tmp = tmp[match(keepIDs, rownames(tmp)), ]
#tmp2 = tmp2[match(keepIDs, rownames(tmp2)), ]
#
#res = apply(tmp, 2, cor, tmp2); rownames(res) = colnames(tmp2)
#breaks=seq(-0.1,0.7, by=0.01)
#heatmap.2(res, trace="none",breaks=breaks,col=colorRampPalette(c("white","steelblue"))(length(breaks)-1), symbreaks=F)
#

require(DeconRNASeq)
load("~/Box Sync/Work/miRNA/DataSets/miRNA_normalised.RData")
miRNA_annotations.all = read.csv("/Users/robk/Box Sync/Work/miRNA/DataSets/run_annotations_full.csv", stringsAsFactors=F)
miRNA_annotations.all = miRNA_annotations.all[miRNA_annotations.all$ID %in% colnames(exprs.miRNA.norm), ]
exprs.miRNA.norm = exprs.miRNA.norm[, match(miRNA_annotations.all$ID, colnames(exprs.miRNA.norm))]
exprs.miRNA.norm = exprs.miRNA.norm[, miRNA_annotations.all$wt == 1]
miRNA_annotations.all = miRNA_annotations.all[miRNA_annotations.all$wt == 1, ]
tissues = unique(miRNA_annotations.all$tissue)
tissueAverages = matrix(0,nrow=nrow(exprs.miRNA.norm),ncol=length(tissues), dimnames=list(rownames(exprs.miRNA.norm), tissues))
for(i in 1:length(tissues)){
  tissueAverages[,i] = rowMeans(exprs.miRNA.norm[, colnames(exprs.miRNA.norm) %in% miRNA_annotations.all[miRNA_annotations.all$tissue %in% tissues[i], ]$ID])
}
hist(log10(apply(tissueAverages, 1, max)), breaks=100)
tissueAverages = tissueAverages[apply(tissueAverages, 1, max) >= 0.005, ]



load("/Users/robk/WORK/YALE_offline/exRNA/TomTuschl/exceRpt_smallRNAQuants_ReadsPerMillion.RData")
load("/Users/robk/Downloads/exceRpt AD CSF first40_exceRpt_smallRNAQuants_ReadsPerMillion.RData")    # Kendall AD CSF
load("/Users/robk/Downloads/exceRpt_illuminaAdapter_2015-3-2-14%3A41%3A45_smallRNAQuants_RPM.RData") # Kendall Urine/Saliva/CSF

tissueAverages.toUse = tissueAverages[rownames(tissueAverages) %in% rownames(exprs.miRNA.rpm), ]
sampleExprs = exprs.miRNA.rpm[match(rownames(tissueAverages.toUse), rownames(exprs.miRNA.rpm)), ]

res = DeconRNASeq(as.data.frame(sampleExprs), as.data.frame(tissueAverages.toUse), known.prop=F, fig=F, checksig=TRUE)
fractions = res$out.all
rownames(fractions) = colnames(sampleExprs)
par(oma=c(15,0,0,10))
breaks=seq(0,max(fractions),length.out=50)
heatmap.2(t(fractions), trace="none",breaks=breaks,col=colorRampPalette(c("white","red","yellow"))(length(breaks)-1), key.xlab="estimated fraction contributed")


total = 20477310
sigs = rbind(c(125260,162195),c(104415,301971))/20477310
solve(rbind(c(15156,12995),c(4590,0)))
res = DeconRNASeq(as.data.frame(rbind(c(15156,12995),c(4590,0))), as.data.frame(rbind(c(1,1),c(1,1))), known.prop=F, fig=F, checksig=TRUE)
require(NMF)
nmf(rbind(c(15156,12995),c(4590,0)), 1)


# 
# 
# require(plyr)
# mirna.old = read.table("/Users/robk/WORK/YALE_offline/exRNA/TESTING/HUMAN_SRR822469_TEST_hg19/mature_sense.grouped",sep="\t",header=T,stringsAsFactors=F)[,c(1,4)]
# mirna.old = ddply(mirna.old,"name",function(mat){ sum(mat[,2]) })
# mirna.new = read.table("/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/readCounts_miRNAmature_sense.txt",sep="\t",header=F,stringsAsFactors=F)
# mirna.new[,1] = sapply(mirna.new[,1], function(str){ strsplit(str,":")[[1]][1] })
# 
# x = log10(mirna.old[match(mirna.new[,1], mirna.old[,1]), 2])
# y = log10(mirna.new[,2])
# plot(x, y, main="mature miRNA read counts", xlab="sRNABench log10(count)", ylab="NEW log10(count)")
# abline(0,1,col=rgb(1,0,0,0.5))
# cor(x, y, use="complete.obs")
# 
# 
# trna.old = read.table("/Users/robk/WORK/YALE_offline/exRNA/TESTING/HUMAN_SRR822469_TEST_hg19/hg19_tRNAs_sense.grouped",sep="\t",header=T,stringsAsFactors=F)[,c(1,4)]
# trna.old = ddply(trna.old,"name",function(mat){ sum(mat[,2]) })
# trna.new = read.table("/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/readCounts_tRNA_sense.txt",sep="\t",header=F,stringsAsFactors=F)
# #trna.new[,1] = sapply(mirna.new[,1], function(str){ strsplit(str,":")[[1]][1] })
# x = log10(trna.old[match(trna.new[,1], trna.old[,1]), 2])
# y = log10(trna.new[,2])
# plot(x, y, main="tRNA read counts", xlab="sRNABench log10(count)", ylab="NEW log10(count)")
# abline(0,1,col=rgb(1,0,0,0.5))
# cor(x, y, use="complete.obs")


