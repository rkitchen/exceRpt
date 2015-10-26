## Read Tuschl:
mapping.stats.Tuschl = t(read.table("~/WORK/YALE_offline/exRNA/TomTuschl/exceRpt_readMappingSummary.txt", sep="\t", header=T, stringsAsFactors=F, row.names=1))
## Read Amy Buck's data
mapping.stats.Buck = t(read.table("~/WORK/YALE_offline/exRNA/AmyBuck/Results/exceRpt_readMappingSummary.txt", sep="\t", header=T, stringsAsFactors=F, row.names=1))

## Read ERCC QC data
dir="/Users/robk/WORK/YALE_offline/exRNA/QC"
files = dir(dir)[grep("stats$",dir(dir))]
files_short = gsub(".stats","",gsub(".fastq.stats","",files))
mapping.stats.QC = matrix(NA,nrow=nrow(mapping.stats.Tuschl),ncol=length(files),dimnames=list(rownames(mapping.stats.Tuschl), files_short))
for(i in 1:length(files)){
  tmp = read.table(paste(dir, files[i], sep="/"),stringsAsFactors=F,header=T,row.names=1)
  mapping.stats.QC[match(rownames(tmp), rownames(mapping.stats.QC)), i] = tmp[,1]
}


## merge
mapping.stats = t(cbind(mapping.stats.QC, mapping.stats.Tuschl, mapping.stats.Buck))
mapping.stats = as.data.frame(mapping.stats[,-grep("input_to_",colnames(mapping.stats))])
mapping.stats = na.omit(mapping.stats[,-5])

##
## Order samples based on similarity of mapping statistics
##
sampleOrder = 1
if(nrow(mapping.stats) > 1){
  h = hclust(dist(1-cor(na.omit(t(mapping.stats)))))
  sampleOrder = h$order
}


#mapping.stats = as.data.frame(mapping.stats.all)


##
## Plot heatmap of mapping percentages through the pipeline
##
toplot = melt(as.matrix(mapping.stats / mapping.stats$input)); colnames(toplot) = c("Sample","Stage","ReadFraction")
toplot$Stage = with(toplot, factor(Stage, levels = rev(levels(Stage))))
toplot$Sample = factor(as.character(toplot$Sample), levels=rownames(mapping.stats)[sampleOrder])
p = ggplot(toplot, aes(x=Sample, y=Stage, group=Sample, fill=ReadFraction, label=sprintf("%1.1f%%",ReadFraction*100))) +geom_tile() +scale_fill_gradient2(low="white",high="yellow",mid="steelblue", midpoint=0.5) +theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +ggtitle("fraction aligned reads (normalised by # input reads)")
if(nrow(mapping.stats) < 50){ p = p +geom_text(size=3) }
p



##
## Plot heatmap of mapping percentages through the pipeline
##
toplot = melt(as.matrix(mapping.stats / mapping.stats$reads_used_for_alignment)[,-c(1:6)]); colnames(toplot) = c("Sample","Stage","ReadFraction")
toplot$Stage = with(toplot, factor(Stage, levels = rev(levels(Stage))))
toplot$Sample = factor(as.character(toplot$Sample), levels=rownames(mapping.stats)[sampleOrder])
p = ggplot(toplot, aes(x=Sample, y=Stage, group=Sample, fill=ReadFraction, label=sprintf("%1.1f%%",ReadFraction*100))) +geom_tile() +scale_fill_gradient2(low="white",high="yellow",mid="steelblue", midpoint=0.5) +theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +ggtitle("fraction aligned reads (normalised by # non-contaminant reads)")
if(nrow(mapping.stats) < 50){ p = p +geom_text(size=3) }
p


