#Clayton Nevins
#PM III 

#Libraries needed for analysis
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("limma")
biocLite("maSigPro")
biocLite("statmod")
biocLite("NOISeq")
biocLite("EDASeq")
library(NOISeq) 
library(EDASeq) 
library(limma)
library(edgeR)
library(DESeq2)

##I need to read in the htseq counts results for each of the 6 samples. Also read in the GC content and length files provided
##Somehow combine the htseq counts files. SOrt them first maybe?

ref_counts_union1 = read.delim("ref_counts_union1.txt", as.is = TRUE)
ref_counts_union2 = read.delim("ref_counts_union2.txt", as.is = TRUE)
ref_counts_union3 = read.delim("ref_counts_union3.txt", as.is = TRUE)
ref_counts_union4 = read.delim("ref_counts_union4.txt", as.is = TRUE)
ref_counts_union5 = read.delim("ref_counts_union5.txt", as.is = TRUE)
ref_counts_union6 = read.delim("ref_counts_union6.txt", as.is = TRUE)
length3 = read.delim("gene_transcript_length.txt", as.is = TRUE)
length33 <- length3[,-2]
head(length33)
colnames(length33) <- c("gene", "length")


GC3 = read.delim("gene_GC_content.txt", as.is = TRUE)
head(GC3)
colnames(GC3) <- c("gene", "GCcontent")

colnames(ref_counts_union1) <- c("g", "A_1")
colnames(ref_counts_union2) <- c("g", "A_2")
colnames(ref_counts_union3) <- c("g", "A_3")
colnames(ref_counts_union4) <- c("g", "B_1")
colnames(ref_counts_union5) <- c("g", "B_2")
colnames(ref_counts_union6) <- c("g", "B_3")

common <- intersect(ref_counts_union1$g, GC3$gene)  
head(common)

newGC3 <- GC3[GC3$gene %in% common, ]
ref_counts_union1 <- ref_counts_union1[ref_counts_union1$g %in% common, ]
ref_counts_union2 <- ref_counts_union2[ref_counts_union2$g %in% common, ]
ref_counts_union3 <- ref_counts_union3[ref_counts_union3$g %in% common, ]
ref_counts_union4 <- ref_counts_union4[ref_counts_union4$g %in% common, ]
ref_counts_union5 <- ref_counts_union5[ref_counts_union5$g %in% common, ]
ref_counts_union6 <- ref_counts_union6[ref_counts_union6$g %in% common, ]
newlength33 <- length33[length33$gene %in% common, ]
newlength33 <- aggregate(newlength33[, -c(1)], by = list(newlength33$gene),
          mean, na.rm = TRUE)
colnames(length33) <- c("gene", "length")

newr1 <- ref_counts_union1[,2]
newr2 <- ref_counts_union2[,2]
newr3 <- ref_counts_union3[,2]
newr4 <- ref_counts_union4[,2]
newr5 <- ref_counts_union5[,2]
newr6 <- ref_counts_union6[,2]

counts <- cbind(newr1, newr2, newr3, newr4, newr4, newr6)
colnames(counts) <- c("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")
head(counts)
row.names(counts) <- ref_counts_union1$g
head(counts)

mfactors = data.frame("treat" = substr(colnames(counts), start = 1, stop = 1),
                       "biorep" = substr(colnames(counts), start = 3, stop = 3))

rownames(mfactors) = colnames(counts)
head(mfactors)
mfactors
######


# Preparing data for NOISeq package
mdata = NOISeq::readData(data = counts, factors = mfactors, gc = newGC3, length = newlength33)

## Length bias
mylengthbias = dat(mdata, factor = "treat", norm = FALSE, type = "lengthbias")

par(mfrow = c(1,2))
explo.plot(mylengthbias, samples = 1)
explo.plot(mylengthbias, samples = 2)

## GC content bias
myGCbias = dat(mdata, factor = "treat", norm = FALSE, type = "GCbias")

par(mfrow = c(1,2))
explo.plot(myGCbias, samples = 1)
explo.plot(myGCbias, samples = 2)


###Normalization

mcounts <- counts
head(mcounts)
head(newGC3)
row.names(newGC3) <- newGC3$gene
newGC33 <- newGC3[with(newGC3, order(gene)), ]

colnames(newlength33) <- c("gene", "length")
row.names(newlength33) <- newlength33$gene
newlength33 <- newlength33[with(newlength33, order(gene)), ]

## Preparing data for EDASeq package
edadata = newSeqExpressionSet(counts=as.matrix(mcounts), featureData=newlength33[,2, drop = FALSE],
                              phenoData=mfactors[,"treat", drop = FALSE])

## Correcting GC content bias with EDASeq package: loess method
edadata <- withinLaneNormalization(edadata,"length", which="full")
mynormdata = edadata@assayData$normalizedCounts

#edadata = newSeqExpressionSet(counts=as.matrix(mcounts), featureData=newlength33[,2, drop = FALSE],
#                           phenoData=mfactors[,"treat", drop = FALSE])

## Correcting GC content bias with EDASeq package: loess method
#edadata <- withinLaneNormalization(edadata,"Length", which="full")
#mynormdata = edadata@assayData$normalizedCounts

## Applying TMM normalization (between-samples) with NOISeq package
mynormdata = tmm(mynormdata)
head(mynormdata)

# More uniform distributions by TMM
par(mfrow = c(1,2))
boxplot(log(as.matrix(mcounts+1)) ~ col(mcounts), main = "Before normalization")
boxplot(log(mynormdata+1) ~ col(mynormdata), main = "After normalization")

# GC bias memoved
#mydata.norm = NOISeq::readData(data = mynormdata, factors = mfactors, gc = newGC33, length = newlength33)
#myGCbias = dat(mydata.norm, factor = "treat", norm = TRUE, type = "GCbias")
#par(mfrow = c(1,2))
#explo.plot(myGCbias, samples = 1)
#explo.plot(myGCbias, samples = 1:2)
#head(mcounts)

#length bias removed
#mylengthbias = dat(mynormdata, factor = "treat", norm = TRUE, type = "lengthbias")
#par(mfrow = c(1,2))
#explo.plot(mylengthbias, samples = 1)
#explo.plot(mylengthbias, samples = 2)

#newpca <- PCA.GENES(mynormdata)
#plot(prcomp(mynormdata, scale. = TRUE))
#biplot(prcomp(mynormdata, scale. = TRUE))



new <- mynormdata[apply(mynormdata[,-1], 1, function(x) !all(x==0)),]

#newCNscaled <- CNScaled[complete.cases(CNScaled), ]

PCA <- princomp(new, cor = TRUE) 
PCA$loadings[,1]
head(PCA$scores)
PCA$sdev
Tv <- cumsum(PCA$sdev^2/sum(PCA$sdev^2))
Tv

PCA.plot <- function (PCA, main, col = rep(c(1:6),each = 3)) {
  PCAl <- PCA$loadings
  Tvar <- cumsum(PCA$sdev^2/sum(PCA$sdev^2))
  vPC1 <- round(Tvar[1] * 100)
  vPC2 <- round((Tvar[2] - Tvar[1])*100) 
  plot (PCAl[,1], PCAl[,2], 
        col = "grey", 
        main = main, 
        xlab = paste("PC1: ", vPC1, "% explained variance.", sep = ""),
        ylab = paste("PC2: ", vPC2, "% explained variance.", sep = "")
  ) 
  text (PCAl[,1], PCAl[,2], rownames(PCAl), offset = 0,
        col = col, cex = 0.75
  )
}

PCA.plot(PCA, main = "Treatment and time response PCA")

library("corpcor")

#Scaling 
CNScaled <- t(apply(new, 1, scale, scale = TRUE))
colnames(CNScaled) <- colnames(new)
newCNscaled <- CNScaled[complete.cases(CNScaled), ]
#newCNscaled <- na.omit(CNScaled)
#make.positive.definite(newCNscaled)
#cov2cor(newCNscaled)
PCA1 <- princomp(newCNscaled)
PCA.plot(PCA1, main = "Treatment response PCA")

save(mcounts, mynormdata, newGC3, newlength33, mfactors, file = "PMthreedata.RData")





#Usinging count data because this is required for deseq and edge r, optional for noiseq
mcond2compare = colnames(mcounts)[c(grep("A", colnames(mcounts)), grep("B", colnames(mcounts)))]
mcond2compare
mygroups = rep(c("A", "B"), each = 3)
mygroups

### DESeq2
mDESeq = DESeqDataSetFromMatrix(mcounts[,mcond2compare], 
                                 colData = data.frame("condi" = mygroups), design = ~condi)
myDESeq = DESeq(mDESeq)
myDESeq = results(myDESeq, contrast = c("condi","A","B"))
head(myDESeq)
summary(myDESeq) #UP-7079, DOWN-5641

myedgeR = DGEList(counts = mcounts[,mcond2compare], group = mygroups)
myedgeR = calcNormFactors(myedgeR)  
myedgeR = estimateCommonDisp(myedgeR)
myedgeR = estimateTagwiseDisp(myedgeR, trend = "movingave")
myedgeR = exactTest(myedgeR)
names(myedgeR)
head(myedgeR$table)
newtab <- myedgeR$table #####################################
head(newtab)
topTags(myedgeR)
myedgeR_decide = decideTestsDGE(myedgeR, adjust.method = "BH", p.value = 0.05, lfc = 0)
summary(myedgeR_decide) #DEGs upreg 3850, downreg 4428  


### NOISeq
mynoiseqbio = NOISeq::readData(mcounts[,mcond2compare], factors = data.frame("condi" = mygroups))
mynoiseqbio = noiseqbio(mynoiseqbio, norm = "n", k = NULL, factor = "condi", r = 30)
mynoiseqbio = degenes(mynoiseqbio, q = 0.95)
head(mynoiseqbio1) ###############
table(mynoiseqbio[,"theta"] < 0) #UP-8819, DOWN-9124


subset(newtab, select= PValue < 0.05)
x <- 0.05
newmyedgeR <- newtab[newtab[,ncol(newtab)]<x,]
#mynoiseqbio
library('VennDiagram')

newmyedgeR <- newmyedgeR[,-2]
newmyedgeR <- newmyedgeR[,2]

cleanedgeR <- newmyedgeR[,0]
cleannoiseq <- mynoiseqbio1[,0]
#colnames(cleanedgeR) <- c('DEGs')
cleandeseq <- myDESeq[,0]

cleanedgeR$DEGs <- rownames(cleanedgeR)
cleannoiseq$DEGs <- rownames(cleannoiseq)

CNSEQ <- as.matrix(cleannoiseq)
write.csv(CNSEQ, file='CNSEQ.txt')

EDGE <- as.matrix(cleanedgeR)
write.csv(EDGE, file='EDGE.txt')

cccc <- intersect(cleanedgeR$DEGs, cleannoiseq$DEGs)
head(cccc)
length(cccc) #11946 DEGs shared
cnts <- vennCounts(cccc)


c3 <- cbind(cleanedgeR$DEGs, cleannoiseq$DEGs)
colnames(c3) <- c('edgeR', 'noiseq')
head(c3)


save(cleanedgeR, cleannoiseq, mynoiseqbio, newmyedgeR, cccc, mcounts, mynormdata, newGC3, newlength33, mfactors, file = "PMthreedata2.RData")



mynoiseqbio.all = degenes(mynoiseqbio, q = 0)
expr.data <- cbind (name = rownames(mynoiseqbio.all), FC = mynoiseqbio.all[,"log2FC"])
head(expr.data)
write.table( expr.data, "expr.data.txt", row.names = F, col.names = T, quote = F, sep = "\t" )
mynoiseqbio.deg = degenes(mynoiseqbio, q = 0.99)
head(mynoiseqbio.deg)
mynoiseqbio.deg <- mynoiseqbio.deg[abs(mynoiseqbio.deg["log2FC"])> 0.3,]
dim(mynoiseqbio.deg)

write.table(rownames(mynoiseqbio.deg), "DE.data.txt", row.names = F, col.names = F, quote = F, sep = "\t" )




