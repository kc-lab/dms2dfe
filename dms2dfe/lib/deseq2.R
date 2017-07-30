#!/usr/bin/Rscript

library("DESeq2")

args = commandArgs(trailingOnly=TRUE)
counts_fh=args[1]
annots_fh=args[2]
tranformORnormalize=args[3]

# counts_fh="/data_lbl/aas_all/all.csv"
# annots_fh="data_fit/aas_annot/KKA2_S__Kan14_L__avg_WRT_KKA2_Bkg__avg"

colData <- read.csv(annots_fh, row.names=1, check.names = FALSE)
rownames(colData) <- gsub(" ", ".", rownames(colData))
countData <- read.csv(counts_fh,sep=",",row.names="mutids", check.names = FALSE)
rownames(colData)
head(countData)
countData <- as.matrix(countData[,rownames(colData)])

#denanrows
# countData <- countData[ -which( is.na(countData), arr.ind=TRUE )[,"row"], ]
colData <- colData[,c("condition","type")]
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
storage.mode(countData) = "integer"

#MAKE DDS
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition,
                             )

# dds$condition <- relevel(dds$condition, ref="unsel")
dds <- DESeq(dds)

if(tranformORnormalize == 2){
    dds$condition <- factor(dds$condition, levels=c("ref","sel"))
    res <- results(dds)
    write.csv(as.data.frame(res),file=paste(annots_fh,".deseq2_res.csv",sep=''))
    disp=mcols(dds)[c("baseMean", "baseVar","dispersion",
                      "dispGeneEst","dispOutlier","dispMAP","deviance","dispFit","dispIter","dispGeneEst")]
    rownames(disp)=rownames(dds)
    write.csv(disp,file=paste(annots_fh,".deseq2_disp.csv",sep=''))
    }
if(tranformORnormalize == 1){
    rld <- rlog(dds, blind=TRUE)
    write.csv(assay(rld),file=paste(annots_fh,".deseq2_rld.csv",sep=''))

    vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
    write.csv(assay(vsd),file=paste(annots_fh,".deseq2_vsd.csv",sep=''))
    }

# write.csv(counts(dds, normalized = TRUE), file=paste(annots_fh,".deseq2_counts_norm.csv",sep=''))
# write.csv(counts(dds, normalized = FALSE), file=paste(annots_fh,".deseq2_counts_raw.csv",sep=''))
