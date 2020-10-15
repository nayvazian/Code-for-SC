# Code-for-SC
The code for R &amp; Linux &amp; Python

## R

### DESeq2
I used the DESeq2 to do the DEG analysis.

#### Import salmon results
```

dir <- "Documents/lab/DEseq2/quants"
list.files(dir)
samples <- paste0("Mouse2_", c(seq(1,4)),"_quant")
files <- file.path(dir, samples, "quant.sf")
names(files) <-paste0(c('Mouse 2-SVZL-P3','Mouse 2-SVZR-P2','R172H SVZ L4P2','R172H SVZ L5P2','R172H SVZ L6P2','181004#1 SVZ P2','181004#4 SVZ P2','181004#5 SVZ P2','181004#2 SVZ P2','181004#3 SVZ P2')) # sample name
all(file.exists(files))
```
#### Reference gene-level annotation package
```
library(EnsDb.Mmusculus.v79)
txdb <- EnsDb.Mmusculus.v79
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
```

#### Associated salmon output with annotation package

```
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
names(txi)
```

#### Using DESeq2 to analyze differential gene expression
```
library(DESeq2)
sampleTable <- data.frame(condition = factor(c(rep("SVZ",2), rep("Control", 8)))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
```

#### Differential expression analysis
```
dds <- DESeq(dds)
res <- results(dds, name="condition_treated_vs_untreated")
res <- results(dds, contrast=c("condition","treated","untreated"))
res
```

#### Log fold change shrinkage(reduce noises)
```
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC
```

#### p-values and adjusted p-values
```
resOrdered <- res[order(res$pvalue),] 
summary(res)
sum(res$padj < 0.1, na.rm=TRUE) 
res01 <- results(dds, alpha=0.1) 
summary(res01)
```

#### MA-plot
Plot counts: examine the counts of reads for a single gene across the group 
```
plotMA(res, ylim=c(-2,2))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
```
![example output](count.pdf)
#### Export results to CVS files
```
write.csv(as.data.frame(resOrdered), 
          file="desktop/SVZ_control.csv")
```
          
#### Extracting transformed values
```
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
```

#### Effects of transformations on the variance
standard deviation of transformed data, across samples, aginst the mean, using the shifted logarithm transformation)
```
ntd <- normTransform(dds)

library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))
```

#### Heatmap: various transformation of data
Based on the DEG analysis, we can create a heatmap to visualize the data or chose the gene you interested in. Then Convert Ensembl ID to gene name. 
```
library("pheatmap")
df <- as.data.frame(colData(dds)[c("condition")])
rownames(df) <- colnames(dds)
pheatmap(assay(ntd)[c('ENSMUSG00000000093',
                      'ENSMUSG00000000103',
                      'ENSMUSG00000000305'),], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df) # the gene you are interested in. 

mat <- assay(ntd)[c('ENSMUSG00000000093',
                      'ENSMUSG00000000103',
                      'ENSMUSG00000000305'),]
                      
library(biomaRt)
mart <- useMart("ensembl","mmusculus_gene_ensembl", host = 'uswest.ensembl.org', ensemblRedirect = FALSE)
gns <- getBM(c("external_gene_name","ensembl_gene_id"), "ensembl_gene_id", row.names(mat), mart)
row.names(mat)[match(gns[,2], row.names(mat))] <- gns[,1] #notice the order of geneiD and gene name
pheatmap(mat, show_rownames=TRUE, annotation_col=df,display_numbers =TRUE)        
```
![example output](Heatmap.png)
#### Heatmap of the sample to sample distances
```
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```
![example output](sampledis.png)

### Seurat
#### Input the data(expression matrix)
```
a <- read.csv('/Users/chensitong/Desktop/sc/dpi25/dpi25.csv',header = T, row.names = 1)
cbmc.rna <- as.sparse(a) 
pbmc <- CreateSeuratObject(counts = cbmc.rna, project = 'dpi25', min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```

#### expression level data without mt number 
```
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
#### Normalizing the data
```
pbmc <- NormalizeData(pbmc)
```
#### Identification of highly variable features (feature selection)
```
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
```
#### Identify the 10 most highly variable genes
```
top10 <- head(VariableFeatures(pbmc), 10)
```

#### plot variable features with and without labels
```
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

#### PCA
```
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
```

#### determine the significant PC
```
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
```
#### clustering 
```
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.25)## resolution determine the number of cluster
```
#### UMAP
```
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
```

####  find markers for every cluster compared to all remaining cells, report only the positive ones
```
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c('Cd3','Cd19'))##marker 
FeaturePlot(pbmc, features = c('Ppp1r14b','Acsl','Pdgfra','Olig1','Olig2','Sox2','Ccnd2','Sox11','Chd7','Tcf4','Set'))
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

new.cluster.ids <- c('RG','Astrocyte','iGC','OPC','pri-OPC','Cycling NB','GABAergic','Glutamatergic','Ependymal')
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
```
### Convert Ensembl ID to gene name. 
input the file.
```
ID <- read.csv('Desktop/late_control.csv')
```
#### Convert mouse ensebl_gene_id to mouse gene name 
```
library('biomaRt')
mart <- useMart("ensembl","mmusculus_gene_ensembl", host = 'uswest.ensembl.org', ensemblRedirect = FALSE)
genes <- ID$X
mgi <- getBM( attributes= c("ensembl_gene_id",'mgi_symbol'),values=genes,mart= mart,filters= "ensembl_gene_id")
```
#### Convert mouse gene name to human gene name
```
library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes = mgi$mgi_symbol
genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes ,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
```

### GSEA prerank input file generation. 
After DEG analysis, I will use GSEA prerank to the enrichment analysis. Here is the code for the GSEA input file.
```
f$fcsign <- sign(f$log2FoldChange)
f$logP <- -log10(f$pvalue)
f$metric = f$logP/f$fcsign
y <- f[,c('hgnc_symbol','metric')]
y <- na.omit(y)
write.table(y,file="desktop/SVZ_control.rnk",quote=F,sep="\t",row.names=F)
```
