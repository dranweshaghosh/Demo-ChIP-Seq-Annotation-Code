library(GenomicRanges)
library(rtracklayer)
library(IRanges)
# Obtaining the peaks' bed data
peaks.438M <- read.table("438M.bed", skip = 0)
peaks.443M <- read.table("443M.bed", skip = 0)
peaks.438M <- RangedData(space=peaks.438M[,1], IRanges(start=peaks.438M[,2], end=peaks.438M[,3]), strand=peaks.438M[,6])
peaks.443M <- RangedData(space=peaks.443M[,1], IRanges(start=peaks.443M[,2], end=peaks.443M[,3]), strand=peaks.443M[,6])
# Generating Histogram of Width of Peaks
width <- unlist(end(ranges(peaks.438M)) - start(ranges(peaks.438M)))
hist(width, col = 3, breaks = 10, main = paste("Width Histogram of 438M"))
# Annotating Peaks and Looking at Features Where Peaks are Found
library(GenomicFeatures)
library(ChIPpeakAnno)
data("TSS.human.GRCh37")
annotatedPeak = annotatePeakInBatch(peaks.438M, AnnotationData=TSS.human.GRCh37)
temp = as.data.frame(annotatedPeak)
write.table(temp, "DS_NPC_annotatedPeaks.xls", sep="\t", row.names=FALSE)
b<- addGeneIDs(annotatedPeak,"org.Hs.eg.db",c("symbol"))
c<- as.data.frame(b)
write.table(c,file="DS_NPC_annotatedPeakList_GeneId.xls", sep="\t",
            col.names=TRUE, row.names=FALSE)
write.table(c,file="DS_NPC_annotatedPeakList_GeneId.bed", sep="\t",
            col.names=TRUE, row.names=FALSE)

annotatedPeak = annotatePeakInBatch(peaks.443M, AnnotationData=TSS.human.GRCh37)
temp = as.data.frame(annotatedPeak)
write.table(temp, "Normal_NPC_annotatedPeaks.xls", sep="\t", row.names=FALSE)
b<- addGeneIDs(annotatedPeak,"org.Hs.eg.db",c("symbol"))
c<- as.data.frame(b)
write.table(c,file="Normal_NPC_annotatedPeakList_GeneId.xls", sep="\t",
            col.names=TRUE, row.names=FALSE)
write.table(c,file="Normal_NPC_annotatedPeakList_GeneId.bed", sep="\t",
            col.names=TRUE, row.names=FALSE)

values(annotatedPeak)
table(temp$insideFeature)
pie(table(temp$insideFeature))
z<-temp[temp$insideFeature=="overlapStart",]
y = z$distancetoFeature [!is.na(z$distancetoFeature)]
hist(y, xlab = "Distance To Nearest TSS", main = "", breaks = 100, xlim = c(min(y)-100, max(y)+100))

dataDirectory <- getwd()
peaks.438M <- import.bed(file.path(dataDirectory, '438M.bed'))
peaks.443M <- import.bed(file.path(dataDirectory, '443M.bed'))
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak(peaks.438M, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
plotAnnoPie(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno, vennpie = TRUE)
plotDistToTSS(peakAnno, title = "Distribution of Pax6 binding loci\nrelative to TSS in DS NPC")
peakAnno <- annotatePeak(peaks.443M, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db")
plotAnnoPie(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno, vennpie = TRUE)
plotDistToTSS(peakAnno, title = "Distribution of Pax6 binding loci\nrelative to TSS in normal NPCs")

library(chipseq)
ovlp = findOverlaps(peaks.438M, peaks.443M)
head(ovlp)
ov = min(length(unique(ovlp@from)), length(unique(ovlp@to)))
library(VennDiagram)
draw.pairwise.venn(
  area1 = length(peaks.438M),
  area2 = length(peaks.443M),
  cross.area = ov,
  category = c(" Down Syndrome NPC ", " Normal NPC "),
  fill = c("steelblue", "blue3"),
  cat.cex = 0.7)

# Enriched Regions = Common peaks identified on both samples
enriched.regions = Reduce(subsetByOverlaps, list(peaks.438M, peaks.443M))
library(biomaRt)
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
               dataset = "hsapiens_gene_ensembl")
listAttributes(mart)[1:3,]
filterlist <- c(1:22, 'X', 'Y')
ds = useDataset('hsapiens_gene_ensembl', mart = mart)
egs = getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position', 'strand'), 
            filters = 'chromosome_name',
            values = filterlist,
            mart = ds)

head(egs)
tail(egs)
egs$TSS = ifelse(egs$strand ==1, egs$start_position, egs$end_position)
head(egs)

# Considering +- 2000 bp around TSS

promoter_regions = GRanges(seqnames = Rle(paste0('chr', egs$chromosome_name)),
                           ranges = IRanges(start = egs$TSS - 2000, end = egs$TSS + 2000),
                           strand = Rle(rep("*", nrow(egs))),
                           gene = egs$external_gene_name)

head(promoter_regions)
common_peaks <- egs[unique(ovlp@from),]
common_peaks[1:3,]
write.table(common_peaks, file = "common_peaks.bed", sep="\t",
            col.names=TRUE, row.names=FALSE)


ovlp2 = findOverlaps(enriched.regions, promoter_regions)
cat(sprintf("%d of %d promoters are overlapped by an enriched region.", 
    length(unique(ovlp2@from)), length(promoter_regions)))
#76 of 60467 promoters are overlapped by an enriched region.
ovlp2 = findOverlaps(peaks.438M, promoter_regions)
cat(sprintf("%d of %d promoters are overlapped by 438M.", 
            length(unique(ovlp2@from)), length(promoter_regions)))
#841 of 60467 promoters are overlapped by 438M.
ovlp2 = findOverlaps(peaks.443M, promoter_regions)
cat(sprintf("%d of %d promoters are overlapped by 443M.", 
            length(unique(ovlp2@from)), length(promoter_regions)))
#269 of 60467 promoters are overlapped by 443M.

ovlp2 = findOverlaps(promoter_regions,enriched.regions)
cat(sprintf("%d of %d enriched regions overlap a promoter.", length(unique(ovlp2@from)), length(enriched.regions)))
#86 of 1134 enriched regions overlap a promoter.
ovlp2 = findOverlaps(promoter_regions,peaks.438M)
cat(sprintf("%d of %d peaks in 438M overlap a promoter.", length(unique(ovlp2@from)), length(peaks.438M)))
#950 of 11620 peaks in 438M overlap a promoter.
ovlp2 = findOverlaps(promoter_regions,peaks.443M)
cat(sprintf("%d of %d peaks in 443M overlap a promoter.", length(unique(ovlp2@from)), length(peaks.443M)))
#309 of 3368 peaks in 443M overlap a promoter.

#Which promoters are overlapped by peaks?
pos.TSS = egs[unique(findOverlaps(promoter_regions, peaks.438M)@from),]
pos.TSS[1:3,]
write.table(pos.TSS, file = "DS_NPC_promoters.bed", sep="\t",
            col.names=TRUE, row.names=FALSE)
pos.TSS = egs[unique(findOverlaps(promoter_regions, peaks.443M)@from),]
pos.TSS[1:3,]
write.table(pos.TSS, file = "Normal_NPC_promoters.bed", sep="\t",
            col.names=TRUE, row.names=FALSE)


# Generating heat maps
library(clusterProfiler)
covplot(peaks.438M)
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(peaks.438M, windows = promoter)
heatmaply(tagMatrix, xlim=c(-3000,3000), scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "white", high = "red", midpoint = 0.5, limits = c(0, 1)))
tagHeatmap(tagMatrix, xlim=c(-3000,3000), color = "red")
plotAvgProf(tagMatrix, xlim = c(-3000,3000), conf = 0.95, resample = 500, xlab = "Genomic Region (5' -> 3')", ylab = "Read Count Frequency")
heatmap(tagMatrix, Colv = NA, scale = "column")

# Annotating peaks and performing Gene Ontology Analysis
files = c("/Users/anweshaghosh/Documents/Sample/Peak Annotation/438M.bed", "/Users/anweshaghosh/Documents/Sample/Peak Annotation/443M.bed")
print(files)
names(files) <- c("DS NPC", "Normal NPC")
tagMatrixList = lapply(files, getTagMatrix, windows = promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000,3000))
plotAvgProf(tagMatrixList, xlim=c(-3000,3000), conf = 0.95, resample=500, facet = "row")
tagHeatmap(tagMatrixList, xlim=c(-3000,3000), color = NULL)
peakAnnoList <- lapply(files, annotatePeak, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster = genes, 
                           fun = "enrichKEGG",
                           pvalueCutoff = 0.3,
                           organism = "hsa")
plot(compKEGG, showCategory = 25, "KEGG Pathway Enrichment Analysis")
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)
compGO<- compareCluster(geneCluster = genes, 
                           fun = "enrichGO",
                           pvalueCutoff = 0.2,
                           OrgDb='org.Hs.eg.db')
plot(compGO, showCategory = 25, "Gene Ontology Analysis")
compGO<- compareCluster(geneCluster = genes, 
                        fun = "groupGO",
                        OrgDb='org.Hs.eg.db')
plot(compGO, showCategory = 25, "Gene Ontology Analysis")

compPATH<- compareCluster(geneCluster = genes, 
                        fun = "enrichDO")
plot(compPATH, showCategory = 25, "Pathway Analysis")

# Comparing peaks
write.table(pos.TSS, file = "promoters_commonoverlap.bed", quote = F, sep = "\t", row.names = F, col.names = F)
pos.TSS = egs[unique(findOverlaps(promoter_regions, peaks.438M)@from),]
pos.TSS[1:3,]
write.table(pos.TSS, file = "promoters_overlap_438M.bed", quote = F, sep = "\t", row.names = F, col.names = F)
pos.TSS = egs[unique(findOverlaps(promoter_regions, peaks.443M)@from),]
pos.TSS[1:3,]
write.table(pos.TSS, file = "promoters_overlap_443M.bed", quote = F, sep = "\t", row.names = F, col.names = F)

df_from_gr <- data.frame(seqnames=seqnames(promoter_regions),
                         starts=start(promoter_regions)-1,
                         ends=end(promoter_regions),
                         names=c(rep(".", length(promoter_regions))),
                         scores=c(rep(".", length(promoter_regions))),
                         strands=strand(promoter_regions))

write.table(df_from_gr, file = "promoters.bed", quote = F, sep = "\t", row.names = F, col.names = F)

# Possible Bedtools alternative
#bedtools intersect -a 443M.bed -b promoters.bed > 443_promoters.bed -u
#bedtools intersect -a 438M.bed -b promoters.bed > 438_promoters.bed -u

promoters.438M <- read.table("438_promoters.bed", skip = 0)
promoters.438M <- promoters.438M[1:6]
colnames(promoters.438M) <- c('chr','start','end','id','score','strand')
granges_promoters_438M <- makeGRangesFromDataFrame(promoters.438M, keep.extra.columns = TRUE)

library(org.Hs.eg.db)
annotatedPeak = annotatePeakInBatch(granges_promoters_438M, AnnotationData = TSS.human.GRCh37)
temp = as.data.frame(annotatedPeak)
write.table(temp, "438M_promoters_annotatedPeaks.xls", sep="\t", row.names=FALSE)
b<- addGeneIDs(annotatedPeak,"org.Hs.eg.db",c("symbol"))
c<- as.data.frame(b)
write.table(c,file="438M_promoters_annotatedPeakList_GeneId.xls", sep="\t",
            col.names=TRUE, row.names=FALSE)
write.table(c,file="438M_promoters_annotatedPeakList_GeneId.bed", sep="\t",
            col.names=TRUE, row.names=FALSE)

promoters.443M <- read.table("443_promoters.bed", skip = 0)
promoters.443M <- promoters.443M[1:6]
colnames(promoters.443M) <- c('chr','start','end','id','score','strand')
granges_promoters_443M <- makeGRangesFromDataFrame(promoters.443M, keep.extra.columns = TRUE)

annotatedPeak = annotatePeakInBatch(granges_promoters_443M, AnnotationData = TSS.human.GRCh37)
temp = as.data.frame(annotatedPeak)
write.table(temp, "443M_promoters_annotatedPeaks.xls", sep="\t", row.names=FALSE)
b<- addGeneIDs(annotatedPeak,"org.Hs.eg.db",c("symbol"))
c<- as.data.frame(b)
write.table(c,file="443M_promoters_annotatedPeakList_GeneId.xls", sep="\t",
            col.names=TRUE, row.names=FALSE)
write.table(c,file="443M_promoters_annotatedPeakList_GeneId.bed", sep="\t",
            col.names=TRUE, row.names=FALSE)



