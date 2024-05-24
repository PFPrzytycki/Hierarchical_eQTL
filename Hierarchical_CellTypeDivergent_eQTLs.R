# load CellWalkR library
library("CellWalkR")

# load scATAC from Ziffra et al
x.sp <- readRDS("/data/AllPrimary_Sub_PCgenes.rds")

# load labeling data from Polioudakis et al
Polioudakis_markers <- fread("/data/Polioudakis_cluster_enriched_genes.csv")

# build network for CellWalk
Polioudakis_markers <- Polioudakis_markers[Polioudakis_markers$Gene %in% colnames(x.sp@gmat),]
labelGenes <- data.frame(gene = Polioudakis_markers$Gene,
             cluster = as.character(Polioudakis_markers$Cluster),
             logFC = Polioudakis_markers$Log2_fold_change, stringsAsFactors = FALSE)
ATACGenePeak <- mapSnapATACToGenes(labelGenes, x.sp, whichMat = "gmat")
labelEdges <- computeLabelEdges(labelGenes, x.sp, ATACGenePeak, whichMat = "gmat")
labelEdgesList <- list(labelEdges)

# precomupted jaccard distance between all cells
cellEdges <- as.matrix(fread("/data/Jaccard_dist_mat.txt"))

# tune edge weights
edgeWeights <- tuneEdgeWeights(cellEdges, 
                              labelEdgesList, 
                              labelEdgeOpts = 10^seq(-5,5,1), 
                              parallel=TRUE, numCores=3, steps=3)

# generate cell walk with optimal s paramter
cellWalk <- walkCells(cellEdges,labelEdgesList, labelEdgeWeights = 100)

# save CellWalk
saveRDS(cellWalk, "/data/CellWalk.rds")

# load CellWalk
cellWalk <- readRDS("/data/CellWalk.rds")

# Run on fetal brain eQTLs 
# utility function to convert eQTL data to GRanges
VarToGRanges <- function(x){
    rangeMat <- t(data.frame(strsplit(x[[2]], "_")))
    GRanges(rangeMat[,1], IRanges(as.numeric(rangeMat[,2]),as.numeric(rangeMat[,2])))
}
hg19_hg38_chain <- import.chain("/data/hg19ToHg38.over.chain")
# load fine-mappped eQTLs
file <- "/data/mixed_ciseqtl_90hcp_perm_purity_filtered.txt"
ciseQTLs <- fread(file)
ciseQTLs_ranges <- VarToGRanges(ciseQTLs)

# run CellWalkR labelBulk on eQTLs
labelScores <- labelBulk(cellWalk, liftOver(ciseQTLs_ranges, chain = hg19_hg38_chain), 
                            x.sp@pmat, 
                            x.sp@peak, 
                            allScores = TRUE)
write.table(labelScores, paste0(file,"_labelScores.txt"), row.names = FALSE) #raw label scores

# total number of scoreable eQTLs
length(which(!is.na(rowSums(labelScores))))
  # 11765
labelScores_filter = labelScores[!is.na(rowSums(labelScores)),]

# plot examples:
maxLabel <- sapply(1:nrow(labelScores_filter), function(i) colnames(labelScores_filter)[order(labelScores_filter[i,], decreasing = TRUE)])
p <- ggplot() + geom_point(aes(labelScores_filter$Mic, labelScores_filter$End)) + 
    geom_point(aes(labelScores_filter$Mic[labelScores_filter$oRG>1 & maxLabel=="oRG"], labelScores_filter$End[labelScores_filter$oRG>1 & maxLabel=="oRG"]), color=rainbow(31)[27]) + 
        geom_point(aes(labelScores_filter$Mic[labelScores_filter$vRG>1 & maxLabel=="vRG"], labelScores_filter$End[labelScores_filter$vRG>1 & maxLabel=="vRG"]), color=rainbow(31)[28]) +  
    geom_point(aes(labelScores_filter$Mic[labelScores_filter$Mic>2 & maxLabel=="Mic"], labelScores_filter$End[labelScores_filter$Mic>2 & maxLabel=="Mic"]), color=rainbow(31)[4]) + 
    geom_point(aes(labelScores_filter$Mic[labelScores_filter$End>2 & maxLabel=="End"], labelScores_filter$End[labelScores_filter$End>2 & maxLabel=="End"]), color=rainbow(31)[6]) + 
    xlab("Mic score") + ylab("End score") + theme_classic()
ggsave(paste0(file, "_Mic_End.png"), p, width=4, height=4)
p <- ggplot() + geom_point(aes(labelScores_filter$vRG, labelScores_filter$oRG)) + 
    geom_point(aes(labelScores_filter$vRG[labelScores_filter$oRG>1 & maxLabel=="oRG"], labelScores_filter$oRG[labelScores_filter$oRG>1 & maxLabel=="oRG"]), color=rainbow(31)[27]) + 
        geom_point(aes(labelScores_filter$vRG[labelScores_filter$vRG>1 & maxLabel=="vRG"], labelScores_filter$oRG[labelScores_filter$vRG>1 & maxLabel=="vRG"]), color=rainbow(31)[28]) + 
    geom_point(aes(labelScores_filter$vRG[labelScores_filter$Mic>2 & maxLabel=="Mic"], labelScores_filter$oRG[labelScores_filter$Mic>2 & maxLabel=="Mic"]), color=rainbow(31)[4]) + 
    geom_point(aes(labelScores_filter$vRG[labelScores_filter$End>2 & maxLabel=="End"], labelScores_filter$oRG[labelScores_filter$End>2 & maxLabel=="End"]), color=rainbow(31)[6]) + 
    xlab("vRG score") + ylab("oRG score") + theme_classic()
ggsave(paste0(file, "_vRG_oRG.png"), p, width=4, height=4)

# plot as upset plots of high confidence and low confidence cell-type eQTLs:
num2LabelScores <- apply(labelScores_filter, 1, function(x) length(which(x>2))) #z-scores >2 are high conf
num1LabelScores <- apply(labelScores_filter, 1, function(x) length(which(x>1 & x<2))) #z-scores between 1 and 2 are low
highConfidence <- labelScores_filter[num2LabelScores>0,]>2 #eQTls with at least 1 high confidence label
lowConfidence <- labelScores_filter[num1LabelScores>0 & num2LabelScores==0,]>1 #eQTls with at least 1 low confidence label
cellTypes <- colnames(labelScores)[c(1,8,9,11,12,13,5,10,14,16,3,2,15,4,7,6)] #consistent order of cell type names
pdf(file=paste0(file, "_highConfidence_upset.pdf"), 4, 6)
UpSetR::upset(data.table(highConfidence+0), sets=cellTypes, keep.order=TRUE, order.by = "freq", mainbar.y.max=3000, sets.bar.color=rev(rainbow(31)[c(2,4,6,31,30,25,28,27,22,21,19,18,11,13,15,16)]), sets.x.label = "Num eQTLs",  mb.ratio = c(0.6, 0.4))
dev.off()
pdf(file=paste0(file, "_lowConfidence_upset.pdf"), 8, 6)
UpSetR::upset(data.table(lowConfidence+0), sets=cellTypes, keep.order=TRUE, order.by = "freq", mainbar.y.max=3000, sets.bar.color=rev(rainbow(31)[c(2,4,6,31,30,25,28,27,22,21,19,18,11,13,15,16)]), main.bar.color = "gray", sets.x.label = "Num eQTLs",  mb.ratio = c(0.6, 0.4))
dev.off()

# Add Hierarchy to CellWalkR (Hierarchy taken from Polioudakis et al)
# manually build hierarchy
cellTypes <- colnames(labelEdges)
allCellTypes <- c(cellTypes, 
    paste(cellTypes[c(8,1)], collapse=" "), #1
    paste(cellTypes[c(8,1,9)], collapse=" "), #2
    paste(cellTypes[c(8,1,9,11)], collapse=" "), #3
    paste(cellTypes[c(12,13)], collapse=" "), #4
    paste(cellTypes[c(8,1,9,11,12,13)], collapse=" "), #5
    paste(cellTypes[c(5,10)], collapse=" "), #6
    paste(cellTypes[c(8,1,9,11,12,13,5,10)], collapse=" "), #7
    paste(cellTypes[c(14,16)], collapse=" "), #8
    paste(cellTypes[c(14,16,3)], collapse=" "), #9
    paste(cellTypes[c(2,15)], collapse=" "), #10
    paste(cellTypes[c(14,16,3,2,15)], collapse=" "), #11
    paste(cellTypes[c(8,1,9,11,12,13,5,10,14,16,3,2,15)], collapse=" "), #12
    paste(cellTypes[c(8,1,9,11,12,13,5,10,14,16,3,2,15,4)], collapse=" "), 
    paste(cellTypes[c(8,1,9,11,12,13,5,10,14,16,3,2,15,4,7)], collapse=" "), 
    paste(cellTypes[c(8,1,9,11,12,13,5,10,14,16,3,2,15,4,7,6)], collapse=" "))
cellTypesH <- cbind(
    c(-8, -9, -11, -12, 3,  -5, 5, -14, -3,  -2,  9,  7, -4, -7, -6),
    c(-1,  1,   2, -13, 4, -10, 6, -16,  8, -15, 10, 11, 12, 13, 14))
treeMatrix <- function(treeMerge) {
    A = Matrix(0, 2*nrow(treeMerge)+1, 2*nrow(treeMerge)+1)
    for(i in 1:nrow(treeMerge)){ #should build names here too
        children = treeMerge[i,]
        children = sapply(children, function(x) ifelse(x<0, abs(x), x+nrow(treeMerge)+1)) 
        A[children, i+nrow(treeMerge)+1] = 1
        A[i+nrow(treeMerge)+1, children] = 1
    }
    A
}
cellTypesH <- treeMatrix(cellTypesH)
colnames(cellTypesH) = rownames(cellTypesH) = allCellTypes

# add tree edges to label edges to make combined graph
l <- length(cellTypes)
l_all <- length(allCellTypes)
expandLabelEdges <- cbind(labelEdges, matrix(0,dim(labelEdges)[1],l_all-l))
weight <- 100
combinedGraph <- rbind(cbind(cellTypesH,t(weight*expandLabelEdges)),
        cbind(weight*expandLabelEdges,cellEdges))

# compute CellWalk on expanded graph
infMat <- randomWalk(combinedGraph)
normMat <- normalizeInfluence(infMat[-(1:l_all),(1:l_all)])
colnames(normMat) <- allCellTypes
cellLabels <- apply(normMat, 1, function(x) colnames(normMat)[order(x, decreasing = TRUE)][1])
cellWalkH <- list(infMat=infMat, normMat=normMat, cellLabels=cellLabels)
class(cellWalkH) <- "cellWalk"

# save CellWalk with hierarchy
saveRDS(cellWalkH, "/data/CellWalk_typeH.rds")

# load CellWalk with hierarchy
cellWalkH <- readRDS("/data/CellWalk_typeH.rds")

# need to run label bulk on hierarchy cell walk
labelScores_typeH <- labelBulk(cellWalkH, liftOver(ciseQTLs_ranges, chain = hg19_hg38_chain), 
                            x.sp@pmat, 
                            x.sp@peak, 
                            allScores = TRUE)
write.table(labelScores_typeH, paste0(file,"_labelScores_typeH.txt"), row.names=FALSE) #raw label scores
labelScores_typeH <- apply(labelScores_typeH, 2, function(x) (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)) #normalize label scores
labelScores_typeH <- data.frame(labelScores_typeH)
labelScores_typeH_filter <- labelScores_typeH[!is.na(rowSums(labelScores)),]
colnames(labelScores_typeH_filter) <- gsub("ExM\\.U", "ExM-U", colnames(labelScores_typeH_filter)) #fix names

# hierarchy summary stats
tree <- list()
tree$merge <- cbind(
    c(-8, -9, -11, -12, 3,  -5, 5, -14, -3,  -2,  9,  7, -4, -7, -6),
    c(-1,  1,   2, -13, 4, -10, 6, -16,  8, -15, 10, 11, 12, 13, 14))
tree$height <- c(1,2,3,1,4,1,5,1,2,1,3,6,7,8,9)
tree$labels <- colnames(labelScores)[1:16]
tree$order <- c(1,8,9,11,12,13,5,10,14,16,3,2,15,4,7,6)
class(tree) <- "hclust"
labelScores_typeH_filter_c <- labelScores_typeH_filter[,1:16]
for(i in 1:nrow(tree$merge)) { 
    index <- sapply(tree$merge[i,], function(x) ifelse(x<0, abs(x), x+16))
    labelScores_typeH_filter_c <- cbind(labelScores_typeH_filter_c, 
                                                    (labelScores_typeH_filter[,index[1]]+
                                                     labelScores_typeH_filter[,index[2]])+
                                                     labelScores_typeH_filter[,i+16]/2)
}
colnames(labelScores_typeH_filter_c) <- colnames(labelScores_typeH_filter)

# check if child nodes of parent node are significant
eQTL_labelFrame = (labelScores_typeH_filter_c>2)+0
for(i in 1:nrow(tree$merge)) { 
    index <- sapply(tree$merge[i,], function(x) ifelse(x<0, abs(x), x+16))
    sigPoint <- which((eQTL_labelFrame[,index[1]]==1 | eQTL_labelFrame[,index[2]]==1))
    higherNodes <- grepl(paste0(colnames(eQTL_labelFrame)[index[1]],"\\.|",colnames(eQTL_labelFrame)[index[2]],"\\."),colnames(eQTL_labelFrame))
    higherNodes[index] <- FALSE
    eQTL_labelFrame[sigPoint,higherNodes] <- 0
}

# check top scoring cell type (non-hierarchical)
topLabel <- sapply(1:nrow(labelScores_filter), function(i) ifelse(num2LabelScores[i]==0,"Broad",colnames(labelScores_filter)[order(labelScores_filter[i,], decreasing=TRUE)]))

# and get sig hierarchical labels
treeLabel <- sapply(1:nrow(labelScores_typeH_filter_c), function(i) ifelse(sum(eQTL_labelFrame[i,])==0,"Broad",colnames(labelScores_typeH_filter_c)[eQTL_labelFrame[i,]==1][order(labelScores_typeH_filter_c[i,eQTL_labelFrame[i,]==1], decreasing = TRUE)[1]]))

# plot tree with top scores
plotScores <- table(factor(treeLabel, levels=colnames(eQTL_labelFrame)))[nrow(tree$merge)+16]
testList = tree$merge[nrow(tree$merge),]
while(length(testList)>0){
    thisTest <- testList[1]
    testList <- testList[-1]
    if(thisTest<0){
      plotScores <- c(plotScores, table(factor(treeLabel, levels=colnames(eQTL_labelFrame)))[-thisTest])
    } else{
      plotScores <- c(plotScores, table(factor(treeLabel, levels=colnames(eQTL_labelFrame)))[ thisTest+16])
      testList <- c(tree$merge[thisTest,], testList)
    }
}
tree_dend <- as.ggdend(as.dendrogram(tree) %>% hang.dendrogram)
tree_dend$nodes$pch <- 15
tree_dend$nodes$cex <- sqrt(plotScores/10) #linear scale for square area
tree_dend$nodes$col <- rainbow(31)
p <- ggplot(tree_dend)
ggsave(paste0(file, "_tree.pdf"), p, width=5, height=6)

# Give descriptive names
treeLabel_named <- treeLabel
treeLabel_named[!treeLabel_named %in% colnames(labelScores_typeH_filter_c)[1:27]] <- "Broad"
treeLabel_named[treeLabel_named == colnames(labelScores_typeH_filter_c)[23]] <- "Broadly Neuronal"
treeLabel_named <- gsub("\\.", "/", treeLabel_named)

# plot transitions
plotFrame <- table(topLabel, treeLabel_named)
plotFrame_melt <- melt(plotFrame)
plotFrame_melt$topLabel <- factor(plotFrame_melt$topLabel, levels=c('Broad','Per','Mic','End','ExM-U','IP','ExM','ExN','InMGE','InCGE','ExDp2','ExDp1','OPC','oRG','vRG','PgG2M','PgS'))
plotFrame_melt$treeLabel <- factor(plotFrame_melt$treeLabel, levels=c('Broad','Broad1','Per','Broad2','Mic','Broad3','End','Broad4','Broadly Neuronal','ExM/ExN/IP/ExM-U/InMGE/InCGE','ExM/ExN/IP/ExM-U','ExM-U','ExM/ExN/IP','IP','ExM/ExN','ExM','ExN','InMGE/InCGE','InMGE','InCGE','ExDp2/ExDp1','ExDp2','ExDp1','oRG/vRG/OPC/PgG2M/PgS','oRG/vRG/OPC','OPC','oRG/vRG','oRG','vRG','PgG2M/PgS','PgG2M','PgS'))
p <- ggplot(plotFrame_melt[!plotFrame_melt$treeLabel=="Broad",]) + geom_histogram(aes(x=treeLabel, y=value, fill=topLabel), stat="identity") + scale_fill_manual(values=c("gray",rainbow(31)[c(4,6,11,15,16,18,19,21,22,25,30)])) + theme_classic() + theme(axis.text.x = element_text(angle = 90,  hjust=1, vjust=.5)) + labs(fill="Non-hierarchical\nCell Type") + xlab("Hierarchical Cell Type") + ylab("Count")
ggsave(paste0(file, "_newLabels.pdf"), p, width=8, height=6)

# plot umaps with hierarchical and non-hierarchical labels
max_label_tree <- apply(labelScores_typeH_filter_c, 1, function(x) colnames(labelScores_typeH_filter_c)[order(x[1:27], decreasing = TRUE)][1])
max_label_tree_strict <- max_label_tree
max_label_tree_strict[apply(labelScores_typeH_filter_c, 1, max)<2] <- "Broad"
max_label_tree_strict[max_label_tree_strict == colnames(labelScores_typeH_filter_c)[23]] <- "Broadly Neuronal"
max_label_tree_strict <- gsub("\\.", "/", max_label_tree_strict)
label_umap <- uwot::umap(labelScores_filter))  
p <- ggplot() + geom_point(aes(label_umap[,1],label_umap[,2], color=factor(max_label_tree_strict, levels=levels(plotFrame_melt$treeLabel)), alpha=apply(labelScores_typeH_filter_c, 1, max)/3)) + scale_color_manual(values=c("gray",rainbow(31)[which(table(factor(max_label_tree_strict, levels(plotFrame_melt$treeLabel)))>0)-1])) + theme_classic() + labs(color="Hierarchical\nCell Type") + xlab("UMAP-1") + ylab("UMAP-2") + guides(alpha=FALSE, color=FALSE)
ggsave(paste0(file, "_UMAP.png"), p, width=8, height=6)
p <- ggplot() + geom_point(aes(label_umap[,1],label_umap[,2], color=factor(topLabel, levels=levels(plotFrame_melt$topLabel)), alpha=apply(labelScores_filter, 1, max)/3)) + scale_color_manual(values=c("gray",rainbow(31)[c(4,6,11,15,16,18,19,21,22,25,30)])) + theme_classic() + labs(color="Hierarchical\nCell Type") + xlab("UMAP-1") + ylab("UMAP-2") + guides(alpha=FALSE, color=FALSE)
ggsave(paste0(file, "_UMAP_noH.png"), p, width=8, height=6)

# cell-type divergent eQTLs
allBulkQTLs <- ciseQTLs[!is.na(rowSums(labelScores)),]
allBulkQTLs$Type <- max_label_tree_strict
gene_eQTLs_tree <- allBulkQTLs %>% filter(Type!="Broad") %>% 
    group_by(phenotype_id) %>% filter(length(unique(Type))>1)
write.table(gene_eQTLs_tree, paste0(file,"_gene_eQTLs_tree.txt"), row.names=FALSE)

# non-tree scores for comparison
allBulkQTLs$topLabel <- topLabel
gene_eQTLs_basic <- allBulkQTLs %>% filter(topLabel!="Broad") %>% 
    group_by(phenotype_id) %>% filter(length(unique(topLabel))>1)

# Calculate distances in high dimensional and umap space between eQTLs belonging to the same hierarhical cell type
allHTypeDist_umap <- c()
allHTypeDist <- c()
for(thisType in unique(treeLabel)){
    if(thisType=="Broad"){next}
    if(length(which(treeLabel==thisType))<2){next}
    dists <- dist(labelScores_filter[treeLabel==thisType,])
    if(max(dists)==0){next}
    umapDists <- dist(label_umap[treeLabel==thisType,])
    allHTypeDist_umap <- c(allHTypeDist_umap, umapDists[dists!=0])
    allHTypeDist <- c(allHTypeDist, dists[dists!=0])
}
# Calculate distances in high dimensional and umap space between eQTLs belonging to the same eGene
allGeneDist_umap <- c()
allGeneDist <- c()
for(thisGene in unique(allBulkQTLs$phenotype_id)){
    if(length(which(allBulkQTLs$phenotype_id==thisGene))<2){next}
    dists <- dist(labelScores_filter[allBulkQTLs$phenotype_id==thisGene,])
    if(max(dists)==0){next}
    umapDists <- dist(label_umap[allBulkQTLs$phenotype_id==thisGene,])
    allGeneDist_umap <- c(allGeneDist_umap, umapDists[dists!=0])
    allGeneDist <- c(allGeneDist, dists[dists!=0])
}
# Check that eQTLs cluster more by assigned cell type than by eGENE
wilcox.test(allGeneDist_umap, allHTypeDist_umap) # p-value < 2.2e-16
wilcox.test(allGeneDist, allHTypeDist) # p-value < 2.2e-16

# Check distances between top 2
uniqueGenes <- unique(gene_eQTLs_tree$phenotype_id)
maxDist <- sapply(uniqueGenes, function(gene) max(sapply(which(allBulkQTLs$phenotype_id==gene), function(i) sapply(which(allBulkQTLs$phenotype_id==gene), function(j) sum(abs(labelScores_typeH_filter_c[i,]-labelScores_typeH_filter_c[j,]))))))
maxUMAPDist <- sapply(uniqueGenes, function(gene) max(sapply(which(allBulkQTLs$phenotype_id==gene), function(i) sapply(which(allBulkQTLs$phenotype_id==gene), function(j) sum(abs(label_umap[i,]-label_umap[j,]))))))
uniqueGenes_filter <- uniqueGenes[maxDist>40 & maxUMAPDist>8]

# and for non-tree
uniqueGenes_basic <- unique(gene_eQTLs_basic$phenotype_id)
maxDist_basic <- sapply(uniqueGenes_basic, function(gene) max(sapply(which(allBulkQTLs$phenotype_id==gene), function(i) sapply(which(allBulkQTLs$phenotype_id==gene), function(j) sum(abs(labelScores_filter[i,]-labelScores_filter[j,]))))))
maxUMAPDist_basic <- sapply(uniqueGenes_basic, function(gene) max(sapply(which(allBulkQTLs$phenotype_id==gene), function(i) sapply(which(allBulkQTLs$phenotype_id==gene), function(j) sum(abs(label_umap[i,]-label_umap[j,]))))))

# count num div eQTLs
gene_eQTLs_tree_counts <- allBulkQTLs %>% filter(phenotype_id %in% uniqueGenes_filter) %>% filter(Type!="Broad") %>% 
    group_by(phenotype_id) %>% summarize(numTypes = length(unique(Type)))
gene_eQTLs_tree_counts_table <- c("2"=length(which(gene_eQTLs_tree_counts$numTypes==2)), "3"=length(which(gene_eQTLs_tree_counts$numTypes==3)), "4+"=length(which(gene_eQTLs_tree_counts$numTypes>3)))

# and without tree
gene_eQTLs_basic_counts <- allBulkQTLs %>% filter(phenotype_id %in% uniqueGenes_basic) %>% filter(topLabel!="Broad") %>% 
    group_by(phenotype_id) %>% summarize(numTypes = length(unique(topLabel)))
gene_eQTLs_basic_counts_table <- c("2"=length(which(gene_eQTLs_basic_counts$numTypes==2)), "3"=length(which(gene_eQTLs_basic_counts$numTypes==3)), "4+"=length(which(gene_eQTLs_basic_counts$numTypes>3)))

p <- ggplot() + geom_histogram(aes(x="With Hierarchical\nCell Types", y=gene_eQTLs_tree_counts_table, 
                fill=names(gene_eQTLs_tree_counts_table)), stat="identity") + 
           geom_histogram(aes(x="No Hierarchical\nCell Types", y=gene_eQTLs_basic_counts_table, 
                fill=names(gene_eQTLs_basic_counts_table)), stat="identity") + 
           theme_classic() + labs(fill="Num eQTLs\nPer Gene") + ylab("Num Genes") + theme(axis.title.x=element_blank())
ggsave(paste0(file, "_div_eQTL_counts.pdf"), p, width=4, height=6)
# only 88 genes had cell-type divergent eQTLs, now 613 genes do

# upset plot
upsetFrame <- data.table(uniqueGenes_filter, t(sapply(uniqueGenes_filter, function(gene) sapply(unique(gene_eQTLs_tree$Type), function(type) as.numeric(type %in% gene_eQTLs_tree$Type[gene_eQTLs_tree$phenotype_id==gene])))))
fullSet <- c('Broad','Broad1','Per','Broad2','Mic','Broad3','End','Broad4','Broadly Neuronal','ExM/ExN/IP/ExM-U/InMGE/InCGE','ExM/ExN/IP/ExM-U','ExM-U','ExM/ExN/IP','IP','ExM/ExN','ExM','ExN','InMGE/InCGE','InMGE','InCGE','ExDp2/ExDp1','ExDp2','ExDp1','oRG/vRG/OPC/PgG2M/PgS','oRG/vRG/OPC','OPC','oRG/vRG','oRG','vRG','PgG2M/PgS','PgG2M','PgS')
pdf(file=paste0(file, "_div_eQTL_upset.pdf"), 8, 6)
UpSetR::upset(upsetFrame, sets=rev(fullSet[fullSet %in% colnames(upsetFrame)]),  keep.order=TRUE, sets.bar.color=rev(c(rainbow(31)[which(fullSet %in% colnames(upsetFrame))-1])), order.by = "freq", sets.x.label = "Num Genes",  mb.ratio = c(0.6, 0.4))
dev.off()

# topLabel upset plot for comparison
upsetFrame <- data.table(uniqueGenes_basic, t(sapply(uniqueGenes_basic, function(gene) sapply(unique(gene_eQTLs_basic$topLabel), function(type) as.numeric(type %in% gene_eQTLs_basic$topLabel[gene_eQTLs_basic$phenotype_id==gene])))))
pdf(file=paste0(file, "_div_eQTL_topLabel_upset.pdf"), 8, 6)
UpSetR::upset(upsetFrame, sets=rev(fullSet[fullSet %in% colnames(upsetFrame)]),  keep.order=TRUE, sets.bar.color=rev(c(rainbow(31)[which(fullSet %in% colnames(upsetFrame))-1])), order.by = "freq", sets.x.label = "Num Genes")
dev.off()

# check multiome data
ATAC_counts <- fread("/data/GSE162170_multiome_atac_counts.tsv.gz")
ATAC_peaks <- fread("/data/GSE162170_multiome_atac_consensus_peaks.txt.gz")
RNA_counts <- as.data.frame(fread("/data/GSE162170_multiome_rna_counts.tsv.gz"))
cell_data <- fread("/data/GSE162170_multiome_cell_metadata.txt.gz")
cell_cluster_names <- fread("/data/GSE162170_multiome_cluster_names.txt.gz")
gene_eQTLs_tree <- fread("/data/gene_eQTLs_typeH_tree.txt", data.table=FALSE)

# build data frame of peak accessibility to gene expression for each cell
peakExpFrame <- data.frame()
for(gene in unique(gene_eQTLs_tree$phenotype_id)){
    if(!gene %in% RNA_counts$V1){next}
    peaks <- c()
    for(whichLocs in which(gene_eQTLs_tree$phenotype_id==gene)){
        liftLoc <- liftOver(GRanges(paste0("chr",gene_eQTLs_tree$chr[whichLocs]),IRanges(gene_eQTLs_tree$pos[whichLocs],gene_eQTLs_tree$pos[whichLocs])), chain = hg19_hg38_chain)
        peaks <- c(peaks, which(ATAC_peaks$seqnames==as.character(seqnames(liftLoc)) & ATAC_peaks$start<as.numeric(start(liftLoc)) & ATAC_peaks$end>as.numeric(start(liftLoc))))
    }
    peaks <- unique(peaks)
    if(length(peaks)<2){next}
    for(i in 1:length(peaks)){
        exp <- as.numeric(RNA_counts[RNA_counts$V1==gene,which(as.numeric(ATAC_counts[peaks[i],])>0 & colSums(ATAC_counts[peaks[-i],])==0)+1])
        if(length(exp)==0){next}
            peakExpFrame <- rbind(peakExpFrame, data.frame(gene, peak=peaks[i], exp=exp, type=cell_cluster_names$Cluster.Name[match(cell_data$seurat_clusters[as.numeric(ATAC_counts[peaks[i],])>0 & colSums(ATAC_counts[peaks[-i],])==0], cell_cluster_names$Cluster.ID[cell_cluster_names$Assay=="Multiome RNA"])]))
    }
    peakExpFrame <- rbind(peakExpFrame, data.frame(gene, peak=0, exp=as.numeric(RNA_counts[RNA_counts$V1==gene,which(colSums(ATAC_counts[peaks,])==0)+1]), type=cell_cluster_names$Cluster.Name[match(cell_data$seurat_clusters[which(colSums(ATAC_counts[peaks,])==0)+1], cell_cluster_names$Cluster.ID[cell_cluster_names$Assay=="Multiome RNA"])]))
}
write.table(peakExpFrame, "/data/gene_eQTLs_typeH_tree_exp_wType.txt", row.names=FALSE)

# check for sig diff exp when is accessible
peakExpFrame_Gene_pval <- peakExpFrame %>% group_by(gene) %>% summarize(pval=wilcox.test(exp[peak==0], exp[peak!=0])$p.value)
peakExpFrame_Gene_min_pval <- peakExpFrame %>% group_by(gene) %>% summarize(pval=min(sapply(unique(peak), function(p) wilcox.test(exp[peak==0], exp[peak==p])$p.value)))
peakExpFrame_Gene_second_pval <- peakExpFrame %>% group_by(gene) %>% summarize(pval=sort(sapply(unique(peak), function(p) wilcox.test(exp[peak==0], exp[peak==p])$p.value))[2])
p <- ggplot() + geom_histogram(aes(peakExpFrame_Gene_min_pval$pval), bins = 20) + theme_classic() + xlab("p-value") + ylim(c(0,170)) + geom_hline(aes(yintercept=length(which(peakExpFrame_Gene_min_pval$pval>=.5))/10), linetype="dashed")
ggsave(paste0(file, "_top_two_eQTLs_pval1.pdf"), p, width=6, height=4)
p <- ggplot() + geom_histogram(aes(peakExpFrame_Gene_second_pval$pval), bins=20) + theme_classic() + xlab("p-value") + ylim(c(0,170)) + geom_hline(aes(yintercept=length(which(peakExpFrame_Gene_second_pval$pval>=.5))/10), linetype="dashed")
ggsave(paste0(file, "_top_two_eQTLs_pval2.pdf"), p, width=6, height=4)

# calculate all pvals and LFC and generate volcano plots
peakExpFrame_logFC <- peakExpFrame %>% filter(peak!=0) %>% group_by(gene, peak) %>% summarize(logFC=log(mean(exp)/mean(peakExpFrame$exp[peakExpFrame$gene == gene & peakExpFrame$peak==0])))
peakExpFrame_pval <- peakExpFrame %>% filter(peak!=0) %>% group_by(gene, peak) %>% summarize(pval=wilcox.test(exp,peakExpFrame$exp[peakExpFrame$gene == gene & peakExpFrame$peak==0])$p.value)

p <- ggplot() + geom_point(aes(peakExpFrame_logFC$logFC, sapply(-log10(peakExpFrame_pval$pval), function(x) ifelse(x>15,15,x)), color=abs(peakExpFrame_logFC$logFC)>log(2) & peakExpFrame_pval$pval < .05), size=1) + theme_classic() + guides(color=FALSE) + xlab("Log-fold Change in Expression") + ylab("-log10 p-value") + xlim(c(-2,3)) + scale_color_manual(values=c("gray","black")) + geom_vline(aes(xintercept=log(2)), linetype="dashed") + geom_vline(aes(xintercept=-log(2)), linetype="dashed") + geom_hline(aes(yintercept=-log10(.05)), linetype="dashed")
ggsave(paste0(file, "_volcano.pdf"), p, width=4, height=4)

peakExpFrame_pval_LFC <- data.frame(peakExpFrame_pval, logFC=peakExpFrame_logFC$logFC)
peakExpFrame_Gene_min_pval_LFC <- peakExpFrame_pval_LFC %>% group_by(gene) %>% summarize(logFC = logFC[order(pval)][1], pval = pval[order(pval)][1])
peakExpFrame_Gene_second_pval_LFC <- peakExpFrame_pval_LFC %>% group_by(gene) %>% summarize(logFC = logFC[order(pval)][2], pval = pval[order(pval)][2])
p <- ggplot() + geom_point(aes(peakExpFrame_Gene_min_pval_LFC$logFC, sapply(-log10(peakExpFrame_Gene_min_pval_LFC$pval), function(x) ifelse(x>15,15,x)), color=abs(peakExpFrame_Gene_min_pval_LFC$logFC)>log(2) & peakExpFrame_Gene_min_pval_LFC$pval < .05), size=1) + theme_classic() + guides(color=FALSE) + xlab("Log-fold Change in Expression") + ylab("-log10 p-value") + xlim(c(-2,3)) + scale_color_manual(values=c("gray","black")) + geom_vline(aes(xintercept=log(2)), linetype="dashed") + geom_vline(aes(xintercept=-log(2)), linetype="dashed") + geom_hline(aes(yintercept=-log10(.05)), linetype="dashed")
ggsave(paste0(file, "top_volcano.pdf"), p, width=4, height=4)
p <- ggplot() + geom_point(aes(peakExpFrame_Gene_second_pval_LFC$logFC, sapply(-log10(peakExpFrame_Gene_second_pval_LFC$pval), function(x) ifelse(x>15,15,x)), color=abs(peakExpFrame_Gene_second_pval_LFC$logFC)>log(2) & peakExpFrame_Gene_second_pval_LFC$pval < .05), size=1) + theme_classic() + guides(color=FALSE) + xlab("Log-fold Change in Expression") + ylab("-log10 p-value") + xlim(c(-2,3)) + scale_color_manual(values=c("gray","black")) + geom_vline(aes(xintercept=log(2)), linetype="dashed") + geom_vline(aes(xintercept=-log(2)), linetype="dashed") + geom_hline(aes(yintercept=-log10(.05)), linetype="dashed")
ggsave(paste0(file, "second_volcano.pdf"), p, width=4, height=4)

# gene name map:
geneTable <- fread("/data/ENSEMBL_to_SYMBOL.txt", header=FALSE)
peakExpFrame$geneName <- geneTable$V2[match(peakExpFrame$gene,geneTable$V1)]

# get peak to eQTL index
peakTypeMatchIndex <- data.frame()
for(gene in unique(gene_eQTLs_tree$phenotype_id)){
    if(!gene %in% RNA_counts$V1){next}
    for(whichLocs in which(gene_eQTLs_tree$phenotype_id==gene)){
        liftLoc <- liftOver(GRanges(paste0("chr",gene_eQTLs_tree$chr[whichLocs]),IRanges(gene_eQTLs_tree$pos[whichLocs],gene_eQTLs_tree$pos[whichLocs])), chain = hg19_hg38_chain)
        peaks <- which(ATAC_peaks$seqnames==as.character(seqnames(liftLoc)) & ATAC_peaks$start<as.numeric(start(liftLoc)) & ATAC_peaks$end>as.numeric(start(liftLoc)))
        if(length(which(peakTypeFrame$gene==gene & peakTypeFrame$peak %in% peaks))==0) {next}
        peakTypeMatchIndex <- rbind(peakTypeMatchIndex, data.frame(gene_eQTLs_tree[whichLocs,], peaks))
    }
}
write.table(peakTypeMatchIndex, "/data/gene_eQTLs_typeH_tree_type_match_index.txt", row.names=FALSE)

# Plot multi-ome expression levels across cell types when each cell-type divergent eQTL is accessible
plotFrame <- peakExpFrame 
plotFrame$variant_id <- sapply(strsplit(peakTypeMatchIndex$variant_id[match(plotFrame$peak,peakTypeMatchIndex$peaks)], "_"), function(x) paste0(x[1],"_",x[2]))
plotFrame$variant_id <- sapply(strsplit(plotFrame$variant_id, "_"), function(x) paste0(x[1],"_",x[2]))
plotFrame$variant_id[plotFrame$peak==0] <- "No Accessible eQTls"
plotFrame$Type <- peakTypeMatchIndex$Type[match(plotFrame$peak,peakTypeMatchIndex$peaks)]
plotFrame$Type <- gsub("ExM\\.U", "ExM-U", plotFrame$Type)
plotFrame$Type <- gsub("\\.", "/", plotFrame$Type)
plotFrame$Type[plotFrame$peak==0] <- ""
plotFrame$Type[plotFrame$Type=="ExM/ExN/IP/ExM-U/InMGE/InCGE/ExDp2/ExDp1"] <- "Broadly Neuronal"
plotFrame$peakName <- paste0(plotFrame$variant_id, "\n", plotFrame$Type) 
p <- ggplot(plotFrame[plotFrame$gene %in% peakExpFrame_Gene_second_pval$gene[p.adjust(peakExpFrame_Gene_second_pval$pval, method = "fdr")<.05],]) + geom_boxplot(aes(x=peakName, y=log(exp+1))) + facet_wrap(~geneName, scales="free", ncol=4) + stat_summary(aes(x=peakName, y=log(exp+1)), fun.y=mean, geom="point", shape=20, size=3, color="red", fill="red") + theme_classic() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x=element_blank())
ggsave(paste0(file, "Multi_QTL_typeH_exp_filterSig.png"), p, width = 10, height=12)

# Create table of multiome cell type and CellWalkR cell type for each peak and gene
peakTypeFrame <- data.frame()
for(gene in unique(gene_eQTLs_tree$phenotype_id)){
    if(!gene %in% RNA_counts$V1){next}
    peaks <- c()
    for(whichLocs in which(gene_eQTLs_tree$phenotype_id==gene)){
        liftLoc <- liftOver(GRanges(paste0("chr",gene_eQTLs_tree$chr[whichLocs]),IRanges(gene_eQTLs_tree$pos[whichLocs],gene_eQTLs_tree$pos[whichLocs])), chain = hg19_hg38_chain)
        peaks <- c(peaks, which(ATAC_peaks$seqnames==as.character(seqnames(liftLoc)) & ATAC_peaks$start<as.numeric(start(liftLoc)) & ATAC_peaks$end>as.numeric(start(liftLoc))))
    }
    peaks <- unique(peaks)
    if(length(peaks)<2){next}
    for(i in 1:length(peaks)){
        peakTypeFrame <- rbind(peakTypeFrame, data.frame(gene, peak=peaks[i], Type=cell_cluster_names$Cluster.Name[cell_cluster_names$Assay=="Multiome RNA"], TypeCounts=as.numeric(table(factor(cell_data$seurat_clusters[as.numeric(ATAC_counts[peaks[i],])>0 & colSums(ATAC_counts[peaks[-i],])==0], levels=cell_cluster_names$Cluster.ID[cell_cluster_names$Assay=="Multiome RNA"])))))
    }
    peakTypeFrame <- rbind(peakTypeFrame, data.frame(gene, peak=0, Type=cell_cluster_names$Cluster.Name[cell_cluster_names$Assay=="Multiome RNA"], TypeCounts=as.numeric(table(factor(cell_data$seurat_clusters[colSums(ATAC_counts[peaks,])==0], levels=cell_cluster_names$Cluster.ID[cell_cluster_names$Assay=="Multiome RNA"])))))
}
write.table(peakTypeFrame, "/data/gene_eQTLs_typeH_tree_type.txt", row.names=FALSE)

# Calculate enrichment for each cell type
peakTypeFrame <- peakTypeFrame %>% group_by(gene,peak,Type) %>% summarize(typeCounts=sum(TypeCounts))
peakTypeFrame$typeCounts <- peakTypeFrame$typeCounts+1 
peakTypeFrame <- peakTypeFrame %>% group_by(gene, peak) %>% mutate(countFrac=typeCounts/sum(typeCounts))
peakTypeFrame <- peakTypeFrame %>% group_by(gene, peak) %>% mutate(enrichment=countFrac/peakTypeFrame$countFrac[peakTypeFrame$peak==0 & peakTypeFrame$gene==gene])

# Plot enrichment of multi-ome cell types for hierarchical annotations
plotFrame <- peakTypeFrame
plotFrame$geneName <- geneTable$V2[match(plotFrame$gene,geneTable$V1)]
plotFrame$variant_id <- sapply(strsplit(peakTypeMatchIndex$variant_id[match(plotFrame$peak,peakTypeMatchIndex$peaks)], "_"), function(x) paste0(x[1],"_",x[2]))
plotFrame$variant_id <- sapply(strsplit(plotFrame$variant_id, "_"), function(x) paste0(x[1],"_",x[2]))
plotFrame$variant_id[plotFrame$peak==0] <- "No Accessible eQTls"
plotFrame$peakType <- peakTypeMatchIndex$Type[match(plotFrame$peak,peakTypeMatchIndex$peaks)]
plotFrame$peakType <- gsub("ExM\\.U", "ExM-U", plotFrame$peakType)
plotFrame$peakType <- gsub("\\.", "/", plotFrame$peakType)
plotFrame$peakType[plotFrame$peak==0] <- ""
plotFrame$peakType[plotFrame$peakType=="ExM/ExN/IP/ExM-U/InMGE/InCGE/ExDp2/ExDp1"] <- "Broadly Neuronal"
plotFrame$peakName <- paste0(plotFrame$variant_id, "\n", plotFrame$peakType) 
plotFrame <- plotFrame[plotFrame$Type != "EC/Peric.",]
p <- ggplot(plotFrame[plotFrame$gene %in% peakExpFrame_Gene_second_pval$gene[p.adjust(peakExpFrame_Gene_second_pval$pval, method = "fdr")<.05] & plotFrame$peak!=0 & plotFrame$enrichment>1,]) + geom_histogram(aes(x=peakName, y=enrichment, fill=Type), stat="identity", position="dodge") + facet_wrap(~geneName, scales="free_x", ncol=4) + scale_y_log10() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank()) + ylab("Cell-type Enrichment") + labs(fill="Multiome scRNA\nCell Type") 
ggsave(paste0(file, "Multi_QTL_typeH_TypeFracEnrichment.png"), p, width = 10, height=12)

# get matching peak types
peakTypeMatch <- data.frame()
for(gene in unique(gene_eQTLs_tree$phenotype_id)){
    if(!gene %in% RNA_counts$V1){next}
    for(whichLocs in which(gene_eQTLs_tree$phenotype_id==gene)){
        liftLoc <- liftOver(GRanges(paste0("chr",gene_eQTLs_tree$chr[whichLocs]),IRanges(gene_eQTLs_tree$pos[whichLocs],gene_eQTLs_tree$pos[whichLocs])), chain = hg19_hg38_chain)
        peaks <- which(ATAC_peaks$seqnames==as.character(seqnames(liftLoc)) & ATAC_peaks$start<as.numeric(start(liftLoc)) & ATAC_peaks$end>as.numeric(start(liftLoc)))
        #peaks in enichment analysis
        if(length(which(peakTypeFrame$gene==gene & peakTypeFrame$peak %in% peaks))==0) {next}
        enrichmentFrame <- data.frame(t(peakTypeFrame$enrichment[peakTypeFrame$gene==gene & peakTypeFrame$peak %in% peaks]))
        colnames(enrichmentFrame) <- peakTypeFrame$Type[peakTypeFrame$gene==gene & peakTypeFrame$peak %in% peaks]
        peakTypeMatch <- rbind(peakTypeMatch, data.frame(gene_eQTLs_tree[whichLocs,], enrichmentFrame))
    }
}
write.table(peakTypeMatch, "/data/gene_eQTLs_typeH_tree_type_match.txt", row.names=FALSE)

# gene name map:
peakTypeMatch$geneName <- geneTable$V2[match(peakTypeMatch$phenotype_id, geneTable$V1)]

# compute which genes had at least 100 reads
eQTL_exp_Genes <- peakExpFrame_wType %>% group_by(gene) %>% summarise(totalReads = sum(exp)) %>% filter(totalReads>100)

# Look at a few common classes of divergent genes
# select just RG/In divergent genes
upsetFrame <- data.table(uniqueGenes_filter, t(sapply(uniqueGenes_filter, function(gene) sapply(unique(gene_eQTLs_tree$Type), function(type) as.numeric(type %in% gene_eQTLs_tree$Type[gene_eQTLs_tree$phenotype_id==gene])))))
eQTL_RG_In_Div_Genes <- upsetFrame$uniqueGenes_filter[upsetFrame$`oRG/vRG`==1 & upsetFrame$`InMGE/InCGE`==1 & rowSums(upsetFrame[,-1])==2]
p = ggplot(peakTypeMatch %>% filter(peakTypeMatch$phenotype_id %in% eQTL_RG_In_Div_Genes & peakTypeMatch$Type %in% c("oRG.vRG","InMGE.InCGE"))) + geom_point(aes(IN1,RG,color=Type, shape=geneName)) + scale_x_log10(lim=c(.13,5.2)) + scale_y_log10(lim=c(.13,5.2)) + theme_classic() + geom_abline(aes(slope=1,intercept=0), linetype="dashed", alpha=.5) + scale_shape_manual(values=0:8) + labs(color="Hierarchical\neQTL Cell Type", shape="Gene") + xlab("Multimodal Interneuron Enrichment") + ylab("Multimodal RG Enrichment")
ggsave(paste0(file, "_RG_v_IN_Enrichment.pdf"), p, width=6, height=4)
# select just Deeplayer/Upper Maturing regulation
eQTL_Dp_UM_Div_Genes <- upsetFrame$uniqueGenes_filter[upsetFrame$`ExM/ExN/IP`==1 & upsetFrame$`ExDp2/ExDp1`==1 & rowSums(upsetFrame[,-1])==2]
p <- ggplot(peakTypeMatch %>% filter(peakTypeMatch$phenotype_id %in% eQTL_Dp_UM_Div_Genes & peakTypeMatch$Type %in% c("ExM.ExN.IP","ExDp2.ExDp1"))) + geom_point(aes(SP,nIPC.GluN1,color=Type, shape=phenotype_id)) + scale_x_log10(lim=c(.22,17.2)) + scale_y_log10(lim=c(.22,17.2)) + theme_classic() + geom_abline(aes(slope=1,intercept=0), linetype="dashed") + scale_shape_manual(values=0:8) + labs(shape="Gene", color="eQTL Cell Type") + xlab("Subplate Enrichment") + ylab("nIPC/GluN Enrichment")
ggsave(paste0(file, "_Dp_v_UM_Enrichment.pdf"), p, width=6, height=4)
# select just RG/Maturing Neural divergent regulation:
eQTL_RG_M_Div_Genes <- upsetFrame$uniqueGenes_filter[upsetFrame$`oRG/vRG`==1 & upsetFrame$`ExM/ExN/IP/ExM-U`==1]
p <- ggplot(peakTypeMatch %>% filter(peakTypeMatch$phenotype_id %in% eQTL_RG_M_Div_Genes & peakTypeMatch$Type %in% c("oRG.vRG","ExM.ExN.IP.ExM.U"))) + geom_point(aes(GluN4,RG,color=Type, shape=geneName)) + scale_x_log10() + scale_y_log10() + theme_classic() + geom_abline(aes(slope=1,intercept=0), linetype="dashed", alpha=.5) + scale_shape_manual(values=0:8) + labs(color="Hierarchical\neQTL Cell Type", shape="Gene") + xlab("Multimodal Excitatory Neuron Enrichment") + ylab("Multimodal RG Enrichment")
ggsave(paste0(file, "_RG_v_M_Enrichment.pdf"), p, width=6, height=4)

# plot expression for FDR < 0.1 genes when eQTL_1 vs eQTL_2 is accessible
peakExpFrame_wType <- fread("/data/gene_eQTLs_typeH_tree_exp_wType.txt")
cellTypeExp <- peakExpFrame_wType %>% filter(gene %in% peakExpFrame_Gene_second_pval$gene[p.adjust(peakExpFrame_Gene_second_pval$pval, method = "fdr")<.1]) %>% group_by(gene,peak) %>% summarize(exp = mean(exp+1))
cellTypeExp <- cellTypeExp %>% group_by(gene,peak) %>% summarize(FC=exp/cellTypeExp$exp[cellTypeExp$gene==gene & cellTypeExp$peak==0])
cellTypeExp <- cellTypeExp %>% filter(peak!=0)
cellTypeExp$type <- peakTypeMatchIndex$Type[match(cellTypeExp$peak,peakTypeMatchIndex$peaks)]
plotFrame <- cellTypeExp %>% filter(peak!=0) %>% group_by(gene) %>% summarize(eQTL1 = FC[order(abs(log10(FC)), decreasing=TRUE)[1]], eQTL2=FC[order(abs(log10(FC)), decreasing=TRUE)][2], eQTL1_type = type[order(abs(log10(FC)), decreasing=TRUE)[1]], eQTL2_type = type[order(FC, decreasing=TRUE)[2]])
plotFrame <- plotFrame[plotFrame$eQTL1_type!=plotFrame$eQTL2_type,]
plotFrame <- plotFrame[sapply(1:length(plotFrame$eQTL1_type), function(x) !grepl(plotFrame$eQTL1_type[x], plotFrame$eQTL2_type[x])),]
plotFrame <- plotFrame[sapply(1:length(plotFrame$eQTL1_type), function(x) !grepl(plotFrame$eQTL2_type[x], plotFrame$eQTL1_type[x])),]
plotFrame$geneName <- geneTable$V2[match(plotFrame$gene, geneTable$V1)]
p <- ggplot(plotFrame) + geom_text_repel(aes(eQTL1, eQTL2,  label=geneName), min.segment.length = unit(0,"cm"), size=3) + theme_classic() + geom_vline(aes(xintercept=1), linetype="dashed", alpha=.5) + geom_hline(aes(yintercept=1), linetype="dashed", alpha=.5) + xlab("Fold-change in expression when\nonly eQTL-1 is accessible") + ylab("Fold-change in expression when\nonly eQTL-2 is accessible") + scale_x_log10()
ggsave(paste0(file, "_topDiff_exp_GeneLabel.pdf"), p, width=6, height=4)

# also make boxplots of expression in all cells where each eQTL is accessible
peakExpFrame_wType$peakType <- peakTypeMatchIndex$Type[match(peakExpFrame_wType$peak,peakTypeMatchIndex$peaks)]
boxplotFrame <- peakExpFrame_wType %>% filter(paste(gene,peakType) %in% c(paste(plotFrame$gene, plotFrame$eQTL1_type), paste(plotFrame$gene, plotFrame$eQTL2_type)) | (gene %in% plotFrame$gene & peak==0))
p <- ggplot(boxplotFrame) + geom_boxplot(aes(as.character(peak), log(exp+1))) + facet_wrap(~gene, scales="free") + theme_classic()
ggsave(paste0(file, "_topDiff_exp_boxplot.pdf"), p, width=12, height=8)

# plot multiome cell type enrhicment:
plotFrame <- melt(data.frame(peakTypeMatch)[,-c(1:17)])
p <- ggplot(plotFrame) + geom_boxplot(aes(Type, value, fill=Type)) + scale_y_log10() + theme_classic() + facet_wrap(~variable, scale="free_y", ncol=3) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position="bottom") + ylab("Multiome Cell Type Enrichment")
ggsave(paste0(file, "_cell_type_Enrichment.pdf"), p, width=12, height=8)

# Calculate entropy following https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009305#sec018
meanTypeExp <- sapply(unique(cell_data$seurat_clusters), function(x) apply(RNA_counts[,which(cell_data$seurat_clusters==x)+1], 1, function(c) mean(c/cellCounts[which(cell_data$seurat_clusters==x)])))
entropy <- apply(meanTypeExp, 1, function(x) -sum((x+1)/1E6*log((x+1)/1E6)))
geneTypeEntropy <- data.frame(gene=RNA_counts$V1,entropy)
plotFrame <- c()
for(x in 1:4){
    geneSelect <- allBulkQTLs %>% filter(Type!="Broad") %>% group_by(phenotype_id) %>% filter(length(unique(Type))==x)
    plotFrame <- rbind(plotFrame, data.frame(eQTL_Count=x, Entropy=geneTypeEntropy$entropy[geneTypeEntropy$gene %in% geneSelect$phenotype_id]))
}
geneSelect <- allBulkQTLs %>% filter(Type!="Broad") %>% group_by(phenotype_id) %>% filter(length(unique(Type))>=5)
plotFrame <- rbind(plotFrame, data.frame(eQTL_Count=5, Entropy=geneTypeEntropy$entropy[geneTypeEntropy$gene %in% geneSelect$phenotype_id]))
p <- ggplot() + geom_boxplot(aes(as.character(eQTL_Count), Entropy), data=plotFrame) + geom_boxplot(aes(x="0", y=geneTypeEntropy$entropy)) + scale_y_log10() + theme_classic() + xlab("Number of cell-type specific eQTLs")
ggsave(paste0(file, "_cell_type_specific_entropy.png"), p, width=6, height=4)

# Check MPRAs from Deng et al.
# primary cells
cts_ratios <- fread("/data/MPRA/cell_type_specific/ratios.csv")
cts_ratios_ranges <- sapply(strsplit(cts_ratios$insert, "::"), function(x) x[2])
cts_ratios <- cts_ratios[grepl("-", cts_ratios_ranges),]
cts_ratios_type <- sapply(strsplit(cts_ratios$insert, "::"), function(x) x[1])
cts_ratios_ranges <- cts_ratios_ranges[grepl("-", cts_ratios_ranges)]
cts_ratios_ranges <- as(cts_ratios_ranges, "GRanges")
# organoids
organoid_ratios <- fread("/data/MPRA/organoid/ratios.csv")
organoid_ratios_ranges <- sapply(strsplit(organoid_ratios$insert, "::"), function(x) x[2])
organoid_ratios <- organoid_ratios[grepl("-", organoid_ratios_ranges),]
organoid_ratios_type <- sapply(strsplit(organoid_ratios$insert, "::"), function(x) x[1])
organoid_ratios_ranges <- organoid_ratios_ranges[grepl("-", organoid_ratios_ranges)]
organoid_ratios_ranges <- as(organoid_ratios_ranges, "GRanges")
# variants
variant_ratios <- fread("/data/MPRA/variant/ratios-corrected.csv")
variant_ratios_ranges <- sapply(strsplit(variant_ratios$insert, "::"), function(x) x[2])
variant_ratios <- variant_ratios[grepl("-", variant_ratios_ranges),]
variant_ratios_type <- sapply(strsplit(variant_ratios$insert, "::"), function(x) x[1])
variant_ratios_ranges <- variant_ratios_ranges[grepl("-", variant_ratios_ranges)]
variant_ratios_ranges <- as(variant_ratios_ranges, "GRanges")

# get all matches
gene_eQTLs_tree_ranges <- liftOver(VarToGRanges(gene_eQTLs_tree), hg19_hg38_chain)
gene_eQTLs_tree_ratio <- cbind(gene_eQTLs_tree[findOverlaps(gene_eQTLs_tree_ranges, cts_ratios_ranges)@from,],ratio=cts_ratios$ratio[findOverlaps(gene_eQTLs_tree_ranges, cts_ratios_ranges)@to])
gene_eQTLs_tree_ratio_organoid <- cbind(gene_eQTLs_tree[findOverlaps(gene_eQTLs_tree_ranges, organoid_ratios_ranges)@from,],ratio=organoid_ratios$ratio[findOverlaps(gene_eQTLs_tree_ranges, organoid_ratios_ranges)@to])
gene_eQTLs_tree_ratio_variant <- cbind(gene_eQTLs_tree[findOverlaps(gene_eQTLs_tree_ranges, variant_ratios_ranges)@from,],ratio=variant_ratios$ratio[findOverlaps(gene_eQTLs_tree_ranges, variant_ratios_ranges)@to], insert=variant_ratios$insert[findOverlaps(gene_eQTLs_tree_ranges, variant_ratios_ranges)@to])

# Plot MPRAs overlapping cell-type divergent eQTLs for Fabp7
p <- ggplot() + geom_boxplot(aes("chr6_123081681", cts_ratios$ratio[findOverlaps(as("chr6:122760436-122760636", "GRanges"), cts_ratios_ranges)@to])) + geom_boxplot(aes("chr6_123140875", cts_ratios$ratio[findOverlaps(as("chr6:122819630-122819830", "GRanges"), cts_ratios_ranges)@to])) + scale_y_log10() + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank()) + ylab("RNA/DNA Ratio")
ggsave(paste0(file, "_Multi_QTL_typeH_Fabp7_MPRA.png"), p, width=4, height=4)

# Plot MPRAs overlapping cell-type divergent eQTLs for Ical1
p1 <- ggplot() + geom_boxplot(aes(x=paste(sapply(strsplit(variant_id, "_"), function(x) paste0(x[1],"_",x[2])),"\n",Type), y=ratio), data=gene_eQTLs_tree_ratio %>% filter(phenotype_id==geneName)) +  scale_y_log10(limits=c(.5,2)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axiaqss.title.x = element_blank()) + ylab("RNA/DNA Ratio") + geom_hline(aes(yintercept=median(cts_ratios$ratio[cts_ratios_type %in% c("positive_vaccarino_lab","positive_white_lab")])), linetype="dashed") + geom_hline(aes(yintercept=median(cts_ratios$ratio[cts_ratios_type %in% c("negative_white_lab")])), linetype="dashed", color="red")
p2 <- ggplot() + geom_boxplot(aes(x=paste(sapply(strsplit(variant_id, "_"), function(x) paste0(x[1],"_",x[2])),"\n",Type), y=ratio), data=gene_eQTLs_tree_ratio_organoid %>% filter(phenotype_id==geneName)) +  scale_y_log10(limits=c(.5,2)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank()) + ylab("RNA/DNA Ratio") + geom_hline(aes(yintercept=median(organoid_ratios$ratio[organoid_ratios_type %in% c("positive_vaccarino_lab","positive_white_lab")])), linetype="dashed") + geom_hline(aes(yintercept=median(organoid_ratios$ratio[organoid_ratios_type %in% c("negative_white_lab")])), linetype="dashed", color="red")
p3 <- ggplot() + geom_boxplot(aes(x=paste(sapply(strsplit(variant_id, "_"), function(x) paste0(x[1],"_",x[2])),"\n",Type), y=ratio, color=insert), data=gene_eQTLs_tree_ratio_variant %>% filter(phenotype_id==geneName & variant_id==variantID)) +  scale_y_log10(limits=c(.5,2)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank()) + ylab("RNA/DNA Ratio") + guides(color=FALSE)
p <- p2 + p1 + p3 #with patchwork
ggsave(paste0(file, "_Multi_QTL_typeH_ICA1L_MPRA.pdf"), p, width=6, height=4)

# For comparison, count overlaps with Ziffra et al. peaks (https://cells-test.gi.ucsc.edu/cortex-atac/hub/):

# cell-type specific peaks
ziffraPeakFiles <- list.files("/data/Ziffra_Tracks/", pattern = "Specificpeaks", full.names = TRUE)
eQTL_peakOverlaps <- c()
for(bedfile in ziffraPeakFiles){
    cellTypePeaks <- fread(bedfile)
    eQTL_peakOverlaps <- cbind(eQTL_peakOverlaps, countOverlaps(
        liftOver(GRanges(paste0("chr",allBulkQTLs$chr), IRanges(allBulkQTLs$pos, allBulkQTLs$pos)), chain = hg19_hg38_chain),
        GRanges(cellTypePeaks$V1, IRanges(cellTypePeaks$V2, cellTypePeaks$V3))))
}
colnames(eQTL_peakOverlaps) <- sapply(tools::file_path_sans_ext(basename(ziffraPeakFiles)), function(x) {
    segments <- strsplit(x, "_")[[1]]
    paste(segments[-length(segments)], collapse="_")
    })
# Overlap with non-hierarchical labels
plotNames <- names(sort(table(topLabel), decreasing = TRUE))[-1]
plotList <- lapply(1:10, function(x) list(query=UpSetR::elements, params=list("type", plotNames[x:length(plotNames)]), active=TRUE, query.name=plotNames[x], color=rainbow(31)[which(fullSet == plotNames[x])-1]))
pdf(file=paste0(file, "_ZiffraSpecific_topLabel.pdf"), 7, 6)
UpSetR::upset(data.frame(eQTL_peakOverlaps, type=topLabel), sets=colnames(eQTL_peakOverlaps), queries=plotList, order.by = "freq", sets.x.label = "Num Specific Peaks", nintersects = 18)
dev.off()
# Overlap with hierarchical labels
plotNames <- names(sort(table(treeLabel_named), decreasing = TRUE))[-1]
plotList <- lapply(1:15, function(x) list(query=UpSetR::elements, params=list("type", plotNames[x:length(plotNames)]), active=TRUE, query.name=plotNames[x], color=rainbow(31)[which(fullSet == plotNames[x])-1]))
pdf(file=paste0(file, "_ZiffraSpecific_treeLabel.pdf"), 7, 6)
UpSetR::upset(data.frame(eQTL_peakOverlaps, type=treeLabel_named), sets=colnames(eQTL_peakOverlaps), queries=plotList, order.by = "freq", sets.x.label = "Num Specific Peaks", nintersects = 18)
dev.off()
# Count number of divergent eQTLs identified with cell-type specific peaks
gene_eQTLs_specific_counts <- data.frame(allBulkQTLs, Specific=sapply(1:nrow(eQTL_peakOverlaps), function(i) ifelse(sum(eQTL_peakOverlaps[i,])==1, colnames(eQTL_peakOverlaps)[which(eQTL_peakOverlaps[i,]==1)], "Broad"))) %>% filter(Specific!="Broad") %>% group_by(phenotype_id) %>% summarize(numTypes = length(unique(Specific)))
gene_eQTLs_specific_counts_table <- c("2"=length(which(gene_eQTLs_specific_counts$numTypes==2)), "3"=length(which(gene_eQTLs_specific_counts$numTypes==3)), "4+"=length(which(gene_eQTLs_specific_counts$numTypes>3)))

# Enhancers
ziffraPeakFiles <- list.files("/data/DevBrain/Ziffra_Tracks/", pattern = "Enhancer", full.names = TRUE)
eQTL_peakOverlaps_enhancer <- c()
for(bedfile in ziffraPeakFiles){
    cellTypePeaks <- fread(bedfile)
    eQTL_peakOverlaps_enhancer <- cbind(eQTL_peakOverlaps_enhancer, countOverlaps(
        liftOver(GRanges(paste0("chr",allBulkQTLs$chr), IRanges(allBulkQTLs$pos, allBulkQTLs$pos)), chain = hg19_hg38_chain),
        GRanges(cellTypePeaks$V1, IRanges(cellTypePeaks$V2, cellTypePeaks$V3))))
}
colnames(eQTL_peakOverlaps_enhancer) <- tools::file_path_sans_ext(basename(ziffraPeakFiles))
colnames(eQTL_peakOverlaps_enhancer) <-  sapply(tools::file_path_sans_ext(basename(ziffraPeakFiles)), function(x) {
    segments <- strsplit(x, "_")[[1]]
    paste(segments[-length(segments)], collapse="_")
    })
# Overlap with non-hierarchical labels
plotNames <- names(sort(table(topLabel), decreasing = TRUE))[-1]
plotList <- lapply(1:7, function(x) list(query=UpSetR::elements, params=list("type", plotNames[x:length(plotNames)]), active=TRUE, query.name=plotNames[x], color=rainbow(31)[which(fullSet == plotNames[x])-1]))
pdf(file=paste0(file, "_ZiffraEnhancer_topLabel.pdf"), 7, 6)
UpSetR::upset(data.frame(eQTL_peakOverlaps_enhancer, type=topLabel), sets=colnames(eQTL_peakOverlaps_enhancer), queries=plotList, order.by = "freq", sets.x.label = "Num Enhancers", nintersects = 18)
dev.off()
# Overlap with hierarchical labels
plotNames <- names(sort(table(treeLabel_named), decreasing = TRUE))[-1]
plotList <- lapply(1:15, function(x) list(query=UpSetR::elements, params=list("type", plotNames[x:length(plotNames)]), active=TRUE, query.name=plotNames[x], color=rainbow(31)[which(fullSet == plotNames[x])-1]))
pdf(file=paste0(file, "_ZiffraEnhancer_treeLabel.pdf"), 7, 6)
UpSetR::upset(data.frame(eQTL_peakOverlaps_enhancer, type=treeLabel_named), sets=colnames(eQTL_peakOverlaps_enhancer), queries=plotList, order.by = "freq", sets.x.label = "Num Enhancers", nintersects = 18)
dev.off()
# count div eQTls with enhancers
gene_eQTLs_enhancer_counts <- data.frame(allBulkQTLs, Enhancer=sapply(1:nrow(eQTL_peakOverlaps_enhancer), function(i) ifelse(sum(eQTL_peakOverlaps_enhancer[i,])==1, colnames(eQTL_peakOverlaps_enhancer)[which(eQTL_peakOverlaps_enhancer[i,]==1)], "Broad"))) %>% filter(Enhancer!="Broad") %>% group_by(phenotype_id) %>% summarize(numTypes = length(unique(Enhancer)))
gene_eQTLs_enhancer_counts_table <- c("2"=length(which(gene_eQTLs_enhancer_counts$numTypes==2)), "3"=length(which(gene_eQTLs_enhancer_counts$numTypes==3)), "4+"=length(which(gene_eQTLs_enhancer_counts$numTypes>3)))

# combined plot of all counts for div eQTLs
p <- ggplot() + geom_histogram(aes(x="Cell-type Specific\nPeaks", y=gene_eQTLs_specific_counts_table, 
                fill=names(gene_eQTLs_specific_counts_table)), stat="identity") + 
           geom_histogram(aes(x="Cell-type Specific\nEnhancers", y=gene_eQTLs_enhancer_counts_table, 
                fill=names(gene_eQTLs_enhancer_counts_table)), stat="identity") + 
          geom_histogram(aes(x="With Hierarchical\nCell Types", y=gene_eQTLs_tree_counts_table, 
                fill=names(gene_eQTLs_tree_counts_table)), stat="identity") + 
           geom_histogram(aes(x="No Hierarchical\nCell Types", y=gene_eQTLs_basic_counts_table, 
                fill=names(gene_eQTLs_basic_counts_table)), stat="identity") + 
           theme_classic() + labs(fill="Num eQTLs\nPer Gene") + ylab("Num Genes") + theme(axis.title.x=element_blank())
ggsave(paste0(file, "_div_eQTL_specific_enhancer_counts.pdf"), p, width=7, height=6)

# Look at MPRA ML predictions using model from Deng et al.:
MPRA_predict <- fread("/data/MPRA Predictions/predictions-hg19.csv")
MPRA_predict$variant_id <- gsub(":","_",MPRA_predict$variant_id)
gene_eQTLs_tree_ratio_predict <- merge(x = gene_eQTLs_tree, y = MPRA_predict, by = "variant_id", all.x = TRUE)
mergeFrame <- peakExpFrame_logFC
mergeFrame$variant_id <- peakTypeMatchIndex$variant_id[match(mergeFrame$peak,peakTypeMatchIndex$peaks)]
mergeFrame$phenotype_id <- mergeFrame$gene
gene_eQTLs_tree_ratio_predict_exp <- merge(x = gene_eQTLs_tree_ratio_predict, y = mergeFrame, by = c("phenotype_id", "variant_id"), all.x = TRUE)
# see how often multiple variants are predicted to be functional
plotFrame <- gene_eQTLs_tree_ratio_predict_exp %>% filter(peak!=0) %>% group_by(gene) %>% summarize(eQTL1 = alt_ref[order(abs(log10(alt_ref)), decreasing=TRUE)[1]], eQTL2=alt_ref[order(abs(log10(alt_ref)), decreasing=TRUE)][2], eQTL1_type = Type[order(abs(log10(alt_ref)), decreasing=TRUE)[1]], eQTL2_type = Type[order(alt_ref, decreasing=TRUE)[2]])
# how many have at leat one signficant eQTL
length(which(abs(log(plotFrame$eQTL1))>.2)) #about 23%
plotFrame <- plotFrame[plotFrame$eQTL1_type!=plotFrame$eQTL2_type,]
plotFrame <- plotFrame[sapply(1:length(plotFrame$eQTL1_type), function(x) !grepl(plotFrame$eQTL1_type[x], plotFrame$eQTL2_type[x])),]
plotFrame <- plotFrame[sapply(1:length(plotFrame$eQTL1_type), function(x) !grepl(plotFrame$eQTL2_type[x], plotFrame$eQTL1_type[x])),]
# how many have a second signficant eQTL in a different cell type
length(which(abs(log(plotFrame$eQTL2))>.2)) #about 16%
plotFrame$geneName <- geneTable$V2[match(plotFrame$gene, geneTable$V1)]
p <- ggplot() + geom_point(aes(plotFrame$eQTL1, plotFrame$eQTL2)) + geom_text_repel(aes(plotFrame$eQTL1[abs(log(plotFrame$eQTL1)) + abs(log(plotFrame$eQTL2)) > .5], plotFrame$eQTL2[abs(log(plotFrame$eQTL1)) + abs(log(plotFrame$eQTL2)) > .5],  label=plotFrame$geneName[abs(log(plotFrame$eQTL1)) + abs(log(plotFrame$eQTL2)) > .5]), min.segment.length = unit(0,"cm"), size=3) + theme_classic() + geom_vline(aes(xintercept=1), linetype="dashed", alpha=.5) + geom_hline(aes(yintercept=1), linetype="dashed", alpha=.5) + xlab("Predicted Alt/Ref RNA/DNA Ratio w/ eQTL-1 Variant") + ylab("Predicted Alt/Ref RNA/DNA\nRatio w/ eQTL-2 Variant") + scale_x_log10()
ggsave(paste0(file, "_Multi_QTL_typeH_predicted_MPRA.pdf"), p, width=6, height=3)

# Look at TFs w/ broken motifs
load("/data/sc_dev_cortex_geschwind/raw_counts_mat.rdata") #raw_counts_mat
norm_counts_mat <- raw_counts_mat %*% Diagonal(x = 1 / colSums(raw_counts_mat))
norm_counts_mat@x <- log(10000*norm_counts_mat@x)
metaData <- fread("/data/sc_dev_cortex_geschwind/cell_metadata.csv")
types <- c('oRG/vRG','ExDp2/ExDp1','PgG2M/PgS','ExM/ExN/IP','Mic','InMGE/InCGE','ExM/ExN/IP/ExM-U','End','ExM/ExN','oRG/vRG/OPC','Broadly Neuronal','oRG/vRG/OPC/PgG2M/PgS','Per','ExM/ExN/IP/ExM-U/InMGE/InCGE')
meanExp <- sapply(types, function(type) {
        thisType = ifelse(type=="Broadly Neuronal", "ExM/ExN/IP/ExM-U/InMGE/InCGE/ExDp2/ExDp1", type)
        rowMeans(norm_counts_mat[,metaData$Cluster %in% strsplit(thisType, "/| ")[[1]]])
    })
motifbreakr <- fread("/data/MPRA Predictions/motifbreakr-hg19.tsv")
motifbreakr$variant_id <- gsub(":","_",motifbreakr$SNP_id)
motifbreakr_var <- merge(motifbreakr, gene_eQTLs_tree, by="variant_id")
motifbreakr_var$exp <- sapply(1:nrow(motifbreakr_var),function(x) ifelse(motifbreakr_var$geneSymbol[x] %in% meanExp$V1, meanExp[[motifbreakr_var$Type[x]]][meanExp$V1==motifbreakr_var$geneSymbol[x]], NA))

# create gene x type x TF table 
gene_eQTLs_GeneTypeTF <- motifbreakr_var %>% filter(phenotype_id %in% uniqueGenes_filter, exp>.01) %>% 
    group_by(phenotype_id) %>% filter(length(unique(Type))>1) %>% summarize(geneType = paste0(phenotype_id, "_", Type), TF=geneSymbol, Type)
gene_eQTLs_GeneTypeTF <- unique(gene_eQTLs_GeneTypeTF)

# number of TFs by Type
gene_eQTLs_GeneTypeTF_Count <- gene_eQTLs_GeneTypeTF %>% group_by(TF,Type) %>% summarize(Count = n()) 
gene_eQTLs_GeneTypeTF_CountFilter <- gene_eQTLs_GeneTypeTF_Count[gene_eQTLs_GeneTypeTF_Count$Count>2,]
gene_eQTLs_GeneTypeTF_Count_mat <- igraph::as_adjacency_matrix(igraph::graph.data.frame(gene_eQTLs_GeneTypeTF_CountFilter, directed = FALSE), attr="Count")
pdf(file=paste0(file, "_cell_type_TFs.pdf"), 12, 8)
heatmap.2(as.matrix(gene_eQTLs_GeneTypeTF_Count_mat)[-(1:length(unique(gene_eQTLs_GeneTypeTF_CountFilter$TF))),1:length(unique(gene_eQTLs_GeneTypeTF_CountFilter$TF))][c(6,5,9,1,3,2,4,10,7,8,11),], scale = "none", Rowv=NA, dendrogram = "column", revC=TRUE, col=colorpanel(18,low="white",high="red"), trace = "none", density.info = "none", sepwidth=c(0.01,0.01),sepcolor="gray", colsep=1:ncol(gene_eQTLs_GeneTypeTF_Count_mat),rowsep=1:nrow(gene_eQTLs_GeneTypeTF_Count_mat))
dev.off()

# find all most common pairs
allPairs <- c()
for(gene in unique(gene_eQTLs_GeneTypeTF$phenotype_id)){
    theseTFs <- gene_eQTLs_GeneTypeTF[gene_eQTLs_GeneTypeTF$phenotype_id==gene,]
    for(TF1 in 1:(nrow(theseTFs)-1)){
        for(TF2 in (TF1+1):nrow(theseTFs)){
            if(theseTFs$TF[TF1]!=theseTFs$TF[TF2] & theseTFs$Type[TF1]!=theseTFs$Type[TF2]){
                allPairs <- rbind(allPairs, data.frame(TF1=theseTFs$TF[TF1],Type1=theseTFs$Type[TF1],TF2=theseTFs$TF[TF2],Type2=theseTFs$Type[TF2]))
            }
        }
    }
}
allPairsCount <- allPairs %>% group_by(TF1,TF2) %>% summarize(Count = n()) 
allPairsCount[order(allPairsCount$Count, decreasing = TRUE)[1:10],]
allPairsCount <- allPairs %>% group_by(TF1,TF2, Type1, Type2) %>% summarize(Count = n()) 
allPairsCount[order(allPairsCount$Count, decreasing = TRUE)[1:10],]
allPairsCountFilter <- allPairsCount[allPairsCount$Count>2,]
allPairsCount_mat <- igraph::as_adjacency_matrix(igraph::graph.data.frame(data.frame(paste(allPairsCountFilter$TF1, allPairsCountFilter$Type1),paste(allPairsCountFilter$TF2, allPairsCountFilter$Type2),Count=allPairsCountFilter$Count), directed = FALSE), attr="Count")
pdf(file=paste0(file, "_cell_type_div_TFs.pdf"), 12, 12)
heatmap.2(as.matrix(allPairsCount_mat), scale = "none", symm = TRUE, revC=TRUE, col=colorpanel(15,low="white",high="red"), trace = "none", density.info = "none")
dev.off()
