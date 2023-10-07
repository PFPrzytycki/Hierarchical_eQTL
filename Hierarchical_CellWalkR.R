## This script appends a tree to CellWalkRs Random Walk

#joins in Polioudakis tree
cellTypesH = cbind(
    c(-8, -9, -11, -12, 3,  -5, 5, -14, -3,  -2,  9,  7, -4, -7, -6),
    c(-1,  1,   2, -13, 4, -10, 6, -16,  8, -15, 10, 11, 12, 13, 14))

#function to turn tree into matrix
treeMatrix = function(treeMerge) {
    A = Matrix(0, 2*nrow(treeMerge)+1, 2*nrow(treeMerge)+1)
    for(i in 1:nrow(treeMerge)){ #should build names here too
        children = treeMerge[i,]
        children = sapply(children, function(x) ifelse(x<0, abs(x), x+nrow(treeMerge)+1)) 
        A[children, i+nrow(treeMerge)+1] = 1
        A[i+nrow(treeMerge)+1, children] = 1
    }
    A
}
cellTypesH = treeMatrix(cellTypesH)
colnames(cellTypesH) = rownames(cellTypesH) = allCellTypes

#add tree to label edges 
l = length(cellTypes) #vector of cell type names
l_all = length(allCellTypes) #vector with cell type names and internal node names
expandLabelEdges = cbind(labelEdges, matrix(0,dim(labelEdges)[1],l_all-l)) #add 0s to labelEdges to allow room for internal nodes 
weight = 100 #a weight tuned withouth hierarchy
combinedGraph = rbind(cbind(cellTypesH,t(weight*expandLabelEdges)), #combined graph w/ internal nodes
        cbind(weight*expandLabelEdges,cellEdges))

#compute CellWalk on expanded graph using CellWalkR functions
infMat <- randomWalk(combinedGraph, steps=5) #5 steps is usually enough
normMat <- normalizeInfluence(infMat[-(1:l_all),(1:l_all)]) #normMat and cellLabels are part of usual cellWalk object
colnames(normMat) <- allCellTypes
cellLabels <- apply(normMat, 1, function(x) colnames(normMat)[order(x, decreasing = TRUE)][1]) #top scoring label
cellWalkH <- list(infMat=infMat, normMat=normMat, cellLabels=cellLabels) #make cellWalk object
class(cellWalkH) <- "cellWalk"